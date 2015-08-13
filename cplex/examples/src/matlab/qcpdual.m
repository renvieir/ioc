function qcpdual()
% Illustrates how to query and analyze dual values of a QCP.
%
% ---------------------------------------------------------------------------
% File: qcpdual.m
% Version 12.6.1
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
% Copyright IBM Corporation 2008, 2014. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
% ---------------------------------------------------------------------------

   % Calculate dual multipliers for quadratic constraints.
   % CPLEX does not give us the dual multipliers for quadratic
   % constraints directly. This is because they may not be properly
   % defined at the cone top and deciding whether we are at the cone
   % top or not involves (problem specific) tolerance issues. CPLEX
   % instead gives us all the values we need in order to compute the
   % dual multipliers if we are not at the cone top.
   % The dual multiplier is the dual slack divided
   % by the derivative evaluated at the optimal solution. If the optimal
   % solution is 0 then the derivative at this point is not defined (we are
   % at the cone top) and we cannot compute a dual multiplier.
   function qpi = getqconstrmultipliers(cplex,tol)
      qdslack = cplex.Solution.qcpdslack;
      if max(abs(cplex.Solution.x)) <= tol
         throw(MException('cplex:conetop', 'Cannot compute duals at cone top'));
      end
      qpi = zeros(length(cplex.Model.qc), 1);
      for q = 1:length(qpi)
         thisslack = qdslack(:,q);
         grad = cplex.Model.qc(q).Q' * cplex.Solution.x + cplex.Model.qc(q).Q * cplex.Solution.x + cplex.Model.qc(q).a;

         [maxval,maxind] = max(abs(grad));
         if maxval <= tol
            qpi(q) = 0;
         else
            qpi(q) = qdslack(maxind,q) / grad(maxind);
         end
      end
      getqconstrmultipliers = qpi;
   end

try
   % Initialize the CPLEX object
   cplex = Cplex('qcpdual');
   
   % Create the model
   createmodel(cplex);
   
   % Optimize the problem
   cplex.Param.barrier.qcpconvergetol.Cur = 1e-10;
   cplex.solve();
   
   % Check KKT conditions and print result.
   result = checkkkt(cplex, 1e-6);
   if result == 0
      fprintf('Testing of KKT conditions failed.\n');
   elseif result == -1
      fprintf('Could not test KKT conditions at cone top.\n');
   else
      fprintf('KKT conditions are satisfied.\n');
   end
catch m
   throw (m);
   clear all;
   clear classes;
end

   % The function creates the following model:
   %   Minimize
   %      obj: 3x1 - x2 + 3x3 + 2x4 + x5 + 2x6 + 4x7
   %   Subject To
   %      c1: x1 + x2 = 4
   %      c2: x1 + x3 >= 3
   %      c3: x6 + x7 <= 5
   %      c4: -x1 + x7 >= -2
   %      q1: [ -x1^2 + x2^2 ] <= 0
   %      q2: [ 4.25x3^2 -2x3*x4 + 4.25x4^2 - 2x4*x5 + 4x5^2  ] + 2 x1 <= 9.0
   %      q3: [ x6^2 - x7^2 ] >= 4
   %   Bounds
   %    0 <= x1 <= 3
   %    x2 Free
   %    0 <= x3 <= 0.5
   %    x4 Free
   %    x5 Free
   %    x7 Free
   %   End
   function createmodel(cplex)
      cplex.Model.sense = 'minimize';
      cplex.addCols([3 -1 3 2 1 2 4]', [], ...
                    [0; -inf; 0; -inf; -inf; 0; -inf], ...
                    [3; inf; 0.5; inf; inf; inf; inf]);
      cplex.addRows([4, 3, -inf, -2]', ...
                    [1  1  0  0  0  0  0;
                     1  0  1  0  0  0  0;
                     0  0  0  0  0  1  1;
                     -1 0  0  0  0  0  1], ...
                    [4, inf, 5, inf]', ...
                    ['c1'; 'c2'; 'c3'; 'c4']);
      cplex.addQCs([0 0 0 0 0 0 0;
                    2 0 0 0 0 0 0;
                    0 0 0 0 0 0 0 ]', ...
                    { [-1  0  0  0  0  0  0;
                        0  1  0  0  0  0  0;
                        0  0  0  0  0  0  0;
                        0  0  0  0  0  0  0;
                        0  0  0  0  0  0  0;
                        0  0  0  0  0  0  0;
			0  0  0  0  0  0  0],...
                      [ 0  0  0  0  0  0  0;
                        0  0  0  0  0  0  0;
                        0  0  4.25 -2  0  0  0;
                        0  0  0  4.25 -2  0  0;
                        0  0  0  0  4  0  0;
                        0  0  0  0  0  0  0;
			0  0  0  0  0  0  0],...
                      [ 0  0  0  0  0  0  0;
                        0  0  0  0  0  0  0;
                        0  0  0  0  0  0  0;
                        0  0  0  0  0  0  0;
                        0  0  0  0  0  0  0;
                        0  0  0  0  0  1  0;
			0  0  0  0  0  0 -1] }, ...
                    'LLG', [0.0, 9.0, 4.0], {'q1', 'q2', 'q3'});
   end

   % The function returns 1 if the tested KKT conditions are satisfied and
   % false otherwise.
   % It returns -1 if dual multipliers for quadratic constraints could not
   % be calculated because the optimal solution is at the cone top.
   function satisfied = checkkkt(cplex, tol)
      satisfied = 0;
      % Query solution.
      x = cplex.Solution.x;
      slack = zeros(length(cplex.Model.lhs), 1);
      for i = 1:length(slack)
         if cplex.Model.rhs(i) < inf
            slack(i) = cplex.Model.rhs(i) - cplex.Solution.ax(i);
         else
            slack(i) = cplex.Solution.ax(i) - cplex.Model.lhs(i);
         end
      end
      qslack = zeros(length(cplex.Model.qc), 1);
      for i = 1:length(cplex.Model.qc)
         qslack(i) = cplex.Model.qc(i).rhs - x' * cplex.Model.qc(i).Q * x - cplex.Model.qc(i).a' * x;
      end
      cpi = cplex.Solution.reducedcost;
      rpi = cplex.Solution.dual;
      qpi = getqconstrmultipliers(cplex, tol);

      % Test primal feasibility.
      % This example illustrates the use of dual vectors returned by CPLEX
      % to verify dual feasibility, so we do not test primal feasibility
      % here.

      % Test dual feasibility.
      % We must have
      % - for all <= constraints the dual multiplier is non-positive,
      % - for all >= constraints the dual multiplier is non-negative,
      for i = 1:length(rpi)
         if cplex.Model.sense(i) == 'L' && rpi(i) > tol
            fprintf('Invalid dual multiplier %f for linear row %d\n', rpi(i), i);
            return;
         elseif cplex.Model.sense(i) == 'G' && rpi(i) < -tol
            fprintf('Invalid dual multiplier %f for linear row %d\n', rpi(i), i);
            return;
         end
      end
      for q = 1:length(qslack)
         if cplex.Model.qc(q).sense == 'L' && qpi(q) > tol
            fprintf('Invalid dual multiplier %f for quadratic constraint %d\n', qpi(q), q);
            return;
         elsif cplex.Model.qc(q).sense == 'G' && qpi(q) < -tol
            fprintf('Invalid dual multiplier %f for quadratic constraint %d\n', qpi(q), q);
            return;
         end
      end

      % Test complementary slackness.
      % For any constraint the product of primal slack and dual multiplier
      % must be 0.
      [maxval,maxind] = max(abs(rpi .* slack));
      if maxval > tol
         fprintf('Complementary slackness not satisfied for row %d (%f, %f)\n', maxind, slack(maxind), rpi(maxind));
         return;
      end
      [maxval,maxind] = max(abs(qpi .* qslack));
      if maxval > tol
         fprintf('Complementary slackness not satisfied for quadratic constraint %d (%f, %f)\n', maxind, qslack(maxind), qpi(maxind));
         return;
      end
      for j = 1:length(x)
         if cplex.Model.ub(j)
            slk = cplex.Model.ub(j) - x(j);
            if cpi(j) < -tol
               dual = cpi(j);
            else
               dual = 0.0;
            end
            if abs(slk * dual) > tol
               fprintf('Complementary slackness not satisfied for ub on %d (%f, %f)\n', j, slk, dual);
               return;
            end
         end
         if cplex.Model.lb(j) > -inf
            slk = x(j) - cplex.Model.lb(j);
            if cpi(j) > tol
               dual = cpi(j);
            else
               dual = 0.0;
            end
            if abs(slk * dual) > tol
               fprintf('Complementary slackness not satisfied for lb on %d (%f, %f)\n', j, slk, dual);
               return;
            end
         end
      end

      % Test stationarity.
      % The difference between objective function and gradient at optimal
      % solution multiplied by dual multipliers must be 0, i.e., for the
      % optimal solution x
      % 0 == c
      %      - sum(r in rows)  r'(x)*rpi[r]
      %      - sum(q in quads) q'(x)*qpi[q]
      %      - sum(c in cols)  b'(x)*cpi[c]
      % where r' and q' are the derivatives of a row or quadratic constraint,
      % x is the optimal solution and rpi[r] and qpi[q] are the dual
      % multipliers for row r and quadratic constraint q.
      % b' is the derivative of a bound constraint and cpi[c] the dual bound
      % multiplier for column c.

      % Objective function.
      kktsum = cplex.Model.obj;

      % Linear constraints.
      % The derivative of a linear constraint ax - b (<)= 0 is just a.
      kktsum = kktsum - cplex.Model.A' * rpi;

      % Quadratic constraints.
      % The derivative of a constraint xQx + ax - b <= 0 is
      % Qx + Q'x + a.
      for q = 1:length(cplex.Model.qc)
         kktsum = kktsum - qpi(q) * (cplex.Model.qc(q).Q' * x + cplex.Model.qc(q).Q * x + cplex.Model.qc(q).a);
      end

      % Bound constraints.
      % The derivative for lower bounds is -1 and that for upper bounds
      % is 1.
      % CPLEX already returns dj with the appropriate sign so there is
      % no need to distinguish between different bound types here.
      kktsum = kktsum - cpi;

      % Test that all entries in kktsum are 0.
      [maxval,maxind] = max(abs(kktsum));
      if maxval > tol
         fprintf('Stationarity test failed at index %d: %f\n', maxind, kktsum(maxind));
         return;
      end

      % Print the solution that satisfies KKT conditions.
      fprintf('Optimal solution satisfies KKT conditions.\n');
      x
      cpi
      rpi
      qpi

      satisfied = 1;
   end
end
