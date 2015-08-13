function socpex1()
% Optimize a second order cone programming problem and test KKT conditions
%
% socpex1.m -- Solve a second order cone program to optimality, fetch
%              the dual values and test that the primal and dual solution
%              vectors returned by CPLEX satisfy the KKT conditions.
%              The problems that this code can handle are second order
%              cone programs in standard form. A second order cone
%              program in standard form is a problem of the following
%              type (c' is the transpose of c):
%                min c1'x1 + ... + cr'xr
%                  A1 x1 + ... + Ar xr = b
%                  xi in second order cone (SOC)
%              where xi is a vector of length ni. Note that the
%              different xi are orthogonal. The constraint "xi in SOC"
%              for xi=(xi[1], ..., xi[ni]) is
%                  xi[1] >= |(xi[2],...,xi[ni])|
%              where |.| denotes the Euclidean norm. In CPLEX such a
%              constraint is formulated as
%                  -xi[1]^2 + xi[2]^2 + ... + xi[ni]^2 <= 0
%                   xi[1]                              >= 0
%                             xi[2], ..., xi[ni] free
%              Note that if ni = 1 then the second order cone constraint
%              reduces to xi[1] >= 0.

% ---------------------------------------------------------------------------
% File: socpex1.m
% Version 12.6.1
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
% Copyright IBM Corporation 2008, 2014. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
% ---------------------------------------------------------------------------

   % Compute dual multipliers for second order cone constraints.
   % Note that for a pure second order cone program the dual
   % multiplier for a second order cone constraint is just the dual slack
   % slack value of that constraint's cone head variable.
   function socppi = getsocpconstrmultipliers(cplex)
      % Compute full dual slack.
      dslack = cplex.Solution.reducedcost;
      for i = 1:length(cplex.Model.qc)
         dslack = dslack + cplex.Solution.qcpdslack(:,i);
      end

      % Find cone head variables and pick up dual slacks for them
      socppi = zeros(length(cplex.Model.qc));
      for i = 1:length(cplex.Model.qc)
         % In a second order cone constraints all coefficients are +/-1.
         % The cone head variable has coefficient -1 so we can easily
         % find it by finding the element with smallest coeffcient.
         [minval,minind] = min(diag(cplex.Model.qc(i).Q));
         socppi(i) = dslack(minind);
      end

      getsocpconstrmultipliers = socppi;
   end

try
   % Initialize the CPLEX object
   cplex = Cplex('socpex1');
   
   % Create the model
   constat = createmodel(cplex);
   
   % Optimize the problem
   cplex.Param.barrier.qcpconvergetol.Cur = 1e-9;
   cplex.solve();
   
   % Check the KKT conditions
   if checkkkt(cplex, constat, 1e-9) == 0
      fprintf('Testing of KKT conditions failed.\n');
   else
      fprintf('KKT conditions are satisfied.\n');
   end
catch m
   throw (m);
   clear all;
   clear classes;
end

   % This function creates the following model:
   %   Minimize
   %    obj: x1 + x2 + x3 + x4 + x5 + x6
   %   Subject To
   %    c1: x1 + x2      + x5      = 8
   %    c2:           x3 + x5 + x6 = 10
   %    q1: [ -x1^2 + x2^2 + x3^2 ] <= 0
   %    q2: [ -x4^2 + x5^2 ] <= 0
   %   Bounds
   %    x2 Free
   %    x3 Free
   %    x5 Free
   %   End
   % which is a second order cone program in standard form.
   % The function returns a constat array that is setup as follows:
   % constat >= 0   Column j is contained in a cone constraint
   %                and is the cone head variable of that
   %                constraint. The index of the respective
   %                quadratic constraint is given by constat[j].
   % constat = -1   Column j is contained in a cone constraint
   %                but is not the cone head variable of that
   %                constraint.
   % constat = -2   Column j is not contained in any cone
   %                constraint.
   function constat = createmodel(cplex)
      cplex.Model.sense = 'minimize';
      cplex.addCols([1 1 1 1 1 1]', [], [0; -inf; -inf; 0; -inf; 0], [inf; inf; inf; inf; inf; inf]);
      cplex.addRows(8,  [1 1 0 0 1 0], 8, 'c1');
      cplex.addRows(10, [0 0 1 0 1 1], 10, 'c1');
      cplex.addQCs([0 0 0 0 0 0]', [-1  0  0  0  0  0;
                                     0  1  0  0  0  0;
                                     0  0  1  0  0  0;
                                     0  0  0  0  0  0;
                                     0  0  0  0  0  0;
                                     0  0  0  0  0  0], 'L', 0.0, {'q1'});
      cplex.addQCs([0 0 0 0 0 0]', [  0  0  0  0  0  0;
                                      0  0  0  0  0  0;
                                      0  0  0  0  0  0;
                                      0  0  0 -1  0  0;
                                      0  0  0  0  1  0;
                                      0  0  0  0  0  0], 'L', 0.0, {'q2'});
      constat = [ 1; -1; -1; 2; -1; -2 ];
   end

   % The function returns 1 if the tested KKT conditions are satisfied and
   % false otherwise.
   function satisfied = checkkkt(cplex, constat, tol)
      satisfied = 0;
      % Read primal and dual values.
      x = cplex.Solution.x;
      linslack = cplex.Model.rhs - cplex.Solution.ax;
      socpslack = zeros(length(cplex.Model.qc), 1);
      for i = 1:length(cplex.Model.qc)
         socpslack(i) = x' * cplex.Model.qc(i).Q * x;
      end
      dslack = cplex.Solution.reducedcost;
      for i = 1:length(cplex.Model.qc)
         dslack = dslack + cplex.Solution.qcpdslack(:,i);
      end
      linpi = cplex.Solution.dual;
      socppi = getsocpconstrmultipliers(cplex);

      % Print out the values we just fetched.
      fprintf('x         = [');
      for i = 1:length(x)
         fprintf(' %+7.3f', x(i));
      end
      fprintf(' ]\n');
      fprintf('dslack    = [');
      for i = 1:length(dslack)
         fprintf(' %+7.3f', dslack(i));
      end
      fprintf(' ]\n');
      fprintf('linpi     = [');
      for i = 1:length(linpi)
         fprintf(' %+7.3f', linpi(i));
      end
      fprintf(' ]\n');
      fprintf('linslack  = [');
      for i = 1:length(linslack)
         fprintf(' %+7.3f', linslack(i));
      end
      fprintf(' ]\n');
      fprintf('socppi    = [');
      for i = 1:length(socppi)
         fprintf(' %+7.3f', socppi(i));
      end
      fprintf(' ]\n');
      fprintf('socpslack = [');
      for i = 1:length(socpslack)
         fprintf(' %+7.3f', socpslack(i));
      end
      fprintf(' ]\n');

      % Test primal feasibility.
      % This example illustrates the use of dual vectors returned by CPLEX
      % to verify dual feasibility, so we do not test primal feasibility
      % here.

      % Test dual feasibility.
      % We must have
      % - for all <= constraints the respective pi value is non-negative,
      % - for all >= constraints the respective pi value is non-positive,
      % - since all quadratic constraints are <= constraints the socppi
      %   value must be non-negative for all quadratic constraints,
      % - the dslack value for all non-cone variables must be non-negative.
      for i = 1:length(linpi)
         if cplex.Model.sense(i) == 'L' && linpi(i) < -tol
            fprintf('Invalid dual multiplier %f for linear row %d\n', linpi(i), i);
            return;
         elseif cplex.Model.sense(i) == 'G' && linpi(i) > tol
            fprintf('Invalid dual multiplier %f for linear row %d\n', linpi(i), i);
            return;
         end
      end
      for i = 1:length(socppi)
        if socppi(i) < -tol
           fprintf('Invalid dual multipler %f for quadratic constraint %d\n', socppi(i), i);
           return;
        end
      end
      for i = 1:length(constat)
         if constat(i) == -2 && dslack(i) < -tol
            fprintf('dslack value for column %d is invalid: %f\n', i, dslack(i));
            return;
         end
      end

      % Test complementary slackness.
      % For each constraint either the constraint must have zero slack or
      % the dual multiplier for the constraint must be 0. For a variable that
      % is not in a cone constraint the dual multiplier is in dslack.
      for i = 1:length(linpi)
         if abs(linpi(i)) > tol && abs(linslack(i)) > tol
            fprintf('Complementary slackness not satisfied for row %d (%f, %f)\n', i, linslack(i), linpi(i));
            return;
         end
      end
      for i = 1:length(socppi)
         if abs(socppi(i)) > tol && abs(socpslack(i)) > tol
            fprintf('Complementary slackness not satisfied for quadratic constraint %d (%f, %f)\n', i, socpslack(i), socppi(i));
            return;
         end
      end
      for i = 1:length(constat)
         if constat(i) == -2 && abs(x(i)) > tol && abs(dslack(i)) > tol
            fprintf('Complementary slackness not satisfied for non-cone variable %d: (%f, %f)\n', i, x(i), dslack(i));
            return;
         end
      end

      % Test stationarity.
      % We must have
      %  c - g[i]'(X)*pi[i] = 0
      % where c is the objective function, g[i] is the i-th constraint of the
      % problem, g[i]'(x) is the derivate of g[i] with respect to x and X is the
      % optimal solution.
      % We need to distinguish the following cases:
      % - linear constraints g(x) = ax - b. The derivative of such a
      %   constraint is g'(x) = a.
      % - second order constraints g(x[1],...,x[n]) = -x[1] + |(x[2],...,x[n])|
      %   the derivative of such a constraint is
      %     g'(x) = (-1, x[2]/|(x[2],...,x[n])|, ..., x[n]/|(x[2],...,x[n])|
      %   (here |.| denotes the Euclidean norm).
      % - bound constraints g(x) = -x for variables that are not explicitly
      %   contained in any second order cone constraint. The derivative for
      %   such a constraint is g'(x) = -1.
      % Note that it may happen that the derivative of a second order cone
      % constraint is not defined at the optimal solution X (this happens if
      % X=0). In this case we just skip the stationarity test.
      val = cplex.Model.obj;

      % Linear constraints.
      val = val - cplex.Model.A' *linpi;

      % Quadratic constraints.
      for q = 1:length(cplex.Model.qc)
         norm = 0;
         Q = cplex.Model.qc(q).Q;
         [m,n] = size(Q);
         for i = 1:m
            if Q(i,i) > 0
               norm = norm + x(i) * x(i);
            end
         end
         norm = sqrt(norm);
         if norm < tol
            fprintf('WARNING: Cannot test stationarity at non-differentiable point.\n');
            satisfied = 1;
            return;
         else
            for i = 1:m
               if Q(i,i) > 0
                  val(i) = val(i) + socppi(q) * x(i) / norm;
               elseif Q(i,i) < 0
                  val(i) = val(i) - socppi(q);
               end
            end
         end
      end
      % Handle variables that are not in quadratic constraints.
      for i = 1:length(constat)
         if constat(i) == -2
            val(i) = val(i) - dslack(i);
         end
      end

      % Test that all entries in val are 0.
      for i = 1:length(val)
         if abs(val(i)) > tol
            fprintf('Stationarity test failed at index %d: %f\n', i, val(i));
            return;
         end
      end

      satisfied = 1;
   end
end
