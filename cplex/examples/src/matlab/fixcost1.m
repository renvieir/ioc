function fixcost1()
% Minimize the cost of producing a product
%
% A company must produce a product on a set of machines.
% Each machine has limited capacity.
% Producing a product on a machine has both a fixed cost and a cost per
% unit of production.
%
% Minimize the sum of fixed and variable costs so that the company exactly
% meets demand.

% ---------------------------------------------------------------------------
% File: fixcost1.m
% Version 12.6.1
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
% Copyright IBM Corporation 2008, 2014. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
% ---------------------------------------------------------------------------

try
   % Define data
   nbMachines = 6;
   cost = [15.0 20.0 45.0 64.0 12.0 56.0]';
   capacity = [100.0, 20.0, 405.0, 264.0, 12.0, 256.0]';
   fixedCost = [1900.0, 820.0, 805.0, 464.0, 3912.00, 556.0]';
   demand = 22.0;
   xlb = zeros (nbMachines, 1);
   xub = ones (nbMachines, 1) * Inf;
   fusedlb = zeros (nbMachines, 1);
   fusedub = ones (nbMachines, 1);
   
   % Biuld model
   cplex=Cplex('fixcost1');
   cplex.Model.sense = 'minimize';
   
   % Objective: minimize the sum of fixed and variable costs
   obj = [cost;fixedCost];
   lb = [xlb;fusedlb];
   ub = [xub;fusedub];
   ctype = 'CCCCCCIIIIII';
   cplex.addCols(obj, [], lb, ub, ctype);
   
   % Constraint: meet demand
   cplex.addRows(demand, [ones(1, nbMachines) zeros(1, nbMachines)], demand);
   
   % Constraint: respect capacity constraint on machine 'i'
   cplex.addRows(ones (nbMachines, 1) * -inf,...
      [eye(nbMachines), zeros(nbMachines, nbMachines)], capacity);
   
   % Constraint: only produce product on machine 'i' if it is 'used'
   % (to capture fixed cost of using machine 'i')
   cplex.addRows(ones (nbMachines, 1) * -inf,...
      [eye(nbMachines), -diag(capacity)],...
       zeros(nbMachines, 1));
   
   % Solve model
   cplex.solve();
   
   % Display solution
   fprintf ('\nSolution status = %s\n',cplex.Solution.statusstring);
   fprintf ('Obj: %f\n', cplex.Solution.objval);
   x = cplex.Solution.x(1:nbMachines);
   fused = cplex.Solution.x(nbMachines+1:nbMachines+nbMachines);
   for i = 1:nbMachines
      if fused(i) > 1e-6
         fprintf ('E%d is used for %d\n', i, x(i));
      end
   end
catch m
   disp(m.message);
end

end
