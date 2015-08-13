function  inout3()
% Solve a two-step production planning problem minimizing cost
%
% A company has to produce 3 products, using 2 resources.
%
% Each resource has a limited capacity.
% Each product consumes a given number of machines.
% Each product has a production cost (the inside cost).
% Both products can also be brought outside the company at a given cost
% (the outside cost)
%
% Minimize external production given product demand, a cost constraint,
% and minimum internal production constraints.

% ---------------------------------------------------------------------------
% File: inout3.m
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
   nbProds = 3;
   consumption0 = [0.5 0.4 0.3];
   consumption1 = [0.2 0.4 0.6];
   capacity = [20 40];
   demand = [100 200 300];
   insidecost = [0.6 0.8 0.3];
   outsidecost = [0.8 0.9 0.4];
   
   % Build model
   cplex=Cplex('inout3');
   cplex.Model.sense = 'minimize';
   obj = [insidecost outsidecost]';
   lb = [10 10 10 -inf -inf -inf]';
   ub = (ones (1, 2*nbProds) * inf)';
   
   cplex.addCols(obj, [], lb, ub);
   
   % Must respect capacity constraint for each resource
   cplex.addRows(-inf, [consumption0 zeros(1,length(consumption1))], ...
                 capacity(1));
   cplex.addRows(-inf, [consumption1 zeros(1,length(consumption1))], ...
                 capacity(2));
   
   % Must meet demand for each product
   for i = 1:nbProds
      constr = zeros (1, length (cplex.Model.obj)) ;
      constr(i) = 1;
      constr(i + nbProds) = 1;
      cplex.addRows(demand(i), constr, demand(i));
   end
   
   % Solve model
   cplex.solve();
   fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
   cost=cplex.Solution.objval;
   
   % New objective: minimize outside production
   cplex.Model.obj = [zeros(1,nbProds) ones(1,nbProds)]';
   
   % New constraint: cost must be no more than 10% over minimum
   cplex.addRows(0, [insidecost outsidecost], 1.1*cost);
   cplex.writeModel('inout3.lp');
   cplex.solve();
   
   % Display solution
   fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
   disp ('insidecost:');
   disp (cplex.Solution.x(1:nbProds));
   disp ('outsidecost:');
   disp (cplex.Solution.x((nbProds+1):nbProds+nbProds));
catch m
   throw (m);
end
end
% insidecost:
%    10.0000
%    10.0000
%    36.6667
%
% outsidecost:
%    90.0000
%   190.0000
%   263.3333
