function  inout1()
% Solve a production planning problem minimizing cost
%
% A company has to produce 3 products, using 2 resources.
%
% Each resource has a limited capacity.
% Each product consumes a given number of machines.
% Each product has a production cost (the inside cost).
% Both products can also be bought outside the company at a given
% cost (the outside cost).
%
% Minimize the total cost so that the company exactly meets the demand.

% ---------------------------------------------------------------------------
% File: inout1.m
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
   nbProds      = 3;
   nbResources  = 2;
   consumption  = [0.5 0.4 0.3
      0.2 0.4 0.6];
   demand       = [100 200 300];
   insideCost   = [0.6 0.8 0.3];
   outsideCost  = [0.8 0.9 0.4];
   capacity     = [20 40];
   
   % Build model
   cplex = Cplex('inout1');
   
   % The first 3 variables are the amount of each product to produce
   % and the next 3 variables are the amount of each product to buy.
   
   % The objective is to minimize cost of producing and buying
   % Variables are nonnegative
   cplex.Model.sense = 'minimize';
   cplex.addCols([insideCost outsideCost]', ...
                 [], zeros (nbResources*nbProds, 1), ...
                 ones (nbResources*nbProds, 1) * inf);
   
   % Meet the demand for each product
   for i = 1:nbProds
      
      % Start with an all-zero constraint
      constr = zeros (1, length (cplex.Model.obj)) ;
      
      % Manufacture product inside
      constr(i) = 1;
      
      % Buy product outside
      constr(i+nbProds) = 1;
      
      % Sum of inside and outside product must equal demand
      cplex.addRows(demand(i), constr, demand(i));
   end
   
   % Respect the capacity constraint for each resource
   cplex.addRows(-inf, [consumption(1,:) zeros(1, nbProds)], capacity(1));
   cplex.addRows(-inf, [consumption(2,:) zeros(1, nbProds)], capacity(2));
   
   % Solve model
   cplex.solve();
   
   % Display solution
   fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
   fprintf ('\nCost: %f\n\n', cplex.Solution.objval);
   
   disp ('Produce inside:');
   disp (cplex.Solution.x(1:nbProds)');
   
   disp ('Buy outside:');
   disp (cplex.Solution.x((nbProds+1):nbResources*nbProds)');
catch m
   throw (m);
end
end
%
% Cost:
%    372
%
% Produce inside:
%     40
%      0
%      0
%
% Buy outside:
%     60
%    200
%    300

