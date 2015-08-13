function diet(m)
% Minimize the cost of a diet meeting nutritional constraints
%
% An argument is required:
%   'r' to build the problem by rows
%   'c' to build the problem by columns
%   'i' to build the problem by arrays and to require integral food values
%
% Data is read from the file '../../data/diet.dat'

% ---------------------------------------------------------------------------
% File: diet.m
% Version 12.6.1
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
% Copyright IBM Corporation 2008, 2014. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
% ---------------------------------------------------------------------------

% Input data:
% foodMin[j]      minimum amount of food j to use
% foodMax[j]      maximum amount of food j to use
% foodCost[j]     cost for one unit of food j
% nutrMin[i]      minimum amount of nutrient i
% nutrMax[i]      maximum amount of nutrient i
% nutrPer[i][j]   nutrition amount of nutrient i in food j
%
% Modeling variables:
% Buy[j]          amount of food j to purchase
% Objective:
% minimize sum(j) Buy[j] * foodCost[j]
%
% Constraints:
% forall foods i: nutrMin[i] <= sum(j) Buy[j] * nutrPer[i][j] <= nutrMax[j]

try
   % Read the data
   [foodCost, foodMin, foodMax, nutrMin, nutrMax, nutrPer] ...
      = inputdata ('../../data/diet.dat');
   
   nFood = length (foodCost);
   
   % Build model
   cplex = Cplex('diet');
   cplex.Model.sense = 'minimize';
   
   % The variables are the amount of each food to purchase
   
   % Populate by column
   if m == 'c'
      cplex.addRows(nutrMin, [], nutrMax);
      for k = 1:nFood;
         cplex.addCols(foodCost(k), nutrPer(:,k), foodMin(k), foodMax(k));
      end
      
      % Populate by row
   elseif m == 'r'
      cplex.addCols(foodCost, [], foodMin, foodMax);
      
      % We can add rows as a set
      cplex.addRows(nutrMin, nutrPer, nutrMax);
      
      % Or, we can add them one by one
      %     for i = 1:nNutr
      %         cplex.addRows(nutrMin(i), nutrPer(i,:), nutrMax(i));
      %     end
      
      % Populate using arrays
      % Also require foods to be integral quantities
   elseif m == 'i'
      cplex.Model.obj   = foodCost;
      cplex.Model.lb    = foodMin;
      cplex.Model.ub    = foodMax;
      cplex.Model.lhs   = nutrMin;
      cplex.Model.rhs   = nutrMax;
      cplex.Model.A     = nutrPer;
      cplex.Model.ctype = char (ones ([1, nFood]) * ('I'));
   end
   
   % Solve model
   cplex.solve();
   
   % Display solution
   fprintf ('\nSolution status = %s\n',cplex.Solution.statusstring);
   fprintf ('Diet Cost: %f\n', cplex.Solution.objval);
   disp ('Buy Food:');
   disp (cplex.Solution.x(1:nFood)');
catch m
   throw (m);
end
