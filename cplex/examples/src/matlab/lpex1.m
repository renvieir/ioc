function  lpex1(varargin)
% Enter and optimize a linear programming problem
%
% Demonstrates different methods for creating a problem.
% The user has to choose the method via an argument:
%
%       lpex1('r')     generates the problem by adding rows
%       lpex1('c')     generates the problem by adding columns
%       lpex1('a')     generates the problem by adding arrays
%
% The LP problem solved in this example is
%   Maximize  x1 + 2 x2 + 3 x3
%   Subject to
%      - x1 +   x2 + x3 <= 20
%        x1 - 3 x2 + x3 <= 30
%   Bounds
%        0 <= x1 <= 40
%        0 <= x2
%        0 <= x3

% ---------------------------------------------------------------------------
% File: lpex1.m
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
   % Check the command line arguments
   if nargin ~= 1
      display ('Usage: lpex1(X)');
      display ('where X is one of the following options: ');
      display ('''r''          generate problem by row');
      display ('''c''          generate problem by column');
      display ('''a''          generate problem by arrays');
      display ('Exiting...')
      return;
   else
      m = varargin{1};
   end
   
   % Initialize the CPLEX object
   cplex = Cplex('lpex1');
   
   % Now populate the problem with the data
   if m == 'c',
      % Use addCols to populate the model
      cplex.Model.sense = 'maximize';
      cplex.addRows([-inf; -inf], [], [20; 30]);
      cplex.addCols(1, [-1;  1], 0, 40);
      cplex.addCols(2, [ 1; -3], 0, inf);
      cplex.addCols(3, [ 1;  1], 0, inf);
      
   elseif m == 'r',
      % Use addRows to populate the model
      cplex.Model.sense = 'maximize';
      cplex.addCols([1; 2; 3], [], [0; 0; 0], [40; inf; inf]);
      cplex.addRows(-inf, [-1  1 1], 20);
      cplex.addRows(-inf, [ 1 -3 1], 30);
   else
      % Use arrays to populate the model
      cplex.Model.sense = 'maximize';
      cplex.Model.obj   = [1; 2; 3];
      cplex.Model.lb    = [0; 0; 0];
      cplex.Model.ub    = [40; inf; inf];
      cplex.Model.A     = [-1 1 1; 1 -3 1];
      cplex.Model.lhs   = [-inf; -inf];
      cplex.Model.rhs   = [20; 30];
   end
   
   % Optimize the problem
   cplex.solve();
   
   % Write the solution
   fprintf ('\nSolution status = %s\n',cplex.Solution.statusstring);
   fprintf ('Solution value = %f\n',cplex.Solution.objval);
   disp ('Values = ');
   disp (cplex.Solution.x');
catch m
   disp(m.message);
end
end
