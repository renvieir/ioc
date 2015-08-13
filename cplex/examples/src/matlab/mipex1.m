function mipex1(varargin)
% Enter and optimize a mixed integer programming problem
%
% Demonstrates different methods for creating a problem.
% The user has to choose the method via a function argument:
%
%       mipex1('r')     generates the problem by adding rows
%       mipex1('c')     generates the problem by adding columns
%       mipex1('a')     generates the problem by adding arrays
%
% The MIP problem solved in this example is
%   Maximize  x1 + 2 x2 + 3 x3 + x4
%   Subject to
%      - x1 +   x2 + x3 + 10 x4 <= 20
%        x1 - 3 x2 + x3         <= 30
%               x2      - 3.5x4  = 0
%   Bounds
%        0 <= x1 <= 40
%        0 <= x2
%        0 <= x3
%        2 <= x4 <= 3
%   Integers
%       x4

%

% ---------------------------------------------------------------------------
% File: mipex1.m
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
      display ('Usage: mipex1(X)');
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
   cplex = Cplex('mipex1');
   cplex.Model.sense = 'maximize';
   
   % Now populate the problem with the data
   if m == 'c',
      % Use addCols to populate model
      cplex.addRows([-inf; -inf; 0], [], [  20;   30; 0]);
      cplex.addCols(1, [-1;  1;  0],   0,  40, 'C');
      cplex.addCols(2, [ 1; -3;  1],   0, inf, 'C');
      cplex.addCols(3, [ 1;  1;  0],   0, inf, 'C');
      cplex.addCols(1, [ 10; 0; -3.5], 2,   3, 'I');
      
   elseif m == 'r',
      % Use addRows to populate model
      cplex.addCols([1; 2; 3; 1], [], [0; 0; 0; 2], [40; inf; inf; 3],'CCCI');
      cplex.addRows(-inf, [-1  1 1   10], 20);
      cplex.addRows(-inf, [ 1 -3 1    0], 30);
      cplex.addRows(   0, [ 0  1 0 -3.5],  0);
      
   else
      % Use arrays to populate the model
      cplex.Model.obj   = [ 1;   2;   3; 1];
      cplex.Model.lb    = [ 0;   0;   0; 2];
      cplex.Model.ub    = [40; inf; inf; 3];
      cplex.Model.ctype = 'CCCI';
      cplex.Model.A =  [-1  1  1  10;
         1 -3  1   0;
         0  1  0  -3.5];
      cplex.Model.lhs = [-inf; -inf; 0];
      cplex.Model.rhs = [  20;   30; 0];
   end
   
   % Optimize the problem
   cplex.solve();
   
   % Write the solution
   fprintf ('\nSolution status = %s \n', cplex.Solution.statusstring);
   fprintf ('Solution value = %f \n', cplex.Solution.objval);
   disp ('Values =');
   disp (cplex.Solution.x);
   disp ('Slacks =');
   disp (cplex.Model.rhs - cplex.Solution.ax);
catch m
   throw (m);
end
end
