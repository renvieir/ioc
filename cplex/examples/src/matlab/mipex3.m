function mipex3(varargin)
% Enter and optimize a MIP problem with SOS sets
%
% Demonstrates different methods for creating a problem.
% The user has to choose the method via a function argument:
%
%       mipex3('r')     generates the problem by adding rows
%       mipex3('c')     generates the problem by adding columns
%       mipex3('n')     generates the problem by adding
%                     a list of coefficients

% ---------------------------------------------------------------------------
% File: mipex3.m
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
      display ('Usage: mipex3(X)');
      display ('where X is one of the following options: ');
      display ('''r''          generate problem by row');
      display ('''c''          generate problem by column');
      display ('''n''          generate problem by nonzero');
      display ('Exiting...')
      return;
   else
      m = varargin{1};
   end
   
   % Initialize the CPLEX object
   cplex = Cplex('mipex3');
   cplex.Model.sense = 'maximize';
   
   % Now populate the problem with the data
   if m == 'c',
      % Use addCols to populate the model
      cplex.addRows([-inf; -inf; 0], [], [20; 30; 0]);
      cplex.addCols(1, [-1;  1 ; 0], 0, 40);
      cplex.addCols(2, [ 1; -3 ; 1], 0, inf, 'I');
      cplex.addCols(3, [ 1;  1 ; 0], 0, inf, 'I');
      cplex.addCols(1, [ 10; 0 ; -3.5], 2, 3, 'I');
   elseif m == 'r',
      % Use addRows to populate the model
      cplex.addCols([1; 2; 3; 1], [], [0; 0; 0; 2], [40; inf; inf; 3],'CIII');
      cplex.addRows(-inf, [-1  1 1  10], 20);
      cplex.addRows(-inf, [ 1 -3 1   0], 30);
      cplex.addRows(   0, [ 0  1 0 -3.5],  0);
   else
      % Use arrays to populate the model
      cplex.Model.obj   = [ 1;   2;   3; 1];
      cplex.Model.lb    = [ 0;   0;   0; 2];
      cplex.Model.ub    = [40; inf; inf; 3];
      cplex.Model.ctype = 'CIII';
      cplex.Model.A     =  [-1  1 1   10;
         1 -3 1    0;
         0  1 0 -3.5];
      cplex.Model.lhs   = [-inf; -inf; 0];
      cplex.Model.rhs   = [  20;   30; 0];
   end
   
   % Use addSOSs to add SOS
   cplex.addSOSs('1', [2 3]', [25 18]', {'sos1(1)'});
   
   % Add Order to model
   cplex.Order.ind = [1;3];
   cplex.Order.pri = [8;7];
   cplex.Order.dir = [1;-1];
   cplex.writeOrder('mipex3.ord');
   
   % Optimize the problem and obtain solution.
   cplex.solve();
   
   % Write the output to the screen.
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
