function lpex3()
% Solve a problem in two stages
%
% Modified example from Chvatal, "Linear Programming", Chapter 26.
%   minimize  c*x
%   subject to  Hx = d
%               Ax = b
%               l <= x <= u
%   where
%
%   H = ( -1  0  1  0  1  0  0  0 )  d = ( -3 )
%       (  1 -1  0  1  0  0  0  0 )      (  1 )
%       (  0  1 -1  0  0  1 -1  0 )      (  4 )
%       (  0  0  0 -1  0 -1  0  1 )      (  3 )
%       (  0  0  0  0 -1  0  1 -1 )      ( -5 )
%
%   A = (  2  1 -2 -1  2 -1 -2 -3 )  b = (  4 )
%       (  1 -3  2  3 -1  2  1  1 )      ( -2 )
%
%   c = ( -9  1  4  2 -8  2  8 12 )
%   l = (  0  0  0  0  0  0  0  0 )
%   u = ( 50 50 50 50 50 50 50 50 )
%
%  Treat the constraints with A as the complicating constraints, and
%  the constraints with H as the "simple" problem.
%
%  The idea is to solve the simple problem first, and then add the
%  constraints for the complicating constraints and solve with dual.

% ---------------------------------------------------------------------------
% File: lpex3.m
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
   c = [-9  1  4  2 -8  2  8 12]';
   H = [-1  0  1  0  1  0  0  0;
      1 -1  0  1  0  0  0  0;
      0  1 -1  0  0  1 -1  0;
      0  0  0 -1  0 -1  0  1;
      0  0  0  0 -1  0  1 -1];
   A = [ 2  1 -2 -1  2 -1 -2 -3;
      1 -3  2  3 -1  2  1  1];
   
   d = [-3  1  4  3 -5]';
   b = [ 4 -2]';
   l = zeros(8, 1);
   u = ones(8, 1) * 50;
   
   % Initialize the CPLEX object
   cplex = Cplex('lpex3');
   
   cplex.Model.sense = 'minimize';
   cplex.addCols(c, [], l, u);
   
   % Now add the rows to the problem
   cplex.addRows(d, H, d);
   
   % Solve model with the network algorithm
   cplex.Param.lpmethod.Cur = 3;
   cplex.solve();
   
   % Display solution
   fprintf ('\nAfter network optimization, objective is : %f\n',...
            cplex.Solution.objval);
   
   % Add the complicating constraints
   cplex.addRows(b, A, b);
   
   % Re-solve the model with the dual simplex method
   cplex.Param.lpmethod.Cur = 2;
   cplex.solve();
   
   % Display solution
   fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
   fprintf ('\nAfter dual optimization, objective is: %f\n',...
            cplex.Solution.objval);
   disp ('Solution values are ');
   disp (cplex.Solution.x');
catch m
   throw (m);
end
end

