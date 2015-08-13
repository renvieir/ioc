function miqpex1()
% Enter and optimize a mixed integer quadratic programming problem
%
% This function fills in the data structures for the mixed integer
% quadratic program:
%
%       Maximize
%        obj: x1 + 2 x2 + 3 x3 + x4
%             - 0.5 ( 33x1*x1 + 22*x2*x2 + 11*x3*x3
%                    -  12*x1*x2 - 23*x2*x3 )
%       Subject To
%        c1: - x1 + x2 + x3 + 10x4  <= 20
%        c2: x1 - 3 x2 + x3         <= 30
%        c3:       x2       - 3.5x4  = 0
%       Bounds
%        0 <= x1 <= 40
%        0 <= x2
%        0 <= x3
%        0 <= x4 <= 3
%       Integers
%         x4
%       End

% ---------------------------------------------------------------------------
% File: miqpex1.m
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
   % Initialize the CPLEX object
   cplex = Cplex('miqpex1');
   
   % Fill in the data for the problem using populatebyrow
   populatebyrow (cplex);
   
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
   disp(m.message);
end
   function populatebyrow(cplex)
      cplex.Model.sense = 'maximize';
      cplex.addCols([1 2 3 1]', [], [0 0 0 0]', [40 inf  inf 3]', 'CCCI');
      cplex.addRows(-inf, [-1 1 1 10], 20);
      cplex.addRows(-inf, [1 -3 1  0], 30);
      cplex.addRows(0, [0 1 0 -3.5], 0);
      cplex.Model.Q = [-33     6     0  0;
         6   -22  11.5  0;
         0  11.5   -11  0;
         0     0     0  0];
   end
end
