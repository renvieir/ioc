function qcpex1()
% Enter and optimize a quadratically constrained programming problem
%
% This function fills in the data structures for the quadratic constraint
% program:
%
%       Maximize
%        obj: x1 + 2 x2 + 3 x3
%               - 0.5 ( 33x1*x1 + 22*x2*x2 + 11*x3*x3
%                    -  12*x1*x2 - 23*x2*x3 )
%       Subject To
%        c1: - x1 + x2 + x3 <= 20
%        c2: x1 - 3 x2 + x3 <= 30
%        q1: [ x1^2 + x2^2 + x3^3 ] <= 1.0
%       Bounds
%        0 <= x1 <= 40
%        0 <= x2
%        0 <= x3
%       End

% ---------------------------------------------------------------------------
% File: qcpex1.m
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
   cplex = Cplex('qcpex1');
   cplex.Model.sense = 'maximize';
   
   % Fill in the data for the problem with populatebyrow
   populatebyrow();
   
   % Optimize the problem
   cplex.solve();
   
   % Write the solution
   fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
   fprintf ('Solution value = %f\n', cplex.Solution.objval);
   disp ('Values = ');
   disp (cplex.Solution.x');
   disp ('Slacks = ');
   disp (cplex.Model.rhs - cplex.Solution.ax);
   
   % Finally, write a copy of the problem to a file
   cplex.writeModel('qc.lp');
catch m
   throw (m);
end

   function populatebyrow()
      cplex.addCols([1 2 3]', [], [0; 0; 0], [40; inf; inf]);
      cplex.Model.Q   = [-33     6   0; ...
         6   -22  11.5; ...
         0  11.5 -11];
      cplex.addRows(-inf, [-1  1 1], 20);
      cplex.addRows(-inf, [ 1 -3 1], 30);
      cplex.addQCs([0 0 0]', [1 0 0;0 1 0;0 0 1], 'L', 1.0);
   end
end
