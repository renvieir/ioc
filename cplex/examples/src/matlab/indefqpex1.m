function indefqpex1()
% Enter and optimize an indefinite quadratic programming problem
%
% This function fills in the data structures for the quadratic program:
%
%       Minimize
%        obj: - 0.5 (-3 * xˆ2 - 3 * yˆ2 - 1 * x * y)
%
%       Subject To
%        c1: -x + y >= 0
%        c2:  x + y >= 0
%       Bounds
%        -1 <= x <= 1
%         0 <= y <= 1
%       End

% ---------------------------------------------------------------------------
% File: indefqpex1.m
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
   cplex = Cplex('indefqpex1');
   cplex.Model.sense = 'minimize';
   
   % Fill in the data for the problem by populatebyrow
   populatebyrow();

   % When a non-convex objective function is present, CPLEX will raise
   % an exception unless the parameter solutiontarget is set
   % to accept first-order optimal solutions
   cplex.Param.solutiontarget.Cur = 2;
   
   % CPLEX may converge to either local optimum 
   solveanddisplay();

   % Add a constraint that cuts off the solution at (-1, 1)
   cplex.addRows(0, [1 0], inf);
   solveanddisplay();


   % Remove the newly added constraint and add a new constraint
   % with the opposite sense to cut off the solution at (1, 1)
   cplex.delRows([3]);
   cplex.addRows(-inf, [1 0], 0);
   solveanddisplay();


   % Finally, write a copy of the problem to a file
   cplex.writeModel('indefqpex1.lp');
catch m
   throw (m);
end

   function populatebyrow()
      cplex.addCols([0 0]', [], [-1 0]', [1 1]');
      cplex.Model.Q   = [-3.0   -0.5; ...
                         -0.5   -3.0];
      cplex.addRows(0, [-1  1], inf);
      cplex.addRows(0, [ 1  1], inf);
   end

   function solveanddisplay()
      % Optimize the problem
      cplex.solve();
      
      % Write the solution
      fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
      fprintf ('Solution value = %f\n', cplex.Solution.objval);
      disp ('Values = ');
      disp (cplex.Solution.x);
      disp ('Slacks =');
      disp (cplex.Model.rhs - cplex.Solution.ax);
      disp ('Duals =');
      disp (cplex.Solution.dual);
      disp ('Reduced Costs =');
      disp (cplex.Solution.reducedcost);
      
   end
end
