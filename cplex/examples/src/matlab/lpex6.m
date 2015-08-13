function  lpex6()
% Use an advanced basis to start an LP optimization

% ---------------------------------------------------------------------------
% File: lpex6.m
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
   cplex = Cplex('lpex6');
   
   % The row statuses are all nonbasic for this problem
   cplex.Start.basis.rowstat = [0 0];
   
   % We assume we know the optimal basis.  Variables 2 and 3 are basic,
   % while variable 1 is at its upper bound
   cplex.Start.basis.colstat = [2 1 1]';
   
   % This function builds by column the linear program:
   %
   % Maximize
   %    obj: x1 + 2 x2 + 3 x3
   %    Subject To
   %    c1: - x1 + x2 + x3 <= 20
   %    c2: x1 - 3 x2 + x3 <= 30
   %    Bounds
   %    0 <= x1 <= 40
   %    0 <= x2
   %    0 <= x3
   %    End
   
   cplex.Model.sense = 'maximize';
   cplex.addRows([-Inf -Inf]', [], [20 30]');
   cplex.addCols([1 2 3]', [-1 1 1; 1 -3 1], [0 0 0]', [40 inf inf]');
   
   % Optimize the problem
   cplex.solve();
   
   % Write the solution
   fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
   fprintf ('Solution value = %f\n', cplex.Solution.objval);
   
   if ~isfield (cplex.Solution, 'itcnt')
      disp ('Iteration count = 0');
   else
      fprintf ('Iteration count = %d\n', cplex.Solution.itcnt);
   end
   
   disp ('Values:');
   disp (cplex.Solution.x');
   disp ('Duals:')
   disp (cplex.Solution.dual');
   disp ('Axs:')
   disp (cplex.Solution.ax');
   disp ('Reducedcosts:')
   disp (cplex.Solution.reducedcost');
catch m
   throw (m);
end
end
