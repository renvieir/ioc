function populate(varargin)
% Read and generate multiple solutions to a MIP problem
%
% To run this example, a function argument is required.
%    i.e.,   populate(filename)
%    where
%      filename is the name of the file, with .mps, .lp, or .sav extension

% ---------------------------------------------------------------------------
% File: populate.m
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
      display ('Usage: populate(filename)');
      display ('where filename is a file with extension: ');
      display ('MPS, SAV, or LP (lower case is allowed)');
      display ('Exiting...')
      return;
   else
      m = varargin{1};
   end
   eps = 1.0E-10;
   
   % Initialize the CPLEX object
   cplex = Cplex();
   
   % Now read the file and copy the data into cplex.Model
   cplex.readModel(m);
   
   % Set the solution pool relative gap parameter to obtain solutions
   % of objective value within 10% of the optimal
   cplex.Param.mip.pool.relgap.Cur = 0.1;
   
   % Optimize the problem and obtain multiple solutions
   cplex.populate();
   
   % Display the solution
   fprintf ('Solution status: %d\n', cplex.Solution.status);
   fprintf ('Incumbent objective value: %g\n', cplex.Solution.objval);
   disp ('Incumbent values');
   disp (cplex.Solution.x);
   
   numsoln = size (cplex.Solution.pool.solution, 1);
   fprintf ('The solution pool contains %d solutions\n', numsoln);
   
   % Some solutions are deleted from the pool because of the solution
   % pool relative gap parameter
   fprintf(...
'%d solutions were removed due to the solution pool relative gap parameter.\n',...
      cplex.Solution.pool.numreplaced);
   
   fprintf ('In total, %d solutions were generated.\n',...
      numsoln + cplex.Solution.pool.numreplaced);
   
   % Get the average objective value of solutions in the solution
   % pool
   fprintf ('The average objective value of the solutions is %g\n', ...
      cplex.Solution.pool.meanobj);
   
   fprintf ('Solution       Solution  Values Differing\n');
   fprintf ('  Number      Objective    from Incumbent\n');
   fprintf ('--------  -------------  ----------------\n');
   for i = 1:numsoln
      fprintf ('%8d  %13g', i, cplex.Solution.pool.solution(i).objval);
      numdiff = 0;
      for j = 1:length(cplex.Solution.x)
         if abs (cplex.Solution.pool.solution(i).x(j) - cplex.Solution.x(j)) ...
             > eps
            numdiff = numdiff + 1;
         end
      end
      fprintf ('  %16d\n', numdiff);
   end
catch m
   throw (m);
end
end
