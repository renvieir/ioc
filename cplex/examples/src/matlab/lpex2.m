function lpex2(varargin)
% Read and optimize a linear programming problem
%
% To run this example, a function argument is required.
%    i.e.,   lpex2(filename)
%    where
%      filename is the name of the problem file, with .mps, .lp, or .sav
%      extension

% ---------------------------------------------------------------------------
% File: lpex2.m
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
      display ('Usage: lpex2(filename)');
      display ('where filename is a problem file with extension: ');
      display ('MPS, SAV, or LP (lower case is allowed)');
      display ('Exiting...')
      return;
   else
      m = varargin{1};
   end
   
   % Initialize the CPLEX object
   cplex = Cplex('lpex2');
   
   % Now read the file and copy the data into the created cplex.Model
   cplex.readModel(m);
   
   % Optimize the problem
   cplex.solve();
   
   % Write the solution
   fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
   fprintf ('Solution value = %f\n', cplex.Solution.objval);
   fprintf ('Values = \n');
   disp (cplex.Solution.x);
catch m
   throw (m);
end
end
