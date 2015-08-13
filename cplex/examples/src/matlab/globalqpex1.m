function globalqpex1(varargin)
% Read and optimize a convex nonconvex (mixed integer) quadratic programming problem
% with convex, first order or global optimizer.
%
% To run this example, a function argument is required.
%    That is:   globalqpex1(filename,solutiontarget)
%    where:
%       filename        is the name of the problem file,
%                       with .mps, .lp, or .sav extension
%       solutiontarget  is one of the characters
%                       c       for convex QP
%                       f       for first order solution (only continuous)
%                       g       for global optimum

% ---------------------------------------------------------------------------
% File: globalqpex1.m
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
   if nargin ~= 2
      display ('Usage: globalqpex1(filename,solutiontarget)');
      display ('where');
      display ('filename        is a problem file with extension ');
      display ('                MPS, SAV, or LP (lower case is allowed)');
      display ('solutiontarget  is one of the characters');
      display ('                c       for convex QP           ');
      display ('                f       for first order solution');
      display ('                g       for global optimum      ');
      display ('Exiting...')
      return;
   else
      m = varargin{1};
      t = varargin{2};
   end
   
   % Initialize the CPLEX object
   cplex = Cplex('globalqpex1');
   
   % Now read the file and copy the data into cplex.Model
   cplex.readModel(m);

   % set the solution target parameter
   if     ( t == 'c' )
      cplex.Param.solutiontarget.Cur = 1;
   elseif ( t == 'f' )
      cplex.Param.solutiontarget.Cur = 2;
   elseif ( t == 'g' )
      cplex.Param.solutiontarget.Cur = 3;
   end
   
   % Try to optimize the problem
   try
      cplex.solve();
   catch m
      if findstr ('5002', m.message) > 0
         disp ('Problem is indefinite.');
         return;
      else
         throw m;
      end
   end
   
   % Write the solution
   fprintf ('\nSolution status = %s \n',cplex.Solution.statusstring);
   fprintf ('Solution value = %f \n',cplex.Solution.objval);
   disp ('Values = ');
   disp (cplex.Solution.x');
catch m
   if ( t ~= 'c' )
      disp (m.message);
end
end
