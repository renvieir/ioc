function cplexmilpex
% Use the function cplexmilp to solve a mixed-integer linear programming problem
%
% The MILP problem solved in this example is
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

% ---------------------------------------------------------------------------
% File: cplexmilpex.m
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
   % Since cplexmilp solves minimization problems and the problem
   % is a maximization problem, negate the objective
   f     = [-1 -2 -3 -1]';
   Aineq = [-1  1  1 10;
      1 -3  1  0];
   bineq = [20 30]';
   
   Aeq   = [0  1  0  -3.5];
   beq   =  0;
   
   lb    = [0;    0;   0; 2];
   ub    = [40; inf; inf; 3];
   ctype = 'CCCI';
   
   options = cplexoptimset;
   options.Display = 'on';
   
   [x, fval, exitflag, output] = cplexmilp (f, Aineq, bineq, Aeq, beq,...
      [ ], [ ], [ ], lb, ub, ctype, [ ], options);
   
   fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
   fprintf ('Solution value = %f \n', fval);
   disp ('Values =');
   disp (x');
catch m
   throw (m);
end
end
