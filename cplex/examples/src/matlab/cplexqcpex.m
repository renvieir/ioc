function cplexqcpex()
% Use the function cplexqcp to solve a quadratically constrained programming problem
%
% The QCP problem solved in this example is
%   Maximize  x1 + 2 x2 + 3 x3
%             - 0.5 ( 33x1*x1 + 22*x2*x2 + 11*x3*x3 - 12*x1*x2 - 23*x2*x3 )
%   Subject To
%      - x1 +   x2 + x3 <= 20
%        x1 - 3 x2 + x3 <= 30
%     x1*x1 + x2*x2 + x3*x3 <= 1.0
%   Bounds
%        0 <= x1 <= 40
%        0 <= x2
%        0 <= x3

% ---------------------------------------------------------------------------
% File: cplexqcpex.m
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
   % Since cplexqcp solves minimization problems and the problem
   % is a maximization problem, negate the objective
   H  = [33   6    0;   ...
      6  22   11.5; ...
      0  11.5 11];
   f  = [-1 -2 -3]';
   
   Aineq  = [-1  1  1;
      1 -3  1];
   bineq  = [20 30]';
   
   l = [0 0 0]';
   r = 1;
   Q = [1  0  0
      0  1  0
      0  0  1];
   
   lb = [ 0   0   0]';
   ub = [40 inf inf]';
   
   options = cplexoptimset;
   options.Display = 'on';
   
   [x, fval, exitflag, output] = ...
      cplexqcp (H, f, Aineq, bineq, [ ], [ ], l, Q, r, lb, ub, [ ], options);
   
   fprintf ('\nSolution status = %s \n', output.cplexstatusstring);
   fprintf ('Solution value = %f \n', fval);
   disp ('Values =');
   disp (x');
   
catch m
   throw (m);
end
end
