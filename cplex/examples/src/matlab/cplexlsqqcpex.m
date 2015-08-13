function  cplexlsqqcpex
% Use the function cplexlsqqcp to solve a constrained least squares problem
% Some variables are binary and one constraint is quadratic.

% ---------------------------------------------------------------------------
% File: cplexlsqqcpex.m
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
   C = [0.9501    0.7620    0.6153    0.4057
      0.2311    0.4564    0.7919    0.9354
      0.6068    0.0185    0.9218    0.9169
      0.4859    0.8214    0.7382    0.4102
      0.8912    0.4447    0.1762    0.8936];
   d = [0.0578
      0.3528
      0.8131
      0.0098
      0.1388];
   Aineq = [0.2027    0.2721    0.7467    0.4659
      0.1987    0.1988    0.4450    0.4186
      0.6037    0.0152    0.9318    0.8462];
   bineq = [0.5251
      0.2026
      0.6721];
   lb = -0.1 * ones(4, 1);
   ub =  2.0 * ones(4, 1);
   
   l = [0 0 0 0]';
   r = 1;
   Q = [1  0  0 0
      0  1  0 0
      0  0  1 0
      0  0  0 1];
   
   options = cplexoptimset;
   options.Display = 'on';
   
   [x, resnorm, residual, exitflag, output] = ...
      cplexlsqqcp (C, d, Aineq, bineq,[ ], [ ], l, Q, r, ...
  lb, ub, [ ], options);
   
   fprintf ('\nSolution status = %s\n', output.cplexstatusstring);
   disp ('Values =');
   disp (x');
   disp ('resnorm  =');
   disp (resnorm);
   disp ('residual =');
   disp (residual');
catch m
   throw (m);
end
end
