function  cplexlsqnonnegmilpex
% Use the function cplexlsqnonnegmilp to solve a nonnegative least squares
% problem. The variables are binary.

% ---------------------------------------------------------------------------
% File: cplexlsqnonnegmilpex.m
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
    C = [0.0372    0.2869
        0.6861    0.7071
        0.6233    0.6245
        0.6344    0.6170];
    d = [0.8587
        0.1781
        0.0747
        0.8405];
    Aineq = [0.2027    0.2721
        0.1987    0.1988
        0.6037    0.0152];
    bineq = [0.5251
        0.2026
        0.6721];
    ctype = 'BB';
    
    options = cplexoptimset;
    options.Display = 'on';
    
    [x, resnorm, residual, exitflag, output] = ...
        cplexlsqnonnegmilp (C, d, Aineq, bineq, [], [], [], [], [], ctype, ...
        [], options);
    
    fprintf ('\nSolution status = %s\n', output.cplexstatusstring);
    disp ('Values =');
    disp (x');
    disp ('resnorm =');
    disp (resnorm);
    disp ('residual =');
    disp (residual');
catch m
    throw (m);
end
end
