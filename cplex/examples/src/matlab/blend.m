function blend(varargin)
% Blend four sources to produce an alloy
%
% An argument may be specified: 'c' for a continuous version of the
% problem and 'i' for an integer version. 'c' is used when no argument
% is specified.
%
% The goal is to blend four sources to produce an alloy: pure metal, raw
% materials, scrap, and ingots.
%
% Ingots are discrete.  To model ingots with integer variables and solve
% as a MIP, use 'i' as an argument.  Otherwise, the problem will be
% solved as an LP and the amount of ingots may be fractional.
%
% Each source has a cost.
% Each source is made up of elements in different proportions.
% Alloy has minimum and maximum proportion of each element.
%
% Minimize the cost of producing the requested quantity of alloy.

% ---------------------------------------------------------------------------
% File: blend.m
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
   % Check for command line argument
   if nargin ~= 1
      m = 'c';
   else
      m = varargin{1};
      if m ~= 'i' &&  m ~= 'c'
         display ('Usage:  blend() or blend(X)');
         display ('where X is one of the following options: ');
         display ('  ''i''  ingots must be integral');
         display ('  ''c''  ingots can be fractional');
         display ('''c'' is assumed when there is no argument');
         display ('Exiting...');
         return;
      end
   end
   
   % Define data
   
   % Amount of alloy required
   alloy = 71;
   
   % Cost data for sources:  pure metal, raw material, scrap and ingots
   cMetals = [22.0 10.0 13.0]';
   cRaw    = [6.0 5.0]';
   cScrap  = [7.0 8.0]';
   cIngots = 9.0;
   
   nbMetals  = length (cMetals);
   nbRaw     = length (cRaw);
   nbScrap   = length (cScrap);
   nbIngots  = length (cIngots);
   
   nbElements = nbMetals;
   
   % Range of minimum and maximum proportions for each element
   pMin = [0.05, 0.30, 0.60];
   pMax = [0.10, 0.40, 0.80];
   
   % Proportion of each element and alloy in each source
   pMetals = [1.0 0.0 0.0;
      0.0 1.0 0.0;
      0.0 0.0 1.0;
      0.0 0.0 0.0];
   pRaw    = [0.20 0.01;
      0.05 0.00;
      0.05 0.30;
      0.0  0.0];
   pScrap  = [0.00 0.01;
      0.60 0.00;
      0.40 0.70;
      0.0  0.0];
   pIngots = [0.10 ;
      0.45 ;
      0.45 ;
      0.0  ] ;
   
   % Used to accumulate amount of each element in alloy
   pAlloy = [-1.0  0.0  0.0;
      0.0 -1.0  0.0;
      0.0  0.0 -1.0;
      1.0  1.0  1.0];
   
   % Lower and upper bounds on sources
   lbMetals = ones (nbMetals, 1) * 0.0;
   ubMetals = ones (nbMetals, 1) * Inf;
   
   lbRaw = ones (nbRaw, 1) * 0.0;
   ubRaw = ones (nbRaw, 1) * Inf;
   
   lbScrap = ones (nbScrap, 1) * 0.0;
   ubScrap = ones (nbScrap, 1) * Inf;
   
   lbIngots = ones (nbIngots, 1) * 0.0;
   ubIngots = ones (nbIngots, 1) * 100000;
   
   % Variables types used for the MIP formulation
   if m == 'i',
      tMetals   = char (ones ([1 nbMetals]) * 'C');
      tRaw      = char (ones ([1 nbRaw]) * 'C');
      tScrap    = char (ones ([1 nbScrap]) *'C');
      ict       = char (ones ([1 nbIngots]) *'I');
      tElements = char (ones ([1 nbElements]) *'C');
   end
   
   % Lower and upper bounds on the amount of each element in the
   % alloy are computed from the  min and max proportion of the
   % element in the alloy and the amount of alloy required
   lbElements = ones (nbElements, 1);
   ubElements = ones (nbElements, 1);
   for i = 1:nbElements
      lbElements(i) = alloy * pMin(i);
      ubElements(i) = alloy * pMax(i);
   end
   
   % Build model
   cplex = Cplex('blend');
   
   % Turn off CPLEX logging
   cplex.DisplayFunc = [];
   
   % Variables: amount of each source material used
   %            amount of each element in the final alloy
   
   % Objective Function: Minimize Cost
   cplex.Model.obj   = [cMetals; cRaw; cScrap; cIngots; zeros(nbElements, 1)];
   cplex.Model.sense = 'minimize';
   
   % Min and max quantity of each element in alloy
   cplex.Model.lb = [lbMetals; lbRaw; lbScrap; lbIngots; lbElements];
   cplex.Model.ub = [ubMetals; ubRaw; ubScrap; ubIngots; ubElements];
   
   % Specify variable types when integrality needed
   if m == 'i'
      cplex.Model.ctype = strcat (tMetals, tRaw, tScrap, ict, tElements);
   end
   
   % Constraints:  Material balance, sum of amount used is amount produced
   cplex.Model.A   = [pMetals pRaw pScrap pIngots pAlloy];
   cplex.Model.lhs = zeros (nbElements, 1);
   cplex.Model.rhs = zeros (nbElements, 1);
   
   % Constraint:  Produce requested quantity of alloy
   cplex.Model.lhs(end+1) = alloy;
   cplex.Model.rhs(end+1) = alloy;
   
   cplex.solve();
   
   % Display solution
   fprintf ('\nSolution status = %s\n',cplex.Solution.statusstring);
   fprintf ('Cost: %0.6f\n', cplex.Solution.objval);
   
   stMetals  = 1;
   endMetals = stMetals+nbMetals-1;
   disp ('Pure metal:');
   disp (cplex.Solution.x(stMetals:endMetals)');
   
   stRaw  = endMetals+1;
   endRaw = stRaw+nbRaw-1;
   disp ('Raw material:');
   disp (cplex.Solution.x(stRaw:endRaw)');
   
   stScrap  = endRaw+1;
   endScrap = stScrap+nbScrap-1;
   disp ('Scrap:');
   disp (cplex.Solution.x(stScrap:endScrap)');
   
   stIngots  = endScrap+1;
   endIngots = stIngots+nbIngots-1;
   disp ('Ingots:');
   disp (cplex.Solution.x(stIngots:endIngots)');
   
   stElements  = endIngots+1;
   endElements = stElements+nbElements-1;
   disp ('Elements in Alloy:');
   disp (cplex.Solution.x(stElements:endElements)');
catch m
   throw (m);
   clear all;
   clear classes;
end
% Example output for the continuous variant
%
% Cost:
%   653.5541
%
% Pure metal:
%      0
%      0
%      0
%
% Raw material:
%      0
%      0
%
% Scrap:
%    17.0590
%    30.2311
%
% Ingots:
%    32.4769
%
% Elements in Alloy:
%     3.5500
%    24.8500
%    42.6000



