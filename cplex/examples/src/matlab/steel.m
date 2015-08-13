function  steel()
% Solve a multi-period production planning problem
%  This example is an implementation of the model called "steelT.mod"
%  on page 58 in the AMPL book by Fourer, Gay and Kernighan.  In the AMPL
%  example, a multiperiod production model is given, with data
%  for 4 weeks.
%
%  The parameters for the model are:
%   nProd              the number of products
%   nTime              the number of time periods
%   rate(p)            rate of production for product p
%   inv0(p)            initial inventoryfor product p
%   avail(t)           hours available in time period t
%   market(p)(t)       demand for product p in time period t
%   prodcost(p)        production cost per unit of product p
%   invcost(p)         inventory cost per unit of product p
%   revenue(p)(t)      revenue per unit of product p in time period t
%
%  The decision variables of the model are:
%    Make(p)(t)         amount produced of product p in time period t
%    Inv(p)(t)          amount inventoried of product p in time period t
%    Sell(p)(t)         amount sold of product p in time period t
%  The objective function is to
%  maximize  sum(over p,t) (  revenue(p)(t) * Sell(p)(t)
%                           - prodcost(p)   * Make(p)(t)
%                           - invcost(p)    * Inv(p)(t)  )
%
%  The constraints are
%
%   For each t:   (time availability constraint)
%       sum(over p)  ( (1/rate(p)) * Make(p)(t) ) <= avail(t)
%   For each p, (t=0): (balance constraint)
%       Make(p)(0) - Sell(p)(0) - Inv(p)(0) = -inv0(p)
%   For each pair (p,t) (t>0): (balance constraint)
%
%   The bounds on the variables are:
%     All variables are nonnegative ( >= 0 )
%     For each (p,t),
%        Sell(p)(t) <= market(p)(t)
%     All other variables have infinite upper bounds.

% ---------------------------------------------------------------------------
% File: steel.m
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
   % Define data
   if nargin>1
      [avail rate inv0 prodCost invCost revenue market] = ...
         inputdata (varargin{1});
   else
      [avail rate inv0 prodCost invCost revenue market] = ...
         inputdata ('../../data/steel.dat');
   end
   
   revenue = reshape (revenue', [], 1);
   market = reshape (market', [], 1);
   
   nTime=length (avail);
   nProd=length (rate);
   % Build model
   cplex=Cplex('steel');
   cplex.Model.sense = 'maximize';
   prcost = zeros (nProd * nTime, 1);
   for i = 1:nProd
      prcost(((i-1)*nTime+1):((i-1)*nTime+nTime)) = prodCost(i)*ones(nTime,1);
   end
   incost = zeros (nProd*nTime, 1);
   for i = 1:nProd
      incost(((i-1)*nTime+1):((i-1)*nTime+nTime)) = invCost(i)*ones(nTime,1);
   end
   obj = [revenue; -prcost; -incost];
   lb = zeros (length(obj), 1);
   ub = [market;ones((length(obj)-nProd*nTime),1)*inf];
   cplex.addCols(obj, [], lb, ub);
   
   % Add time availability constraint
   AA = [];
   lhs = [];
   for i = 1:nTime
      A = zeros (1, length (cplex.Model.obj));
      A(nProd*nTime+i) = 1/rate(1);
      A(nProd*nTime+nTime+i) = 1/rate(2);
      AA = [AA; A];
      lhs = [lhs; -inf];
   end
   cplex.addRows(lhs, AA, avail);
   
   % Add balance constraint
   AA = [];
   for i = 1:nProd
      A = zeros (1, length (cplex.Model.obj));
      A(i) = -1;
      A(i+nProd*nTime) = 1;
      A(i+2*nProd*nTime+1) = -1;
      AA = [AA;A];
   end
   cplex.addRows(-inv0, AA, -inv0);
   AA = [];
   rhs = zeros (nProd*nTime-1, 1);
   lhs = zeros (nProd*nTime-1, 1);
   for i = 2:nProd*nTime
      A= zeros (1, length (cplex.Model.obj));
      A(i) = -1;
      A(i+nProd*nTime) =1;
      A(i+2*nProd*nTime-1) =1;
      A(i+2*nProd*nTime) =-1;
      AA = [AA;A];
   end
   cplex.addRows(lhs,AA,rhs);
   
   % Solve model
   cplex.solve();
   
   % Display solution
   fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
   fprintf ('\nprofit: %.1f\n',cplex.Solution.objval);
   
   for i=1:nProd
      fprintf ('product %d amount produced : \n',i);
      fprintf ('%.1f ,', ...
               cplex.Solution.x(((i-1)*nTime+1):((i-1)*nTime+nTime))');
      fprintf ('\namount sold : \n');
      fprintf ('%.1f ,', ...
               cplex.Solution.x(((i+1)*nTime+1):((i+1)*nTime+nTime))') ;
      fprintf ('\namount inventoried : \n');
      fprintf ('%.1f ,', ...
               cplex.Solution.x(((i+3)*nTime+1):((i+3)*nTime+nTime))') ;
      fprintf ('\n');
   end
catch m
   throw (m);
end
end




