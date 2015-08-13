function cutstock(varargin)
% Use column generation to solve a cutting-stock problem
%
% The cutting stock problem in this example is sometimes known in math
% programming terms as a knapsack problem with reduced costs in the
% objective function. Generally, a cutting stock problem begins with a
% supply of rolls of material of fixed length (the stock). Strips are cut
% from these rolls. All the strips cut from one roll are known together as
% a pattern. The point of this example is to use as few rolls of stock as
% possible to satisfy some specified demand of strips. By convention, it is
% assumed that only one pattern is laid out across the stock; consequently,
% only one dimension (the width) of each roll of stock is important.

% ---------------------------------------------------------------------------
% File: cutstock.m
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
   % Read the data
   [rollWidth, size, amount] = inputdata ('../../data/cutstock.dat');
   
   % Build Model
   % Cutting-optimization problem
   RC_EPS = 1.0e-6;
   cutOpt = Cplex();
   cutOpt.Model.sense = 'minimize';
   nWdth = length (size);
   lhs = amount;
   rhs = ones (nWdth, 1) * inf;
   cutOpt.addRows(lhs, [], rhs);
   for j = 1:nWdth
      fill = zeros(nWdth,1);
      fill(j) = floor ((rollWidth/size(j)));
      cutOpt.addCols(1, fill, 0, inf);
   end
   cutOpt.solve();
   
   report1(cutOpt);
   
   % Pattern-generation problem
   patGen = Cplex();
   patGen.Model.sense = 'minimize';
   lhs = 0;
   rhs = rollWidth;
   patGen.addRows(lhs, [], rhs);
   
   % Column-generation procedure
   i = 1;
   while  i
      cutOpt.solve();
      report1 (cutOpt);
      obj = [1; -cutOpt.Solution.dual];
      if i == 1
         patGen.addCols(obj, [0 size'], [1;zeros(nWdth,1)], [1;ones(nWdth,1)...
            *inf]);
      else
         patGen.Model.obj = obj;
      end
      patGen.Model.ctype = char (ones ([1, length(patGen.Model.obj)]) * ('I'));
      patGen.solve();
      report2 (patGen);
      if  patGen.Solution.objval > -RC_EPS,
         break,
      end
      cutOpt.addCols(1, newPatt, 0, inf);
      i = i+1;
   end
   
   cutOpt.Model.ctype = char (ones ([1, length(cutOpt.Model.obj)]) * ('I'));
   cutOpt.solve();
   fprintf ('\nSolution status = %s\n',cutOpt.Solution.statusstring);
   report3 (cutOpt);
catch m
   throw (m);
end
   function report1(cutOpt)
      fprintf ('Using %f rolls\n', cutOpt.Solution.objval);
      disp ('Cut=');
      disp (cutOpt.Solution.x');
      disp ('Fill=');
      disp (cutOpt.Solution.dual');
   end
   function report2(patGen)
      fprintf ('Reduced cost is %f \n', patGen.Solution.objval);
      disp ('Use=');
      disp (patGen.Solution.x(2:(nWdth+1))');
      newPatt = patGen.Solution.x(2:(nWdth+1));
   end
   function report3(cutOpt)
      fprintf ('Best integer Solution uses %f rolls\n', ...
               cutOpt.Solution.objval);
      disp ('cut=');
      disp (cutOpt.Solution.x');
   end
end
