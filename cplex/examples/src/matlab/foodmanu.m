function foodmanu
% A formulation of the food manufacturing problem using indicator
% constraints and semi-continuous variables 

% ---------------------------------------------------------------------------
% File: foodmanu.m
% Version 12.6.1
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
% Copyright IBM Corporation 2008, 2014. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
% ---------------------------------------------------------------------------

NUMPRODUCTS = 5;
VEGOIL1  = 0;
VEGOIL2  = 1;
OIL1     = 2;
OIL2     = 3;
OIL3     = 4;

NUMMONTHS   = 6;

NUMVARS     = 4;
BUY         = 0;
USE         = 1;
STORE       = 2;
IS_USED     = 3;

cost = [110.0, 120.0, 130.0, 110.0, 115.0, ...% Cost for January
   130.0, 130.0, 110.0,  90.0, 115.0, ...% Cost for February
   110.0, 140.0, 130.0, 100.0,  95.0, ...% Cost for March
   120.0, 110.0, 120.0, 120.0, 125.0, ...% Cost for April
   100.0, 120.0, 150.0, 110.0, 105.0, ...% Cost for May
   90.0, 100.0, 140.0,  80.0, 135.0];    % Cost for July

hardness = [8.8, 6.1, 2.0, 4.2, 5.0]; % Hardness of each product

% Initialize the CPLEX environment
cplex = Cplex('foodmanufact');

% Build the model
buildmodel;

% Write a copy of the problem to a file.
cplex.writeModel('foodmanu.lp');

% Optimize the problem and obtain solution.
cplex.solve();

% Display solution
fprintf ('\nSolution status = %s\n',cplex.Solution.statusstring);
fprintf('\n Maximum profit: %f\n',cplex.Solution.objval);
for i = 0:NUMMONTHS-1
   fprintf('\nMonth%d ', i);
   fprintf('\n. buy : ');
   for j = 0:NUMPRODUCTS-1
      fprintf ('%f\t', cplex.Solution.x(varindex(i, j, BUY)));
   end
   fprintf ('\n');
   
   fprintf('\n. use : ');
   for j = 0:NUMPRODUCTS-1
      fprintf ('%f\t', cplex.Solution.x(varindex(i, j, USE)));
   end
   fprintf ('\n');
   
   fprintf('\n. store : ');
   for j = 0:NUMPRODUCTS-1
      fprintf ('%f\t', cplex.Solution.x(varindex(i, j, STORE)));
   end
   fprintf ('\n');
end

disp (' ');

   function out = varindex(m, p, whichvar)
      % A function to map month, product and which variable to an index
      out = m*NUMVARS*NUMPRODUCTS + p*NUMVARS + whichvar + 1;
   end

   function buildmodel
      % Build food manufacturing model using indicator constraints and
      % semi-continuous variables
      cplex.Model.sense = 'maximize';
      colcnt = NUMVARS*NUMMONTHS*NUMPRODUCTS;
      obj     = zeros (colcnt, 1);
      lb      = zeros (colcnt, 1);
      ub      = zeros (colcnt, 1);
      ctype   = char(ones(1, colcnt)* 'C');
      
      % Create variables. For each month and each product, we have 3
      % variables corresponding to the quantity used (semi-continuous),
      % stored (continuous) and bought (continuous) and one binary
      % variable indicating whether or not the product is used during
      % this month.
      for m  = 0:NUMMONTHS-1
         for p = 0 :NUMPRODUCTS-1
            % The quantity bought is a continuous variable. It has a cost
            obj(varindex(m, p, BUY)) = -cost(m*NUMPRODUCTS + p + 1);
            lb(varindex(m, p, BUY))  = 0.0;
            ub(varindex(m, p, BUY))  = Inf;
            ctype(varindex(m, p, BUY))  = 'C';
            % When an oil is used, the quantity must be at least 20
            % tons. This is modeled as a semi-continuous variable.
            obj(varindex(m, p, USE)) = 0.0;
            lb(varindex(m, p, USE))  = 20.0;
            ub(varindex(m, p, USE))  = Inf;
            ctype(varindex(m, p, USE))  = 'S';
            % It is possible to store up to 1000 tons of each
            % product. There are storage costs.
            obj(varindex(m, p, STORE)) = -5.0;
            lb(varindex(m, p, STORE))  = 0.0;
            ub(varindex(m, p, STORE))  = 1000.0;
            ctype(varindex(m, p, STORE))  = 'C';
            % At the end, we must have exactly 500 tons of each
            % product in storage.
            if m == NUMMONTHS - 1
               lb(varindex (m, p, STORE)) = 500.0;
               ub(varindex (m, p, STORE)) = 500.0;
            end
            % The variable indicating whether or not a product is
            % used during a month is a binary variable.
            obj(varindex (m, p, IS_USED))   = 0.0;
            lb(varindex (m, p, IS_USED))    = 0.0;
            ub(varindex (m, p, IS_USED))    = 1.0;
            ctype(varindex (m, p, IS_USED))    = 'B';
         end
      end
      cplex.addCols(obj, [], lb, ub, ctype);
      
      % Constraints for each month
      for m  = 0:NUMMONTHS-1
         for p = 0:NUMPRODUCTS-1
            % For each product, create an indicator constraint linking the
            % quantity used and the binary variable indicating whether or
            % not the product is used
            a = sparse(length(cplex.Model.obj),1);
            a(varindex(m, p, USE)) = 1.0;
            rhs         = 0.0;
            sense        = 'L';
            indvar       = varindex(m, p, IS_USED);
            complemented = 1;
            cplex.addIndicators(indvar, complemented, a, sense, rhs);
         end
         
         % Not more than 200 tons of vegetable oil can be refined
         A1 = sparse(1, length(cplex.Model.obj));
         A1(varindex (m, VEGOIL1, USE)) = 1.0;
         A1(varindex (m, VEGOIL2, USE)) = 1.0;
         cplex.addRows(-inf, A1, 200.0);
         
         % Not more than 250 tons of non-vegetable oil can be refined
         A2 = sparse(1, length(cplex.Model.obj));
         A2(varindex (m, OIL1, USE)) = 1.0;
         A2(varindex (m, OIL2, USE)) = 1.0;
         A2(varindex (m, OIL3, USE)) = 1.0;
         cplex.addRows(-inf, A2, 250.0);
         
         % Constraint on food composition
         % Add a variable corresponding to total quantity produced
         % in a month
         cplex.addCols(150.0, [], 0.0, Inf, 'C');
         totalindex = length(cplex.Model.obj);
         
         % Total quantity = sum (quantities)
         A3 = sparse(1, length(cplex.Model.obj));
         for p = 0:NUMPRODUCTS-1
            A3(varindex (m, p, USE)) = 1.0;
         end
         A3(totalindex) = -1.0;
         cplex.addRows(0, A3, 0);
         
         % Hardness constraints
         % sum (quantity * hardness) >= 3 * total quantity
         % sum (quantity * hardness) <= 6 * total quantity
         A4 = sparse(1, length(cplex.Model.obj));
         for p = 0:NUMPRODUCTS-1
            A4(varindex (m, p, USE)) = hardness(p+1);
         end
         A4(totalindex) = -3.0;
         cplex.addRows(0, A4, Inf);
         
         A5 = sparse(1, length(cplex.Model.obj));
         for p = 0:NUMPRODUCTS-1
            A5(varindex (m, p, USE)) = hardness(p+1);
         end
         A5(totalindex) = -6.0;
         cplex.addRows(-inf, A5, 0);
         
         % The food may never be made up of more than three oils
         A6 = sparse(1, length(cplex.Model.obj));
         for p = 0:NUMPRODUCTS-1
            A6(varindex (m, p, IS_USED)) = 1.0;
         end
         cplex.addRows(-inf, A6, 3);
         
         % If product veg 1 or veg 2 is used then oil 3 must be used
         a = sparse(length(cplex.Model.obj),1);
         a(varindex (m, OIL3, USE)) = 1.0;
         rhs         = 20;
         sense        = 'G';
         indvar       = varindex (m, VEGOIL1, IS_USED);
         complemented = 0;
         cplex.addIndicators(indvar, complemented, a, sense, rhs);
         indvar       = varindex (m, VEGOIL2, IS_USED);
         cplex.addIndicators(indvar, complemented, a, sense, rhs);
         
         % We can store each product from one month to the next,
         % starting with a stock of 500 tons
         for p = 0:NUMPRODUCTS-1
            A7 = sparse(1, length(cplex.Model.obj));
            if m ~= 0
               A7(varindex (m-1, p, STORE)) = 1.0;
               A7(varindex (m, p, BUY)) = 1.0;
               A7(varindex (m, p, USE)) = -1.0;
               A7(varindex (m, p, STORE)) = -1.0;
               lhs = 0;
               rhs = 0;
            else
               A7(varindex (m, p, BUY)) = 1.0;
               A7(varindex (m, p, USE)) = -1.0;
               A7(varindex (m, p, STORE)) = -1.0;
               lhs = -500.0;
               rhs = -500.0;
            end
            cplex.addRows(lhs, A7, rhs);
         end
      end
   end
end
