function transport(varargin)
% Model piecewise linear cost coefficients with special ordered sets (SOS)
%
% The problem is a simple transportation model. Set the function argument
% to 0 for a convex piecewise linear model and to 1 for a concave
% piecewise linear model.
%

% ---------------------------------------------------------------------------
% File: transport.m
% Version 12.6.1
% ---------------------------------------------------------------------------
% Licensed Materials - Property of IBM
% 5725-A06 5725-A29 5724-Y48 5724-Y49 5724-Y54 5724-Y55 5655-Y21
% Copyright IBM Corporation 2008, 2014. All Rights Reserved.
%
% US Government Users Restricted Rights - Use, duplication or
% disclosure restricted by GSA ADP Schedule Contract with IBM Corp.
% ---------------------------------------------------------------------------

% Define data
if nargin == 0
    disp ('Specify an argument to choose between convex and concave problems');
    disp ('Usage: transport model')
    disp (' model = 0 -> convex piecewise linear model');
    disp (' model = 1 -> concave piecewise linear model. [default]');
end

if nargin ~=  1
    convex = 0;
else
    if (varargin{1} == num2str(0))
        convex = 1;
    else
        convex = 0;
    end
end

try
    supply   = [1000.0 850.0 1250.0]';
    nbSupply = length (supply);
    demand   = [ 900.0 1200.0 600.0 400]';
    nbDemand = length (demand);
    n = nbSupply*nbDemand;
    
    % The x coordinate of the last break point of pwl
    k = 1;
    pwl_x = zeros(n, 4);
    pwl_y = zeros(n, 4);

    for i = 1:nbSupply
        for j = 1:nbDemand
            if supply(i) < demand(j)
                midval = supply(i);
            else
                midval = demand(j);
            end
            
            pwl_x(k,:) = [0  200 400 midval];
            
            if convex == 1
                pwl_slope = [30 80 130];
            else
                pwl_slope = [120 80 50];
            end
            
            pwl_y(k, 1) = 0;
            pwl_y(k, 2) = 200 * pwl_slope(1);
            pwl_y(k, 3) = pwl_y(k, 2) + pwl_slope(2) * (400 - 200);
            pwl_y(k, 4) = pwl_y(k, 3) + pwl_slope(3) * (midval - 400);
            
            k = k + 1;
        end
    end
    
    % Build model
    cplex = Cplex('transport');
    cplex.Model.sense = 'minimize';
    
    % x(varind(i, j)) is the amount that is shipped from supplier i to recipient j.
    obj = zeros (n, 1);
    lb  = zeros (n, 1);
    ub  = ones (n, 1) * inf;
    cplex.addCols(obj, [], lb, ub, []);
    
    % y(varind(i, j)) is used to model the PWL cost associated with this shipment.
    colname_y = cell ([n, 1]);
    for k = 1:n
        colname_y{k} = ['y' num2str(k)];
    end
    colname_y = char (colname_y);
    cplex.addCols(ones(n, 1), [], lb, ub, [], colname_y);
    
    % Generate colnames for lambdas.
    colname_id = cell([4*(nbDemand)*nbSupply, 1]);
    for k = 1:4*(nbDemand)*nbSupply
        colname_id{k} = ['lambda_' num2str(k)];
    end
    colname_id = char (colname_id);

    % Add columns for lambda.
    cplex.addCols(zeros(4*(nbDemand)*nbSupply, 1), [], [], [], [], colname_id);
    
    ncols = length(cplex.Model.obj);
    
    % Supply must meet demand.
    for i = 1:nbSupply
        A1 = zeros (1, ncols);
        for j = 1:nbDemand
            A1(varindex (i, j)) = 1;
        end
        cplex.addRows(supply(i), A1, supply(i));
    end
    
    % Demand must meet supply
    for i = 1:nbDemand
        A2 = zeros (1, ncols);
        for j = 1:nbSupply
            A2(varindex (j, i)) = 1;
        end
        cplex.addRows(demand(i), A2, demand(i));
    end
    
    % Add constrains about lambda
    % x = SUM(lambda_i*x_i) 
    % y = SUM(lambda_i*y_i)
    % SUM(lambda_i * x_i) = 1
    for i = 0:n-1
        index  = 2*n + i*4 + 1: 2*n + i*4 + 4;
        
        A = zeros (1, ncols);
        A(index) = pwl_x(i + 1, :);
        A(i+1) = -1;
        cplex.addRows(0, A, 0);
        
        A = zeros (1, ncols);
        A(index) = pwl_y(i + 1, :);
        A(n+i+1) = -1;
        cplex.addRows(0, A, 0);
        
        % SUM(lambda_i * x_i) = 1
        A = zeros (1, ncols);
        A(index) = ones (1, 4);
        cplex.addRows(1, A, 1);
        
        cplex.addSOSs('2', index', [pwl_x(i + 1, 1:3) pwl_x(i + 1, 4) + 1]');
    end
    
    % Solve model
    cplex.solve();
    cplex.writeModel('transport.lp');
    
    % Display solution
    fprintf ('\nSolution status = %s\n', cplex.Solution.statusstring);
    disp (' - Solution:');
    for i = 1:nbSupply
       fprintf ('\n   %d : ', i);
       for j = 1:nbDemand
         fprintf('%d\t', cplex.Solution.x(varindex(i, j)));
       end
    end
    fprintf    ('\n   Cost = %f\n', cplex.Solution.objval);

catch m
    throw (m);
end
    function out = varindex(m, n)
        out = (m - 1) * nbDemand + n;
    end
end
