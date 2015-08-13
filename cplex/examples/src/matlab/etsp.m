function etsp(varargin)
% Model an earliness-tardiness scheduling problem with indicator constraints
%

% ---------------------------------------------------------------------------
% File: etsp.m
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
   if nargin > 1
      [activityOnAResource, duration, jobDueDate, jobEarlinessCost, ...
         jobTardinessCost] = inputdata (varargin{1});
   else
      [activityOnAResource, duration, jobDueDate, jobEarlinessCost, ...
         jobTardinessCost] = inputdata ('../../data/etsp.dat');
   end
   
   nbJob = length (jobDueDate);
   nbResource = size (activityOnAResource, 2);
   
   % Build model
   cplex = Cplex();
   cplex.Model.sense = 'minimize';
   
   % Add activity start time variables
   nbStart = nbJob * nbResource;
   lb      = zeros (nbStart, 1);
   ub      = ones (nbStart, 1)*10000;
   cplex.addCols(zeros (nbStart, 1), [], lb, ub);
   
   % State precedence constraints
   % starttime(i, j) - starttime(i, j-1) >= duration(i, j-1)
   for i = 1:nbJob
      for j = 2:nbResource
         A = sparse ...
            ([1 1], [starttime(i, j) starttime(i, j-1)], [1 -1], 1, nbStart);
         lhs = duration(i, j-1);
         cplex.addRows(lhs, A, Inf);
      end
   end
   
   % Add indicator variables
   nIndicatorVars = nbResource * nbJob * (nbJob - 1);
   colname_ind = cell ([nIndicatorVars, 1]);
   for k = 1:nIndicatorVars
      colname_ind{k} = ['ind' num2str(k)];
   end
   colname_ind = char (colname_ind);
   cplex.addCols(zeros (nIndicatorVars, 1), [], [], [], ...
                 char (ones (1, 448) * 'B'), colname_ind);
              
   % Add ind1 + ind2 >= 1
   %     ind3 + ind4 >= 1
   %     ind5 + ind6 >= 1
   %     ...
   % constrains
   lhs = ones(nIndicatorVars/2, 1);
   rhs = ones(nIndicatorVars/2, 1)*inf;
   AA = zeros(nIndicatorVars/2, length(cplex.Model.obj));
   j =  nbStart + 1;
   for i = 1:nIndicatorVars/2
      AA(i, :) = sparse ([1 1], [j j+1], [1 1], 1, length (cplex.Model.obj));
      j = j + 2;
   end
   
   cplex.addRows(lhs, AA, rhs);
   
   % Add indicator constrains
   % i1 = 1 <-> c1
   % i2 = 1 <-> c2
   index = 1;
   for i = 1:nbResource
      e = nbJob - 1;
      for j = 1:e
         activity1 = activityOnAResource(i, j) + 1;
         for k = j + 1:nbJob
            activity2 = activityOnAResource(i, k) + 1;
            variable = nbStart + index;
            a = sparse(length(cplex.Model.obj), 1);
            a(starttime(j, activity1)) = 1;
            a(starttime(k, activity2)) = -1;
            rhs = duration(k, activity2);

            % ind(nbStart + index) = 1 -> ...
            % starttime(j, activity1) - starttime(k, activity2) >=
            % duration(k, activity2)
            cplex.addIndicators(variable,...
               0,...
               a,...
               'G',...
               rhs);

            % ind(nbStart + index) = 0 -> ...
            % starttime(j, activity1) - starttime(k, activity2) <=
            % duration(k, activity2)
            cplex.addIndicators(variable,...
               1,...
               a,...
               'L',...
               rhs);
            
            variable = nbStart + index + 1;
            a = sparse(length(cplex.Model.obj), 1);
            a(starttime(k, activity2)) = 1;
            a(starttime(j, activity1)) = -1;
            rhs = duration(j, activity1);

            % ind(nbStart + index + 1)  = 1 -> ...
            % starttime(k, activity2) - starttime(j, activity1) >=
            % duration(j, activity1)
            cplex.addIndicators(variable,...
               0,...
               a,...
               'G',...
               rhs);

           % ind(nbStart + index + 1)  = 0 -> ...
           % starttime(k, activity2) - starttime(j, activity1) <=
           % duration(j, activity1)
            cplex.addIndicators(variable,...
               1,...
               a,...
               'L',...
               rhs);
            
            index = index + 2;

         end
      end
   end
   
   % Add Objective function
   % each job have a cost which contain jobEarlinessCost and
   % jobTardinessCost, and get the index of Earliness variables,
   % Tardiness variables and Endness variables
   indexOfEarlinessVar = ...
       length(cplex.Model.obj) + 1:length(cplex.Model.obj) + nbJob;
   cplex.addCols(jobEarlinessCost, [], [], [], []);
   
   indexOfTardinessVar = indexOfEarlinessVar + nbJob;
   cplex.addCols(jobTardinessCost, [], [], [], []);
   
   % Add finished time variable
   indexOfEndnessVar = indexOfTardinessVar + nbJob;
   cplex.addCols(zeros(nbJob, 1));
   
   % Add constrains for each Job
   % indexOfEndnessVar(i) - starttime(i, nbResource) = duration(i, nbResource)
   A = sparse(nbJob, length(cplex.Model.obj));
   lhs = ones(nbJob, 1);
   for i = 1:nbJob
      A(i, indexOfEndnessVar(i)) = 1;
      A(i, starttime(i, nbResource)) = -1;
      lhs(i) = duration(i, nbResource);
   end
   cplex.addRows(lhs, A, lhs);

   % Add constrains for each Job
   % jobDueDate(i) = ...
   % indexOfEndnessVar(i) + indexOfEarlinessVar(i) - indexOfTardinessVar(i)
   A = sparse(nbJob, length(cplex.Model.obj));
   lhs = ones(nbJob, 1);
   for i = 1:nbJob
      A(i, indexOfEndnessVar(i)) = 1;
      A(i, indexOfEarlinessVar(i)) = 1;
      A(i, indexOfTardinessVar(i)) = -1;
      lhs(i) = jobDueDate(i);
   end
   cplex.addRows(lhs, A, lhs);
   
   cplex.Param.emphasis.mip.Cur = 4;
   cplex.solve();
   cplex.writeModel('etsp.sav');
   
   % Display solution
   fprintf ('\nSolution status = %s\n',cplex.Solution.statusstring);
   fprintf('Optimal Value =  %0.5f \n',cplex.Solution.objval);
catch m
   disp(m.message);
end
   % Get index of start time variables
   function t = starttime(job, res)
      t = (job - 1) * nbJob + res;
   end
end

