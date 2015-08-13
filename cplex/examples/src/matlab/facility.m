function facility()
% Solve a simple version of a facility location problem.
%
% Reads data from the file ../../data/facility.dat.

% ---------------------------------------------------------------------------
% File: facility.m
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
   [capacity fixedCost cost] = inputdata ('../../data/facility.dat');
   cost = reshape (cost', 40, 1);
   
   nbLocations = length (capacity);
   nbClients   = length(cost) / nbLocations;
   
   % Build model
   cplex = Cplex('facility');
   
   % The first nbLocations variables model the decision to open a facility.
   % The next nbClients*nbLocations variables model the decision to
   % assign a client to that location.
   
   % The objective is to minimize the cost of supplying all clients.
   % There is a fixed cost for opening a each location and a cost
   % for using a particular client-location assignment.
   % All variables are binary
   cplex.Model.sense = 'minimize';
   
   obj   = [fixedCost;cost];
   lb    = zeros (nbLocations + nbLocations*nbClients, 1);
   ub    = ones (nbLocations + nbLocations*nbClients, 1);
   ctype = char (ones (1, (nbLocations + nbLocations*nbClients))*('B'));
   
   cplex.addCols(obj, [], lb, ub, ctype);
   
   % Each client must be assigned to exactly one location
   for i = 1:nbClients
      supply = zeros (1, nbLocations + nbLocations*nbClients);
      supply((i*nbLocations+1):(i*nbLocations+nbLocations)) = ...
         ones (1, nbLocations);
      cplex.addRows(1, supply, 1);
   end
   
   % The number of clients assigned to a location must be less than the
   % capacity of the location.  The capacity is multiplied by the open
   % variable to model that the capacity is zero if the location is not
   % opened.
   for i = 1:nbLocations
      v    = zeros (1, nbLocations + nbLocations*nbClients);
      v(i) = -capacity(i);
      v(i + nbLocations:nbLocations:i+nbClients*nbLocations) = ...
         ones(1, nbClients);
      cplex.addRows(-inf, v, 0);
   end
   
   cplex.solve();
   cplex.writeModel('facility.lp');
   
   % Display solution
   fprintf ('\nSolution status = %s\n',cplex.Solution.statusstring);
   if cplex.Solution.status == 101
      fprintf ('\nMinimum Cost: %f \n', cplex.Solution.objval);
      open = cplex.Solution.x(1:nbLocations);
      supply = ...
         reshape (cplex.Solution.x(nbLocations+1:end), nbLocations, nbClients);
      for i= 1:nbLocations
         if open(i) ~= 0
            fprintf ...
               ('Facility %d is open and serves the following clients:\n', i);
            for j = 1:nbClients
               if supply(i,j) >= cplex.Param.mip.tolerances.integrality.Cur
                  fprintf ('%d ',j);
               end
            end
            fprintf ('\n');
         end
      end
   end
catch m
   throw (m);
end
end
