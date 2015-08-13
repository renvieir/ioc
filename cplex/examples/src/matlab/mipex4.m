function mipex4(varargin)
%  Read and optimize a MIP problem using a callback

% ---------------------------------------------------------------------------
% File: mipex4.m
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
   % Check the command line arguments
   if nargin ~= 1
      display ('Usage: mipex4(filename)');
      display ('where filename is a file with extension: ');
      display ('MPS, SAV, or LP (lower case is allowed)');
      display ('Exiting...')
      return;
   else
      file = varargin{1};
   end
   % Initialize the CPLEX object
   cplex = Cplex('mipex4');
   
   % Now read the file and copy the data into cplex.Model
   cplex.readModel(file);
   
   % Set InfoCallback
   cplex.InfoCallback.func = @mipex4cb;
   cplex.InfoCallback.data.nodeLimit     = 1000;
   cplex.InfoCallback.data.acceptableGap = 10;
   
   % Optimize the problem
   cplex.solve();
   fprintf ('\nSolution status = %d\n', cplex.Solution.status);
catch m
   throw (m);
end

% Callback function
   function stop = mipex4cb(info,data)
      if ~isempty (info.IncObj),
         gap = info.MipGap * 100.0;
         if (info.NumNodes > data.nodeLimit) && ...
            (gap < data.acceptableGap) ,
            display (sprintf ('good enough: stopping with gap = %g', gap));
            display          '          and incumbent vector ='
            sparse (info.IncX)
            stop = true;
         else
            stop = false;
         end
      end
   end
end



