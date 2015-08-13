function distmipex2(varargin)
% Read and optimize a mixed integer programming problem while
% monitoring bounds with an info callback
%
% To run this example, a function argument is required.
%    i.e.,   distmipex2(vmconfig, filename)
%    where
%       filename is the name of the problem file, with .mps, .lp, or .sav
%       extension
%    Example:
%       distmipex2('process.vmc', 'mexample.mps')

% ---------------------------------------------------------------------------
% File: distmipex2.m
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
   if nargin ~= 2
      display ('Usage: distmipex2(vmconfig, filename)');
      display ('where filename is a file with extension: ');
      display ('MPS, SAV, or LP (lower case is allowed)');
      display ('Exiting...')
      return;
      else
      vmconfig = varargin{1};
      model    = varargin{2};
   end

   % Initialize the CPLEX object
   cpx = Cplex('distmipex2');

   % Now read the file and copy the data into cpx.Model
   cpx.readModel(model);

   % Now read and setup the virtual machine configuration
   cpx.readVMConfig(vmconfig);
   
   % Set InfoCallback
   cpx.InfoCallback.func = @distmipex2cb;
   cpx.InfoCallback.data.acceptableGap = 10;

   % Optimize the problem
   cpx.solve();

   % Write the solution
   fprintf ('\nSolution status = %s \n', cpx.Solution.statusstring);
   fprintf ('Solution value = %f \n', cpx.Solution.objval);
   disp ('Values =');
catch m
   throw (m);
end

% Callback function
   function stop = distmipex2cb(info,data)
      if ~isempty (info.IncObj),
         gap = info.MipGap * 100.0;
         if (gap < data.acceptableGap) 
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



