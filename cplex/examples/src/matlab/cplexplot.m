function cplexplot(varargin)
% Plot while optimizing a MIP problem
%
% To run this example, function arguments are required,
%    i.e.,   cplexplot(filename)
% where filename is the name of the problem file to read, with .mps, .lp,
% or .sav extension.

% ---------------------------------------------------------------------------
% File: cplexplot.m
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
      display ('Usage: cplexplot(filename)');
      display ('where filename is a problem file with extension: ');
      display ('MPS, SAV, or LP (lower case is allowed)');
      display ('Exiting...')
      return;
   else
      m = varargin{1};
   end
   
   % Initialize the CPLEX object
   cplex = Cplex('cplexplot');
   cplex.Param.threads.Cur = 2;
   
   % Now read the file and copy the data into the created cpx.Model
   cplex.readModel(m);
   
   % Define InfoCallback
   cplex.InfoCallback.func = @cplexplotcb;
   cplex.InfoCallback.data.nodeLimit     = 1000;
   cplex.InfoCallback.data.acceptableGap = 10;
   
   % Declare NumIters BestObj IncObj and f and assign to empty
   NumIters = [];
   BestObj = [];
   IncObj = [];
   f = [];
   
   % Optimize the problem and obtain solution.
   cplex.solve();
   
catch m
   throw (m);
end
   function stop = cplexplotcb(info, data)
      if ~isempty (info.IncObj)
         % Get NumIters and BestObj from info
         drawnow;
         NumIters = [NumIters info.NumIters];
         BestObj = [BestObj info.BestObj];
         IncObj = [IncObj info.IncObj];
         
         % Draw figure
         if isempty (f)
            f = figure;
            set (f, 'name', ['CPLEX solving ' cplex.Model.name]);
         end
         
         % Plot
         subplot (2, 1, 1), plot (NumIters, BestObj, '--rs', 'LineWidth', 2,...
            'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'g',...
            'MarkerSize', 5);
         
         % Set xlabel, ylabel and title of plot
         xlabel ('NumIters')
         ylabel ('BestObj')
         title ('Plot of NumIters and BestObj');
         
         % Turn on grid lines for this plot
         grid on
         
         % Plot
         subplot (2, 1, 2), plot (NumIters, IncObj, '--rs', 'LineWidth', 2,...
            'MarkerEdgeColor', 'k',...
            'MarkerFaceColor', 'g',...
            'MarkerSize', 5);
         
         % Set xlabel, ylabel and title of plot
         xlabel ('NumIters')
         ylabel ('IncObj')
         title ('Plot of NumIters and IncObj');
         
         % Turn on grid lines for this plot
         grid on
         
         % If gap is good enough stop CPLEX
         gap = info.MipGap * 100.0;
         if (info.NumNodes > data.nodeLimit) && (gap < data.acceptableGap) ,
            display (sprintf ('good enouch: stopping with gap = %g', gap));
            display          '          and incumbent vector ='
            sparse (info.IncX)
            stop = true;
         else
            stop = false;
         end
      end
   end
end
