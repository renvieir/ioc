function varargout = inputdata(varargin)
% inputdata function is used to input data from .dat file to example.
% see the usage in diet.m

% ---------------------------------------------------------------------------
% File: inputdata.m
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
      display ('Usage: inputdata(filename)');
      display ('where filename is a problem file with extension: .dat');
      display ('Exiting...')
      return;
   else
      m = varargin{1};
   end
   
   k = 1;
   fid = fopen (m);
   while(1)
      % Read directly number or array or matrix from fid
      varargout{k} = fscanf (fid, '%lg',1);
      if isempty (varargout{k})
         varargout{k} = readarray (fid);
         if isempty (varargout{k})
            varargout{k} = readmatrix (fid);
            if isempty (varargout{k})
               break;
            end
         end
      else
         k = k + 1;
         continue;
      end
      k = k + 1;
   end
   fclose (fid);
catch m
   fclose (fid);
   throw (m);
end
   function array = readarray(in)
      % Read array from fid which returned by fopen function.
      num = 1;
      i = 1;
      while(1)
         ch = fscanf (in, '%c',1);
         if isempty (ch)
            array = [];
            return;
         end
         if strcmp(ch, '\t') ||...
            strcmp(ch, '\r') ||...
            strcmp(ch, ' ')  ||...
            strcmp(ch, '\n'),
            i = i + 1;
            continue;
         end
         if  ch == '['
            i = i + 1;
            break;
         end
      end
      
      while(1)
         item = fscanf (in, '%lg',num);
         if isempty (item)
            array = [];
            return;
         end
         array(num) = item;
         num = num + 1;
         
         ch = fscanf (in, '%c', 1);
         if strcmp (ch, '\t') || ...
            strcmp (ch, '\r') || ...
            strcmp (ch, ' ')  || ...
            strcmp (ch, '\n'),
            while (1)
               ch = fscanf (in, '%c', 1);
               if strcmp (ch, '\t') || ...
                  strcmp (ch, '\r') || ...
                  strcmp (ch, ' ')  || ...
                  strcmp (ch, '\n'),
                  i = i + 1;s
                  continue;
               else
                  i = i + 1;
                  break;
               end
            end
         end
         
         if  ch == ']'
            break;
         elseif  ch ~= ','
            array = [];
            error ('item in array should be seprated by "," ')
         end
         i = i + 1;
         
      end
      array = array';
   end
   function matrix = readmatrix(in)
      % Read matrix from fid which returned by fopen function.
      j = 1;
      while(1)
         array = readarray (in);
         if isempty (array)
            matrix = [];
            return;
         end
         matrix(j,:) = array;
         j = j + 1;
         ch = fscanf (in, '%c', 1);
         if strcmp (ch, '\t') ||...
            strcmp (ch, '\r') ||...
            strcmp (ch, ' ')  ||...
            strcmp (ch, '\n')
            continue;
         end
         if  ch == ']'
            break;
         end
         if  ch ~= ','
            error ('array in matrix should be seprated by ","');
            matrix = [];
         end
      end
   end
end

