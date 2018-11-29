function [min_x,min_y,X,Y]=gridsearch(Min,Max,Grid,nloops,fun,varargin)
% GRIDSEARCH Function minimization using grid search.
%
% Synopsis:
%  [min_x,min_y,X,Y]=gridsearch(Min,Max,Grid,nloops,fun)
%  [min_x,min_y,X,Y]=gridsearch(Min,Max,Grid,nloops,fun,varargin)
%
% Descritpion:
%  This function implements the grid search to find the minimum
%  of a given function y = fun(x), x \in X.
%  The domain X is discretized to the grid with minimal and 
%  maximal point given by Min and Max respectively. Number of points 
%  of the grid is given by Grid. The Max, Min and Grid are vectors
%  their entries correspond to individual dimensions of X.
%  After the mininum is found then the grid search is recursively 
%  repeated with a finer grid. Number of nested loops is given by nloops. 
%  The string fun determines function to be minimized. The function is 
%   called as y=feval(fun,x,varargin{:}).
%  
% Input:
%  Min [dim x 1] Minimum point of the grid.
%  Max [dim x 1] Maximum point of the grid.
%  Grid [dim x 1] Number of point in the grid, i.e. grid density.
%  nloops [1x1] Number of nested loops of the grid search.
%  fun [string] Identifies the minimized function.
%  varargin [cell] Additional arguments of the minimized function.
%
% Output:
%  min_x [dim x 1] Found minimum.
%  min_y [1x1] min_y = fun(min_x).
%  X [dim x n] Points which the grid search checked through.
%  Y [1 x n] Y(i)=fun(X(:,i)).
%
% Example:
%  [min_x,min_y,X,Y]=gridsearch(0,10,10,3,'sin');
%  figure; grid on; hold on;
%  plot(X,Y,'.');
%  plot(min_x,min_y,'+r');
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
%  18-July-2003, VF
%  17-July-2003, VF

num_variables=length(Min);

if num_variables > length(Grid),
   Grid = Grid*ones(size(Min));
end

if Grid > 1, step=(Max-Min)./(Grid-1); else step = 0; end
x=Min;

min_y = inf;
stop = 0;

X = [];
Y = [];

while ~stop,
   
   y = feval(fun, x, varargin{:});

   X = [X,x(:)];
   Y = [Y,y];
   
   if y < min_y,
      min_y = y;
      min_x = x;
   end
   
   x(1) = x(1) + step(1);

   if Grid > 1,
   
     for i=1:num_variables,
      if x(i)-step(i)/2 > Max(i),
         x(i)=Min(i);
         if i+1 <= length(x),
           x(i+1)=x(i+1)+step(i+1);
         else
            stop=1;
         end
      end
     end
   else 
     stop = 1;
   end
end

if nloops > 1,
%   Min = min_x - step/2;
%   Max = min_x + step/2;

   tmp=2*step./Grid;
   Min = min_x - step + tmp;
   Max = min_x + step - tmp;    
   
   [x, y, tmpX, tmpY]=gridsearch( Min, Max, Grid, nloops-1, fun, varargin{:});
      
   if y < min_y,
      min_y = y;
      min_x = x;
   end
   
   X = [X,tmpX];
   Y = [Y,tmpY];

end

return;
% EOF