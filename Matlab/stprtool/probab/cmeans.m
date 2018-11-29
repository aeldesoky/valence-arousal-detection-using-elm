function [model,y] = cmeans(X,num_centers,Init_centers)
% CMEANS K-means clustering algorithm.
% 
% Synopsis:
%  [model,y] = cmeans(X,num_centers)
%  [model,y] = cmeans(X,num_centers,Init_centers)
%
% Description:
%  [model,y] = cmeans(X,num_centers) runs C-means clustering 
%   where inital centers are randomly selected from the 
%   input vectors X. The output are found centers stored in 
%   structure model.
%   
%  [model,y] = cmeans(X,num_centers,Init_centers) uses
%   init_centers as the starting point.
%
% Input:
%  X [dim x num_data] Input vectors.
%  num_centers [1x1] Number of centers.
%  Init_centers [1x1] Starting point of the algorithm.
%    
% Output:
%  model [struct] Found clustering:
%   .X [dim x num_centers] Found centers.
%
%   .y [1 x num_centers] Implicitly added labels 1..num_centers.
%   .t [1x1] Number of iterations.
%   .MsErr [1xt] Mean-Square error at each iteration.
%
%  y [1 x num_data] Labels assigned to data according to 
%   the nearest center.
%
% Example:
%  data = load('riply_trn');
%  [model,data.y] = cmeans( data.X, 4 );
%  figure; ppatterns(data); 
%  ppatterns(model,12); pboundary( model );
%
% See also 
%  EMGMM, KNNCLASS.
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 18-dec-2008, VF, fix: Init_centers argument used as the initial solution
% 17-jun-2007, VF, renamed from kmeans to cmeans to avoid conflicts with stats toolbox
% 12-may-2004, VF

[dim,num_data] = size(X);

% random inicialization of class centers
%-----------------------------------------------
if nargin < 3,
  inx=randperm(num_data);
  model.X = X(:,inx(1:num_centers));
  model.y = 1:num_centers;
  model.K = 1;
else
  if num_centers ~= size(Init_centers,2)
      error('Argument num_centers must be equal to the number of columns of Init_centers.'); 
  end
  model.X = Init_centers;
  model.y = 1:num_centers'
end

model.fun = 'knnclass';

old_y = zeros(1,num_data);
t = 0;

% main loop
%-------------------------
while 1,
  
  t = t+1;
  
  % classificitation
  y = knnclass( X, model );
  
  % computation of class centers
  err = 0;
  for i=1:num_centers,
    inx = find(y == i);

    if ~isempty(inx),
      
      % compute approximation error
      err = err + sum(sum((X(:,inx) - model.X(:,i)*ones(1,length(inx)) ).^2));
      
      % compute new centers
      model.X(:,i) = sum(X(:,inx),2)/length(inx);
    end
  end

  % Number of iterations and Mean-Square Error 
  model.t = t;
  model.MsErr(t) = err/num_data;
  
  if sum( abs(y - old_y) ) == 0,
    return;
  end

  old_y = y;
end

return;
% EOF
