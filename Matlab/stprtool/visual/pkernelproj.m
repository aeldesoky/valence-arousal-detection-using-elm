function pkernelproj(model,options)
% PKERNELPROJ Plots isolines of kernel projection.
%
% Synopsis:
%  pkernelproj(model)
%  pkernelproj(model,options)
%
% Description:
%  This function plots isolines corresponding to points
%  having equal value of projection onto kernel basis
%  (e.g. onto principal component in Kernel PCA). 
%  See 'help kernelproj' for more help on kernel projection.
%
% Input:
%  model [struct] Defines kernel data projection. 
%
% options [struct] Controls apperance:
%  .contours [1x1] Number of plotted contours (Default 10).
%  .xgrid [1x1] Density of sampling in x-axis (default 25).x
%  .ygrid [1x1] Density of sampling in y-axis (default 25).
% 
% Example:
%  help kpca;
%  
% See also 
%  KERNELPROJ, KPCA, GDA.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 4-may-2004, VF
% 22-jan-2003, VF
% 8-July-2001, V.Franc 

if nargin < 2, options = []; end
if ~isfield(options,'contours'), options.contours=10; end
if ~isfield(options,'xgrid'), options.xgrid=25; end
if ~isfield(options,'ygrid'), options.ygrid=25; end

oldhold=ishold; 
hold on;

% Testing points - grid [xgrid x ygrid].
w = axis;
xrange=w(1):(w(2)-w(1))/options.xgrid:w(2);
yrange=w(3):(w(4)-w(3))/options.ygrid:w(4);

[X,Y] = meshgrid(xrange,yrange);
Xtst=[reshape(X,1,prod(size(X)));reshape(Y,1,prod(size(Y)))];
% 
Z = kernelproj(Xtst,model);

% Plot contours of selected features.
for i=1:size(Z,1),
  map = reshape(Z(i,:), length(yrange), length(xrange));
  contour(xrange, yrange, map, options.contours, marker_color(i));
end

% 
if ~oldhold, hold off; end

return;
