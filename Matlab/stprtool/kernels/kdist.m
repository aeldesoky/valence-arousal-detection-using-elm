function d=kdist(X,model)
% KDIST Computes distance between vectors in kernel space.
%
% Synopsis:
%  d = kdist(X,model)
%
% Description:
%  It computes distance between vectors mapped into the feature 
%  space induced by the kernel function (model.options.ker,
%  model.options.arg). The distance is computed between images
%  of vectors X [dim x num_data] mapped into feature space
%  and a point in the feature space given by model:
%
%   d(i) = kernel(X(:,i),X(:,i)) 
%          - 2*kernel(X(:,i),models.sv.X)*model.Alpha + b,
%
%  where b [1x1] is assumed to be equal to 
%   model.b = model.Alpha'*kernel(model.sv.X)*model.Alpha.
%
% Input:
%  X [dim x num_data] Input vectors.
%  model [struct] Deternines a point of the feature space:
%   .Alpha [nsv x 1] Multipliers.
%   .sv.X [dim x nsv] Vectors.
%   .b [1x1] Bias.
%   .options.ker [string] Kernel identifier (see 'help kernel').
%   .options.arg [1 x nargs] Kernel argument(s).
%
% Output:
%  d [num_data x 1] Distance between vectors in the feature space.
%
% Example:
%  data = load('riply_trn');
%  model.Alpha = dualmean(size(data.X,2));
%  model.sv.X = data.X;
%  model.options.ker = 'rbf';
%  model.options.arg = 0.25;
%  model.b = model.Alpha'*kernel(data.X,'rbf',0.25)*model.Alpha;
%  [Ax,Ay] = meshgrid(linspace(-5,5,100), linspace(-5,5,100));
%  dist = kdist([Ax(:)';Ay(:)'],model);
%  figure; hold on; 
%  ppatterns(data.X); contour( Ax, Ay, reshape(dist,100,100));
% 
% See also 
%  MINBALL.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 25-aug-2004, VF, MINBALL added to See also 
% 16-may-2004, VF
% 26-feb-2003, VF
% 13-sep-2002, VF
% 15-jun-2002, VF

[dim,num_data]=size(X);

x2 = diagker( X, model.options.ker, model.options.arg);

Ksvx = kernel( X, model.sv.X, model.options.ker, model.options.arg);

d = sqrt( x2 - 2*Ksvx*model.Alpha(:) + model.b*ones(num_data,1) );

return;
% EOF