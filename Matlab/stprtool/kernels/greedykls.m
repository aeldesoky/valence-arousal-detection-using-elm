function [model,Z]=greedykpca(X,y,options)
% GREEDYKLS Greedy Regularized Kernel Least Squares.
%
% Synopsis:
%  model = greedykls(X)
%  model = greedykls(X,options)
%
% Description:
%  This function approximates input vectors X in the feature
%  space using GREEDYKPCA. Then the regularized least squares
%  are applied on the approximated data. 
%
%  See help of KLS for more info about regularize least squares.
%  See help of GREEDYKPCA for more info on approximation of data
%  in the feature space.
%  
% Input:
%  X [dim x num_data] Input column vectors.
%  y [num_data x 1] Output values.
%  
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier. See HELP KERNEL for more info.
%   .arg [1 x narg] Kernel argument.
%   .m [1x1] Maximal number of base vectors (Default m=0.25*num_data).
%   .p [1x1] Depth of search for the best basis vector (Default p=m).
%   .mserr [1x1] Desired mean squared reconstruction errors of approximation.
%   .maxerr [1x1] Desired maximal reconstruction error of approximation.
%     See 'help greedyappx' for more info about the stopping conditions.
%   .verb [1x1] If 1 then some info is displayed (default 0).
% 
% Output:
%  model [struct] Kernel projection:
%   .Alpha [nsv x new_dim] Multipliers defining kernel projection.
%   .sv.X [dim x num_data] Selected subset of the training vectors.
%   .nsv [1x1] Number of basis vectors.
%   .kercnt [1x1] Number of kernel evaluations.
%   .MaxErr [1 x nsv] Maximal reconstruction error for corresponding
%     number of base vectors.
%   .MsErr [1 x nsv] Mean square reconstruction error for corresponding
%     number of base vectors.
% 
% Example:
%  x = [0:0.05:2*pi]; y = sin(x) + 0.1*randn(size(x));
%  model = greedykls(x,y(:),struct('ker','rbf','arg',1,'lambda',0.001));
%  y_est = kernelproj(x,model);
%  figure; hold on;
%  plot(x,y,'+k'); plot(x,y_est,'b'); 
%  plot(x,sin(x),'r'); plot(x(model.sv.inx),y(model.sv.inx),'ob');
%
% See also 
%   KERNELPROJ, KPCA, GREEDYKPCA.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2005, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 01-mar-2005, VF
% 22-feb-2005, VF

start_time = cputime;
[dim,num_data]=size(X);

% process input arguments
%------------------------------------
if nargin < 2, options = []; else options=c2s(options); end
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'m'), options.m = fix(0.25*num_data); end
if ~isfield(options,'p'), options.p = options.m; end
if ~isfield(options,'maxerr'), options.maxerr = 1e-6; end
if ~isfield(options,'mserr'), options.mserr = 1e-6; end
if ~isfield(options,'verb'), options.verb = 0; end
if ~isfield(options,'lambda'), options.lambda = 0.001; end

% greedy algorithm to select subset of training data
%-------------------------------------------------------

[inx,Alpha,Z,kercnt,MsErr,MaxErr] = ...
  greedyappx(X,options.ker,options.arg,...
            options.m,options.p,options.mserr,options.maxerr,options.verb); 
  
% apply ordinary linear least squares
%------------------------------
w = inv( Z*Z' + options.lambda*num_data*eye(size(Z,1))) * Z*y;

% fill up the output model
%-------------------------------------
model.Alpha = Alpha'*w;
model.nsv = length(Alpha);  
model.b = 0;
model.sv.X= X(:,inx);
model.sv.inx = inx;
model.kercnt = kercnt;
model.GreedyMaxErr = MaxErr;
model.GreedyMsErr = MsErr;
model.options = options;
model.cputime = cputime - start_time;
model.fun = 'kernelproj';

return;
% EOF
