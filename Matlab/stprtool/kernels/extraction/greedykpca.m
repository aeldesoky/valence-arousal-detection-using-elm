function [model,Z]=greedykpca(X,options)
% GREEDYKPCA Greedy Kernel Principal Component Analysis.
%
% Synopsis:
%  model = greedykpca(X)
%  model = greedykpca(X,options)
%
% Description:
%  This function implements a greedy kernel PCA algorithm. 
%  The input data X are first approximated by GREEDYKPCA in the 
%  feature space and second the ordinary PCA is applyed on the 
%  approximated data. This algorithm has the same objective function 
%  as the ordinary Kernel PCA but, in addition, the number of data in 
%  the resulting kernel expansion is limited. 
%
%  For more info refer to V.Franc: Optimization Algorithms for Kernel 
%  Methods. Research report. CTU-CMP-2005-22. CTU FEL Prague. 2005.
%  ftp://cmp.felk.cvut.cz/pub/cmp/articles/franc/Franc-PhD.pdf .
%  
% Input:
%  X [dim x num_data] Input column vectors.
%  
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier. See 'help kernel' for more info.
%   .arg [1 x narg] Kernel argument.
%   .m [1x1] Maximal number of base vectors (Default m=0.25*num_data).
%   .p [1x1] Depth of search for the best basis vector (p=m).
%   .mserr [1x1] Desired mean squared reconstruction errors of approximation.
%   .maxerr [1x1] Desired maximal reconstruction error of approximation.
%     See 'help greedyappx' for more info about the stopping conditions.
%   .verb [1x1] If 1 then some info is displayed (default 0).
% 
% Output:
%  model [struct] Kernel projection:
%   .Alpha [nsv x new_dim] Multipliers defining kernel projection.
%   .b [new_dim x 1] Bias the kernel projection.
%   .sv.X [dim x num_data] Seleted subset of the training vectors..
%   .nsv [1x1] Number of basis vectors.
%   .kercnt [1x1] Number of kernel evaluations.
%   .MaxErr [1 x nsv] Maximal reconstruction error for corresponding
%     number of base vectors.
%   .MsErr [1 x nsv] Mean square reconstruction error for corresponding
%     number of base vectors.
% 
% Example:
%  X = gencircledata([1;1],5,250,1);
%  model = greedykpca(X,struct('ker','rbf','arg',4,'new_dim',2));
%  X_rec = kpcarec(X,model);             
%  figure; 
%  ppatterns(X); ppatterns(X_rec,'+r');
%  ppatterns(model.sv.X,'ob',12);
%
% See also 
%   KERNELPROJ, KPCA, GREEDYAPPX.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 09-sep-2005, VF
% 19-feb-2005, VF
% 10-jun-2004, VF
% 05-may-2004, VF
% 14-mar-2004, VF

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

% greedy algorithm to select subset of training data
%-------------------------------------------------------

[inx,Alpha,Z,kercnt,MsErr,MaxErr] = ...
  greedyappx(X,options.ker,options.arg,...
            options.m,options.p,options.mserr,options.maxerr,options.verb); 
  
% apply ordinary PCA
%------------------------------
mu = sum(Z,2)/num_data;
Z=Z-mu*ones(1,num_data);

S = Z*Z';
[U,D,V]=svd(S);

model.eigval=diag(D);
sum_eig = triu(ones(size(Z,1),size(Z,1)),1)*model.eigval;
model.MsErr = MsErr(end)+sum_eig/num_data;

options.new_dim = min([options.new_dim,size(Z,1)]);

V = V(:,1:options.new_dim);

% fill up the output model
%-------------------------------------
model.Alpha = Alpha'*V;
model.nsv = length(inx);  
model.b = -V'*mu;
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
