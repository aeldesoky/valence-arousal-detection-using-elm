function data = gmmsamp(model,num_data)
% GMMSAMP Generates sample from Gaussian mixture model.
% 
% Synopsis:
%  data = gmmsamp(model,num_data)
%
% Description:
%  This function generates num_data samples from a Gaussian 
%  mixture given by structure model. It returnes samples X 
%  and a vector y of Gaussian component responsible for 
%  generating corresponding sample.
%
% Input:
%  model
%   .Mean [dim x ncomp] Mean vectors.
%   .Cov [dim x dim x ncomp] Covariance matrices. In the case of 
%     univariate mixture (dim=0) the variances can enter 
%     as a vector Cov=[var1 var2 ... var_ncomp].
%   .Prior [ncomp x 1] Weighting coefficients of Gaussians.
%  num_data [int] Number of samples.
%
% Output:
%  data.X [dim x num_data] Generated sample data.
%  data.y [1 x num_data] Identifier of Gaussian which generated 
%   given vector.
%
% Example:
%  model = struct('Mean',[-2 3],'Cov',[1 0.5],'Prior',[0.4 0.6]);
%  figure; hold on; 
%  plot([-4:0.1:5], pdfgmm([-4:0.1:5],model),'r');
%  sample = gmmsamp(model,500);
%  [Y,X] = hist(sample.X,10);
%  bar(X,Y/500);
%
% See also 
%  PDFGMM, GSAMP.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 28-apr-2004, VF

% allow input argument to be cell array
model = c2s( model);

% get dimensions
[dim,ncomp]=size(model.Mean);

% univariate variances can be given as a vector
if size(model.Cov,1) ~= size(model.Cov,2), 
  model.Cov = reshape(model.Cov,1,1,ncomp);
end

% generates randomly Gaussian componenents responsible
% for generation of samples
cump = repmat(cumsum( model.Prior(:) ),1,num_data);
rnd = ones(ncomp,1)*rand(1,num_data);

data.y = ncomp-sum(cump > rnd, 1)+1;
data.X = zeros(dim,num_data);

% generate samples
for i=1:ncomp,
  inx = find(data.y==i);
  data.X(:,inx) = gsamp(model.Mean(:,i),model.Cov(:,:,i),length(inx));
end

return;
