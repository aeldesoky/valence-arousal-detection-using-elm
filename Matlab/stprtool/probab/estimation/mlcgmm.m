function model=mlcgmm(data,cov_type)
% MLCGMM Maximal Likelihood estimation of Gaussian mixture model.
% 
% Synopsis:
%  model = mlcgmm(X)
%  model = mlcgmm(X,cov_type)
%  model = mlcgmm(data)
%  model = mlcgmm(data,cov_type)
% 
% Description:
%  It computes Maximum Likelihood estimation of parameters
%  of Gaussian mixture model for given labeled data sample
%  (complete data).
%
%  model = mlcgmm(X) computes parameters (model.Mean,model.Cov)
%   of a single Gaussian distribution for given sample of column 
%   vectors X (all labels are assumed to be 1).
%
%  model = mlcgmm(X,cov_type) specifies shape of covariance matrix:
%   cov_type = 'full'      full covariance matrix (default)
%   cov_type = 'diag'      diagonal covarinace matrix
%   cov_type = 'spherical' spherical covariance matrix
%
%  model = mlcgmm(data) computes parameters of a Gaussian mixture model
%   from a given labeled data sample
%     data.X ... samples,
%     data.y .. labels.
%   It estimates parameters of ncomp=max(data.y) Gaussians and
%   a priory probabilities Prior [1 x ncomp] using Maximum-Likelihood 
%   principle.
%
% Input:
%  X [dim x num_data] Data sample.
%  data.X [dim x num_data] Data sample.
%  data.y [1 x num_data] Data labels.
%  cov_type [string] Type of covariacne matrix (see above).
%
% Output:
%  model [struct] Estimated Gaussian mixture model:
%   .Mean [dim x ncomp] Mean vectors.
%   .Cov [dim x dim x ncomp] Covariance matrices.
%   .Prior [1 x ncomp] Estimated a priory probabilities.
%  
% Example:
%  data = load('riply_trn');
%  model = mlcgmm( data );
%  figure; hold on; ppatterns(data); pgauss( model );
%  figure; hold on; ppatterns(data); pgmm( model );
%
% See also 
%  EMGMM, MMGAUSS, PDFGMM.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 17-aug-2004, VF, labels y do not have to form a sequence 1,2,...,max_y
% 2-may-2004, VF
% 29-apr-2004, VF
% 19-sep-2003, VF
% 27-feb-2003, VF

% processing of  inputs
data=c2s(data);

if ~isstruct(data),
  data.X = data;
  data.y = ones(1,size(data.X,2));
end
 
if nargin < 2, cov_type = 'full'; end

[dim,num_data] = size(data.X);

labels = unique(data.y);
model.Mean = zeros(dim,length(labels));
model.Cov = zeros(dim,dim,length(labels));
for i=1:length(labels),
   
   inx = find(data.y==labels(i));
   n = length(inx);

   model.Mean(:,i) = sum(data.X(:,inx),2)/n;

   XC=data.X(:,inx)-model.Mean(:,i)*ones(1,n);

   switch cov_type,
     case 'full', 
       model.Cov(:,:,i) = XC*XC'/n;
     case 'diag', 
       model.Cov(:,:,i) = diag(sum(XC.^2,2)/n);
     case 'spherical'
       model.Cov(:,:,i) = eye(dim,dim)*sum(sum(XC.^2))/(n*dim);
     otherwise
       error('Wrong cov_type.');
   end
   
   model.Prior(i) = n/num_data;
   model.y(i) = labels(i);
end

model.cov_type = cov_type;
model.fun = 'pdfgmm';

return; 
