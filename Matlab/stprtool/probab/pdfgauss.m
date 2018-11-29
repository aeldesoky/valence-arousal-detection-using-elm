function y = pdfgauss(X, arg1, arg2 )
% PDFGAUSS Evaluates multivariate Gaussian distribution.
%
% Synopsis:
%  y = pdfgauss(X, Mean, Cov)
%  y = pdfgauss(X, model )
%
% Description:
%  y = pdfgauss(X, Mean, Cov) evaluates a multi-variate Gaussian 
%  probability density function(s) for given input column vectors in X.
%  Mean [dim x ncomp] and Cov [dim x dim x ncomp] describe a set of 
%  ncomp Gaussian distributions to be evaluted such that
%
%  y(i,j) = exp(-0.5(mahalan(X(:,j),Mean(:,i),Cov(:,:,i) )))/norm_const
%
%  where i=1:ncomp and j=1:size(X,2). If the Gaussians are
%  uni-variate then the covariaves can be given as a vector
%  Cov = [Cov_1, Cov_2, ..., Cov_comp].
%
%  y = pdfgauss( X, model ) takes Gaussian parameters from structure
%  fields model.Mean and model.Cov.
%
% Input:
%  X [dim x num_data] Input matrix of column vectors.
%  Mean [dim x ncomp] Means of Gaussians.
%  Cov [dim x dim x ncomp] Covarince matrices.
%
% Output:
%  y [ncomp x num_data] Values of probability density function.
%
% Example:
% 
% Univariate case
%  x = linspace(-5,5,100);
%  y = pdfgauss(x,0,1);
%  figure; plot(x,y)
%
% Multivariate case
%  [Ax,Ay] = meshgrid(linspace(-5,5,100), linspace(-5,5,100));
%  y = pdfgauss([Ax(:)';Ay(:)'],[0;0],[1 0.5; 0.5 1]);
%  figure; surf( Ax, Ay, reshape(y,100,100)); shading interp;
%
% See also 
%  GSAMP, PDFGMM.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 28-apr-2004, VF

% process input arguments
if nargin < 3,
  arg1 = c2s(arg1);
  Mean = arg1.Mean;
  Cov =  arg1.Cov;
else
  Mean = arg1;
  Cov =  arg2;
end

% get dimensions
[dim,num_data] = size(X);
ncomp = size(Mean,2);

% univariate variances can be given as a vector
if size(Cov,1) ~= size(Cov,2), Cov = reshape(Cov,1,1,ncomp); end

% alloc memory
y = zeros(ncomp,num_data);

% evaluate pdf for each component
for i=1:ncomp,
  dist = mahalan(X,Mean(:,i),Cov(:,:,i));
  y(i,:) = exp(-0.5*dist)/sqrt((2*pi)^dim*det(Cov(:,:,i)));
end

return;
