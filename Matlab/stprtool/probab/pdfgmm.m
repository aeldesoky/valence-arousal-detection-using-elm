function y = pdfgmm(X, model )
% PDFGMM Evaluates gaussian mixture model.
%
% Synopsis:
%  y = pdfgmm(X, model )
%
% Description:
%  This function evaluates a probability density function 
%  determined by Gaussian mixture model (GMM) for given input column 
%  vectors in X. The GMM is defined as
%         
%  y(i) = sum model.Prior(j)*pdfgauss(X(:,i),model.Mean(:,j),model.Cov(:,:,j))
%       j=1:ncomp
%
%  for all i=1:size(X,2).
% 
% Input:
%  X [dim x num_data] Input matrix of column vectors.
%  model.Mean [dim x ncomp] Means of Gaussians.
%  model.Cov [dim x dim x ncomp] Covarince matrices.
%  model.Prior [ncomp x 1] Weights of components.
%
% Output:
%  y [1 x num_data] Values of probability density function.
%
% Example:
%
% Univariate case
%  x = linspace(-5,5,100);
%  distrib = struct('Mean',[-2 3],'Cov',[1 0.5],'Prior',[0.4 0.6]);
%  y = pdfgmm(x,distrib);
%  figure; plot(x,y);
%
% Multivariate case
%  model.Mean(:,1) = [-1;-1]; model.Cov(:,:,1) = [1,0.5;0.5,1]; 
%  model.Mean(:,2) = [1;1]; model.Cov(:,:,2) = [1,-0.5;-0.5,1]; 
%  model.Prior = [0.4 0.6];
%  [Ax,Ay] = meshgrid(linspace(-5,5,100), linspace(-5,5,100));
%  y = pdfgmm([Ax(:)';Ay(:)'],model);
%  figure; surf( Ax, Ay, reshape(y,100,100)); shading interp;
%
% See also 
%  GMMSAMP, PDFGAUSS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 28-apr-2004, VF

% allows inputs to be given in cell array
model = c2s(model);

y = model.Prior(:)'*pdfgauss(X,model);

return;
