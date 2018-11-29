function X=gsamp(varargin)
% GSAMP Generates sample from Gaussian distribution.
% 
% Synopsis:
%  X = gsamp( Mean, Cov, num_data )
%  X = gsamp( model, num_data )
%
% Description:
%  X = gsamp(Mean,Cov,num_data) generates num_data samples from 
%  a multi-variate Gassian distribution given by mean vector 
%  Mean [dim x 1] and covariance matrix Cov [dim x dim]. 
%
%  X = gsamp(model,num_data) assumes that parameters of Gaussian
%  are given in structure with fields model.Mean a model.Cov.
%  
% Example:
%  model = struct('Mean',1,'Cov',2);
%  figure; hold on;
%  plot([-4:0.1:5],pdfgauss([-4:0.1:5],model),'r');
%  [Y,X] = hist(gsamp(model,500),10);
%  bar(X,Y/500);
%
% See also 
%  PDFGAUSS, GMMSAMP.

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 28-apr-2004, VF, adopted from P.Krizek 

if nargin > 2,
  Mean = varargin{1};
  Cov = varargin{2};
  num_data = varargin{3};
else
  Mean = varargin{1}.Mean;
  Cov = varargin{1}.Cov;
  num_data = varargin{2};
end

% get dimension
dim = length(Mean);

% compute eigen values and vectors
[U,L] = eig(Cov);

% dewhitening transform
X = inv(U')*sqrt(L)*randn(dim,num_data)+repmat(Mean,1,num_data);


return;

