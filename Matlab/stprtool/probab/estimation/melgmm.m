function model=melgmm(X,Alpha,cov_type)
% MELGMM Maximizes Expectation of Log-Likelihood for Gaussian mixture.
% 
% Synopsis:
%  model = melgmm(X,Alpha)
%  model = melgmm(X,Alpha,cov_type)
% 
% Description:
%  model = melgmm(X,Alpha) maximizes expectation of log-likelihood 
%  function for Gaussian mixture model
%                        
%   (Mean,Cov,Prior) =  argmax  F(Mean,Cov,Prior)
%                    Mean,Cov,Prior 
%
%  where
%   F = sum sum Alpha(j,i)*log(pdfgauss(X(:,i),Mean(:,y),Cov(:,:,y)))
%        y   i 
%
%  The solution is returned in the structure model with fields
%  Mean [dim x ncomp], Cov [dim x dim x ncomp] and Prior [1 x ncomp].
%
%  model = melgmm(X,Alpha,cov_type) specifies covariance matrix:
%   cov_type = 'full'      full covariance matrix (default)
%   cov_type = 'diag'      diagonal covarinace matrix
%   cov_type = 'spherical' spherical covariance matrix
%
% Input:
%  X [dim x num_data] Data sample.
%  Alpha [ncomp x num_data] Distribution of hidden state given sample.
%  cov_type [string] Type of covariacne matrix (see above).
%
% Output:
%  model [struct] Gaussian mixture model:
%   .Mean [dim x ncomp] Mean vectors.
%   .Cov [dim x dim x ncomp] Covariance matrices.
%   .Prior [1 x ncomp] Distribution of hidden state.
%
% See also 
%  EMGMM, MLCGMM.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 30-apr-2004, VF
% 19-sep-2003, VF
% 27-feb-2003, VF

% Processing of input arguments 
%----------------------------------------
if nargin < 3, cov_type = 'full'; end
[dim, num_data] = size( X );
  
% ------------------------------------
ncomp = size(Alpha,1);

model.Mean = zeros(dim,ncomp);
model.Cov = zeros(dim,dim,ncomp);

for i=1:ncomp,

  nconst = sum( Alpha(i,:) );
  if ~nconst,
    model.Mean(:,i) = NaN*ones(dim,1);
    model.Cov(:,:,i) = NaN*ones(dim,dim);
    model.Prior(i) = 0;
  else
    model.Mean(:,i) = X*Alpha(i,:)'/nconst;
    XC = X - model.Mean(:,i)*ones(1,num_data);
    
    switch cov_type,
      case 'full'
        model.Cov(:,:,i) = (XC.*(repmat(Alpha(i,:),dim,1)))*XC'/nconst;
      case 'diag'
        model.Cov(:,:,i)=diag(sum(XC.*(ones(dim,1)*Alpha(i,:)).*XC,2))/nconst;
      case 'spherical'        
        model.Cov(:,:,i) = eye(dim)*...
          sum(sum(XC.^2.*(ones(dim,1)*Alpha(i,:)) ))/(nconst*dim);
      otherwise
        error('Wrong cov_type.');
    end
  
    model.Prior(i) = nconst/num_data;
  end
end  
  

return; 
