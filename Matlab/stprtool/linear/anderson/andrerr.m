function [err,r,inx] = andrerr( model, distrib )
% ANDRERR Classification error of the Generalized Anderson's task.
%
% Synopsis:
%  [err,r,inx] = andrerr( model, distrib )
%
% Description:  
%  This function computes the classification error of
%  the given linear classifier and underlying set of Gaussian 
%  distributions as defined in the Generalized Anderson's 
%  task [SH10].
%
% Input:
%  model [struct] Linear classifier:
%   .W [dim x 1] Normal vector the separating hyperplane.
%   .b [real] Bias the hyperplane.
%  
%  distrib [struct] Set of Gaussians with assigned binary labels:
%   .Mean [dim x ncomp] Mean vectors.
%   .Cov [dim x dim x ncomp] Covariance matrices.
%   .y [1 x ncomp] Lables of Gaussians (1 or 2).
%  
% Output:
%  err [real] Probability of misclassification.
%  r [real] Mahalanobis distance of the cloasest Gaussian.
%  inx [int] Index of the cloasest Gaussian.
%
% Example:
%  distrib = load('mars');
%  model = eanders(distrib,{'err',0.06'});
%  figure; pandr( model, distrib );
%  error = andrerr( model, distrib )
%
% See also 
%  ANDRORIG, GANDERS, EANDERS, GGRADANDR.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 4-may-2004, VF
% 17-sep-2003, VF

if ~isfield(distrib,'y'), distrib.y = [1,2]; end
[dim,ncomp] = size(distrib.Mean);

Radius = zeros(ncomp,1);

for i=1:ncomp,
  
  if distrib.y(i) == 1,
    Radius(i) = (model.W'*distrib.Mean(:,i)+model.b)/...
        sqrt(model.W'*distrib.Cov(:,:,i)*model.W);
  else
    Radius(i) = -(model.W'*distrib.Mean(:,i)+model.b)/...
        sqrt(model.W'*distrib.Cov(:,:,i)*model.W);
  end
    
end

[r,inx] = min( Radius );
err=1-cdf('norm',r,0,1);

return;
