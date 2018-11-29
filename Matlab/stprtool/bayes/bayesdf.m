function quad_model=bayesdf(model)
% BAYESDF Computes decision boundary of Bayesian classifier.
%
% Synopsis:
%  quad_model = bayesdf(model)
%
% Description:
%  This function computes parameters of decision boundary
%  of the Bayesian classifier with the following assumptions:
%   - 1/0 loss function (risk = expectation of misclassification).
%   - Binary classification.
%   - Class conditional probabilities are multivariate Gaussians.
%
%  In this case the Bayesian classifier has the quadratic 
%  discriminant function
%              f(x) = x'*A*x + B'*x + C,
%  
%  where the classification strategy is
%  q(x) = 1  if f(x) >= 0,
%       = 2  if f(x) < 0.
%
% Input:
%  model [struct] Two multi-variate Gaussians:
%   .Mean [dim x 2] Mean values.
%   .Cov [dim x dim x 2] Covariances.
%   .Prior [1x2] A priory probabilities.
%
% Output: 
%  quad_model.A [dim x dim] Quadratic term.
%  quad_model.B [dim x 1] Linear term.
%  quad_model.C [1x1] Bias.
%
% Example:
%  trn = load('riply_trn');
%  tst = load('riply_trn');
%  gauss_model = mlcgmm(trn);
%  quad_model = bayesdf(gauss_model);
%  ypred = quadclass(tst.X,quad_model);
%  cerror(ypred,tst.y)
%  figure; ppatterns(trn); pboundary(quad_model); 
%
% See also 
%  BAYESCLS, QUADCLASS
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 18-oct-2005, VF, dealing with Cov given as vector repared
% 01-may-2004, VF
% 19-sep-2003, VF
% 24. 6.00 V. Hlavac, comments into English.

% allow input to be a cell
model = c2s(model);

% univariate variances can be given as a vector
if size(model.Cov,1) == 1 && length(size(model.Cov)) < 3, 
  model.Cov = reshape(model.Cov,1,1,length(model.Cov)); 
end

M1=model.Mean(:,1);
M2=model.Mean(:,2);
C1=model.Cov(:,:,1);
C2=model.Cov(:,:,2);
P1=model.Prior(1);
P2=model.Prior(2);

quad_model.A=(1/2)*(inv(C2)-inv(C1));
quad_model.B=(M1'*inv(C1)-M2'*inv(C2))';

% Treatment of the case when apriori probabilities are zero.
% log(0)=-inf;
if P1==0,
   quad_model.C=-inf;
elseif P2==0,
   quad_model.C=inf;
else
   quad_model.C=(1/2)*(M2'*inv(C2)*M2-M1'*inv(C1)*M1)+...
      log(sqrt(det(C2)))-log(sqrt(det(C1)))+log(P1)-log(P2);
end

quad_model.fun = 'quadclass';

return;
