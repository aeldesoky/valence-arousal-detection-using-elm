function [y,dfce]=quadclass( X, model)
% QUADCLASS Quadratic classifier.
%
% Synopsis:
%  [y,dfce] = quadclass(X,model)
%
% Description:
%  This function classifies input data X using quadratic
%  discriminant function:
%
%  y(i) = argmax X(:,i)'*A(:,:,y)*X(:,i) + X(:,i)'*B(:,y) + C(y)
%           y
%
%  where parameters A [dim x dim x nfun], B [dim x nfun]
%  and model C [1 x nfun] are given in model and nfun is
%  number of discriminant functions.
%
%  In the binary case (nfun=1) the classification rule is following
%    y(i) = 1 if X(:,i)'*A*X(:,i) + X(:,i)'*B + C >= 0
%           2 if X(:,i)'*A*X(:,i) + X(:,i)'*B + C < 0
%  
%  where A [dim x dim], B [dim x 1] and C [1x1] are parameters
%  given in model.
%
% Input:
%  X [dim x num_data] Data to be classified.
%
%  model [struct] Describes quadratic classifier:
%   .A [dim x dim x nfun] Parameter of quadratic term.
%   .B [dim x nfun] Parameter of linear term.
%   .C [1 x nfun] Bias.
%
% Output:
%  y [1 x num_data] Predicted labels.
%  dfce [nfun x num_data] Values of discriminat function.
%
% Example:
%  trn = load('riply_trn');
%  tst = load('riply_tst');
%  gauss_model = mlcgmm(trn);
%  quad_model = bayesdf(gauss_model);
%  ypred = quadclass(tst.X, quad_model);
%  cerror(ypred, tst.y)
%  figure; ppatterns(trn); pboundary(quad_model);
%
% See also 
%  QMAP, LIN2QUAD, LINCLASS
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 2-may-2004, VF

% allow model to be gievn as a cell
model = c2s(model);

[dim, num_data] = size(X);

nfun = size(model.A,3);

if nfun == 1,
  % binary case
  dfce = sum((model.A*X).*X,1) + model.B'*X + model.C;
  y = ones(1,num_data);
  y(find(dfce< 0)) = 2;
else
  % multi-class case 
  dfce = zeros( nfun, num_data );
  
  for i=1:nfun,
    dfce(i,:) = sum((model.A(:,:,i)*X).*X,1) + model.B(:,i)'*X + model.C(i);
  end
  
  [dummy,y] = max( dfce );
  
end

return;
% EOF