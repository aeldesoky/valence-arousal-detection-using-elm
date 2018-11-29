function [y,dfce]=linclass( X, model)
% LINCLASS Linear classifier.
%
% Synopsis:
%  [y,dfce] = linclass( X, model)
%
% Description:
%  This function classifies input data X using linear
%  discriminant function:
%
%  y(i) = argmax W(:,y)'*X(:,i) + b(y)
%           y
%
%  where parameters W [dim x nfun] and b [1 x nfun] are given 
%  in model and nfun is number of discriminant functions.
%
%  In the binary case (nfun=1) the classification rule is following
%    y(i) = 1 if W'*X(:,i) + b >= 0
%           2 if W'*X(:,i) + b < 0
%  
%  where W [dim x 1], b [1x1] are parameters given in model.
%
% Input:
%  X [dim x num_data] Data to be classified.
%
%  model [struct] Parameters of linear classifier:
%   .W [dim x nfun] Linear term.
%   .b [nfun x 1] Bias.
%
% Output:
%  y [1 x num_data] Predicted labels.
%  dfce [nfun x num_data] Values of discriminat function.
%
% Examples:
%  trn = load('riply_trn');
%  tst = load('riply_tst');
%  model = fld( trn );
%  ypred = linclass( tst.X, model );
%  cerror( ypred, tst.y )
%  figure; ppatterns( trn ); pline( model );
%
% See also 
%  PERCEPTRON, MPERCEPTRON, FLD, ANDERSON.
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

nfun = size(model.W,2);

if nfun == 1,
  % binary case
  dfce = model.W'*X + model.b;
  y = ones(1,num_data);
  y(find(dfce < 0)) = 2;
else
  % multi-class case 
  dfce = model.W'*X + model.b(:)*ones(1,num_data);
  [dummy,y] = max( dfce );
end

return;
