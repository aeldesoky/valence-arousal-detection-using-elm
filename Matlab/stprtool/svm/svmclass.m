function [y,dfce] = svmclass(X,model)
% SVMCLASS Support Vector Machines Classifier.
%
% Synopsis:
%  [y,dfce] = svmclass( X, model )
%
% Description:
%  [y,dfce] = svmclass( X, model ) classifies input vectors X
%    into classes using the multi-class SVM classifier
%      y(i) = argmax f_j(X(:,i))
%            j=1..nfun
%    where f_j are linear functions in the feature space given 
%    by the prescribed kernel function (options.ker, options.arg). 
%    The discriminant functions f_j are determined by 
%      .Alpha [nsv x nfun] ... multipliers associated to SV
%      .b [nclass] ... biases of discriminant functions.
%      .sv.X [dim x nsv] ... support vectors.
% 
%    See 'help kernelproj' for more info about valuation of the 
%    discriminant functions f_j.
%
%    In the binary case nfun=1 the binary SVM classifier is used
%      y(i) = 1 if f(X(:,i) >= 0
%           = 2 if f(X(:,i) < 0
%    where f is the disrimiant function given by Alpha [nsv x 1],
%    b [1x1] and support vectors sv.X.
%      
% Input:
%  X [dim x num_data] Input vectors to be classified.
%
%  model [struct] SVM classifier:
%   .Alpha [nsv x nfun] Multipliers associated to suport vectors.
%   .b [nfun x 1] Biases.
%   .sv.X [dim x nsv] Support vectors.
%   .options.ker [string] Kernel identifier.
%   .options.arg [1 x nargs] Kernel argument(s).
%
% Output:
%  y [1 x num_data] Predicted labels.
%  dfce [nfun x num_data] Values of discriminant functions.
%
% Example:
%  trn = load('riply_trn');
%  model = smo(trn,struct('ker','rbf','arg',1,'C',10));
%  tst = load('riply_tst');
%  ypred = svmclass( tst.X, model );
%  cerror( ypred, tst.y )
% 
% See also 
%  SMO, SVMLIGHT, SVMQUADPROG, KFD, KFDQP, MVSVMCLASS.  
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 14-may-2004, VF
% 09-May-2003, VF
% 14-Jan-2003, VF

% allows model to be given in cell
model=c2s(model);

dfce = kernelproj(X, model);
nfun = size(dfce,1);

if nfun == 1,
  % Binary case
  %-------------------------------

  y = ones(size(dfce));
  y( find( dfce < 0 )) = 2;

else  
  % Multi-class case
  %-------------------------------

  [dummy,y] = max( dfce );
end

return;
% EOF

