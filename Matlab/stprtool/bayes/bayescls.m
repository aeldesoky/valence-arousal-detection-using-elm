function [y, dfce] = bayescls( X, model )
% BAYESCLS Bayesian classifier with reject option.
% 
% Synopsis:
%  [y, dfce] = bayescls(X,model)
%
% Description:
%  This function implements the classifier minimizing the Bayesian risk 
%  with 0/1-loss function. It corresponds to the minimization of 
%  probability of misclassification. The input vectors X are classified 
%  into classes with the highest a posterior probabilities computed from 
%  given model.
% 
%  The model contains parameters of conditional class probabilities
%  in model.Pclass [cell 1 x num_classes] and a priory probabilities
%  in model.Prior [1 x num_classes]. 
%
%  The function
%    p = feval(model.Pclass{i}.fun, X, model.pclass{i})
%  is called to evaluate the i-the class conditional probability of X.
%  
%  It returns class labels y [1 x num_data] for each input vector
%  and matrix dfce [num_class x num_data] of unnormalized a posterior
%  probabilities
%    dfce(y,i) = Conditional_probability(X(:,i)|y)*Prior(y).
%
%  If the field model.eps exists then the Bayesian classifier 
%  with the reject option is used. The eps is penalty for the 
%  decision "don't know" which is indicated by label y = 0.
%   
% Input:
%  X [dim x num_data] Vectors to be classified.
%
%  model [struct] Describes probabilistic model:
%   .Pclass [cell 1 x num_classes] Class conditional probabilities.
%   .Prior [1 x num_classes] A priory probabilities.
%
%   .eps [1x1] (optional) Penalty of decision "don't know". 
%
% Output:
%  y [1 x num_data] Labels (1 to num_classes); 0 for "don't know".
%  dfce [num_classes x num_data] Unnormalized a posterior 
%   probabilities (see above).
%
% Example:
%  trn = load('riply_trn');
%  tst = load('riply_tst');
%  inx1 = find(trn.y==1);
%  inx2 = find(trn.y==2);
%  model.Pclass{1} = mlcgmm(trn.X(:,inx1));
%  model.Pclass{2} = mlcgmm(trn.X(:,inx2));
%  model.Prior = [length(inx1) length(inx2)]/(length(inx1)+length(inx2));
%  ypred = bayescls(tst.X,model);
%  cerror(ypred,tst.y)
% 
% See also 
%  BAYESDF, BAYESERR.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 09-jun-2004, VF
% 01-may-2004, VF
% 11-mar-2004, VF, "don't" know decision added.
% 19-sep-2003, VF

[dim,num_data]=size(X);
num_classes = length( model.Pclass );

dfce=zeros(num_classes,num_data);
% compute unnormalized a posterior probabilities
for i=1:num_classes,
  dfce(i,:) = model.Prior(i)*feval(model.Pclass{i}.fun,X,model.Pclass{i});
end

% take maximum
[tmp,y] = max(dfce);

% reject options
if isfield(model, 'eps'),
  perror = 1-tmp./sum(dfce,1);
  
  inx = find( perror > model.eps);
  y(inx) = 0;
end

return;
% EOF
