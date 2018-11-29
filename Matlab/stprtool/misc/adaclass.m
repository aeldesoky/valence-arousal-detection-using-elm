function [y,dfce] = adaclass(X,model)
% ADACLASS AdaBoost classifier.
%
% Synopsis:
%  [y,dfce] = adaclass(X,model)
%
% Description:
%  This function implements the AdaBoost classifier which
%  its discriminant function is composed of a weighted sum
%  of binary rules. It is assumed here that the binary rules
%  respond with label 1 or 2 (not 1 and -1 as used in 
%  AdaBoost literature).
%
% Input:
%  X [dim x num_data] Vectors to be classified.
%  model [struct] AdaBoost classifier:
%   .rule [cell 1 x T] Binary weak rules.
%   .Alpha [1 x T] Weights of the weak rules.
%   .fun = 'adaclass' (optinal).
%
% Output:
%  y [1 x num_data] Predicted labels.
%  dfce [1 x num_data] Values of weighted sum of 
%   weak rules; y(i) = 1 if dfce(i) >= 0 and
%   y(i) = 2 if dfce(i) < 0.
%
% Example:
%  trn_data = load('riply_trn');
%  tst_data = load('riply_tst');
%  options.learner = 'weaklearner';
%  options.max_rules = 50;
%  options.verb = 1;
%  model = adaboost(trn_data, options);
%  ypred1 = adaclass(trn_data.X,model);
%  ypred2 = adaclass(tst_data.X,model);
%  trn_err = cerror(ypred1,trn_data.y)
%  tst_err = cerror(ypred2,tst_data.y)
%
% See also: 
%  ADABOOST, WEAKLEARNER.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2004, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 25-aug-2004, VF
% 11-aug-2004, VF

dfce = [];
for i=1:length(model.rule),

  curr_y = feval(model.rule{i}.fun,X,model.rule{i});
  curr_y = 3-2*curr_y;
  
  if isempty(dfce),
    dfce = curr_y*model.Alpha(i);
  else
    dfce = dfce + curr_y*model.Alpha(i);
  end
end

y = zeros(size(dfce));
y(find(dfce>=0)) = 1;
y(find(dfce<0)) = 2;

return;
% EOF