function model = adaboost(data,options)
% ADABOOST AdaBoost algorithm.
%
% Synopsis:
%  model = adaboost(data,options)
%
% Description:
%  This function implements the AdaBoost algorithm which
%  produces a classifier composed from a set of weak rules.
%  The weak rules are learned by a weak learner which is
%  specified in options.learner. The task of the weak learner
%  is to produce a rule with weighted error less then 0.5.
%  The Adaboost algorithm calls in each stage the weak learner
%
%     rule{t} = feval(options.learner,weight_data)
%   
%  where the structure weight_data contains
%    .X [dim x num_data] Training vectors.
%    .y [1 x num_data] Labels of training vectos (1 or 2).
%    .D [1 x num_data] Distribution (weights) over training 
%      data which defines the weighted error.
%  
%  The item rule{t}.fun must contain name of function
%  which classifies vector X by 
%   y = feval( rule{t}.fun, X, rule{t}).
%
%  It is assumed that the weak rule responds with labels 
%  1 or 2 (not 1,-1 as used in AdaBoost literature).
%
% Input:
%  data [struct] Input training data:
%   .X [dim x num_data] Training vectors.
%   .y [1 x num_data] Labels of training vectos (1 or 2).
%
%  options [struct] Parameters of the AdaBoost:
%   .learner [string] Name of the weak learner.
%   .max_rules [1x1] Maximal number of weak rules (defaul 100).
%    This paramater defines a stopping condition.
%   .err_bound [1x1] AdaBoost stops if the upper bound on the 
%    empirical error drops below the err_bound (default 0.001).
%   .learner_options Additinal options used when the weak learner
%    is called.
%   .verb [1x1] If 1 then some info is displayed.
%
% Output:
%  model [struct] AdaBoost classifier:
%   .rule [cell 1 x T] Weak classification rules.
%   .Alpha [1 x T] Weights of the rules.
%   .WeightedErr [1 x T] Weighted errors of the weak rules.
%   .Z [1 x T] Normalization constants of the distribution D.
%   .ErrBound [1 x T] Upper bounds on the empirical error.
%
% Example:
%  data = load('riply_trn');
%  options.learner = 'weaklearner';
%  options.max_rules = 100;
%  options.verb = 1;
%  model = adaboost(data,options);
%  figure; ppatterns(data); pboundary(model);
%  figure; hold on; plot(model.ErrBound,'r'); 
%  plot(model.WeightedErr);
%
% See also: 
%  ADACLASS, WEAKLEARNER.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2004, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 11-aug-2004, VF

if nargin < 2, options = []; else options = c2s(options); end
if ~isfield(options,'max_rules'), options.max_rules = 100; end
if ~isfield(options,'err_bound'), options.err_bound = 0.001; end
if ~isfield(options,'learner'), options.learner = 'weaklearner'; end
if ~isfield(options,'learner_options'), options.learner_options = []; end
if ~isfield(options,'verb'), options.verb = 0; end

% take data dimensions
[dim,num_data] = size(data.X);

% initial distribution over training samples
data.D = ones(num_data,1)/num_data;

model.Alpha =[];
model.Z = [];
model.WeightedErr = [];
model.ErrBound = [];

t = 0;
go = 1;
while go,
  t = t + 1;

  if options.verb, fprintf('rule %d: ', t); end

  % call weak learner
  if ~isempty(options.learner_options),
    rule = feval(options.learner,data,options.learner_options);
  else
    rule = feval(options.learner,data);
  end    
  
  y = feval(rule.fun,data.X,rule);

  werr = (y(:)~=data.y(:))'*data.D(:);
  if options.verb, fprintf('werr=%f', werr); end
  
  if werr < 0.5,
  
    alpha = 0.5*log((1-werr)/werr);

    % yh(i) = +1 for data.y(i) == y(i)
    % yh(i) = -1 for data.y(i) ~= y(i)
    yh = 2*(y(:) == data.y(:))-1;
    weights = data.D.*exp(-alpha*yh(:));

    % normalization constant
    Z = sum(weights);
    data.D = weights/Z;

    % upper bound on the training error
    err_bound = prod(model.Z);
    
    % store variables
    model.Z = [model.Z; Z];
    model.Alpha = [model.Alpha;alpha];
    model.rule{t} = rule;
    model.ErrBound = [model.ErrBound; err_bound];
    
    % stopping conditions
    if t >= options.max_rules,
      go = 0;
      model.exitflag = 1;
    elseif err_bound <= options.err_bound,
      go = 0;
      model.exitflag = 2;
    end
    
    if options.verb,
     fprintf(', alpha=%f, err_bound=%f\n',alpha,err_bound);
    end
    
  else
    % the weighted error is over 0.5
    if options.verb, fprintf('>= 0.5 thus stop.\n'); end
    
    go = 0;
    model.exitflag = 0;
  end

  model.WeightedErr = [model.WeightedErr; werr];

end

model.fun = 'adaclass';

return;
% EOF