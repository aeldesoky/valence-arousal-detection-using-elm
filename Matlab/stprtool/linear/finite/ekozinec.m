function model=ekozinec(data,options,init_model)
% EKOZINEC Kozinec's algorithm for eps-optimal separating hyperplane.
%
% Synopsis:
%  model = ekozinec(data)
%  model = ekozinec(data,options)
%  model = ekozinec(data,options,init_model)
%
% Description:
%  This function is implementation of the Kozinec's algorithm
%  with eps-optimality stopping condition [SH10]. The algorithm 
%  finds the eps-optimal separating hyperplane.
% 
%  model=ekozinec(data) the Kozinec's rule is used to find the closest 
%   points w1, w2 from the convex hulls of the vectors from the first and 
%   the second class. The found points determine the optimal separating 
%   hyperplane. 
% 
%  model=ekozinec(data,options) specifies stopping conditions of
%   the algorithm in structure options:
%    .eps [1x1] ... controls how close is the found solution to
%     the optimal hyperplane in terms of margin 
%     (default eps=0.01). The options for eps are: 
%       eps > 0 ... eps-optimal hyperplane is sought for.
%       eps == 0 ... algorithm converges to the optimal hyperplane (but it
%                    does not have to stop in finite number of iterations).
%       eps < 0 ... algorithm stops when the separating hyperplane 
%                   is found (zero training error) regardless the margin 
%                   so it solves the same task as the ordinary Perceptron.
%    .tmax [1x1]... maximal number of iterations.
%
%  model = ekozinec(data,options,init_model) specifies initial model
%   which must contain:
%    .W1 [dim x 1] ... Vector from the first convex hull.
%    .W2 [dim x 1] ... Vector from the second convex hull.
%
% Input:
%  data [struct] Labeled (binary) training data. 
%   .X [dim x num_data] Input vectors.
%   .y [1 x num_data] Labels (1 or 2).
%
%  options [struct] 
%   .eps [real] Controls how closeness to the optimal hypeprlane (see above).
%   .tmax [1x1] Maximal number of iterations (default tmax=inf).
%  
%  init_model [struct] Initial model; must contain items
%    .W1 [dim x 1], .W2 [dim x 1] see above.
%
% Output:
%  model [struct] Binary linear classifier:
%   .W [dim x 1] Normal vector of hyperplane.
%   .b [1x1] Bias of hyperplane.
%  
%   .W1 [dim x 1] The nearest vector of the first convex hull.
%   .W2 [dim x 1] The nearest vector of the second convex hull.
%   .margin [1x1] Margin of the found hyperplane.
%   .exitflag [1x1] 1 ... eps-optimality condition satisfied or separating
%                         hyperplane has been found 
%                   0 ... number of iterations exceeded tmax.
%   .t [1x1] Number of iterations.
%
% Example:
%  data = genlsdata( 2, 50, 1);
%  model = ekozinec(data, struct('eps',0.01));
%  figure; ppatterns(data); pline(model); 
%
% See also 
%  PERCEPTRON, MPERCEPTRON, LINCLASS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 19-may-2004, VF
% 3-may-2004, VF
% 17-Sep-2003, VF
% 17-Feb-2003, VF
% 16-Feb-2003, VF
% 21-apr-2001, V.Franc, created

% get data dimensions
[dim,num_data] = size(data.X);
inx1=find(data.y==1);
inx2=find(data.y==2);

% Process input arguments
% --------------------------
if nargin < 2,  options = []; else options = c2s(options); end
if ~isfield(options,'tmax'), options.tmax = inf; end
if ~isfield(options,'eps'), options.eps = 0.01; end
if ~isfield(options,'verb'), options.verb = 0; end

if nargin < 3,
  % creates init model
  % --------------------------
  inx1=find(data.y==1);   inx2=find(data.y==2);
  model.W1 = data.X(:,inx1(1));
  model.W2 = data.X(:,inx2(1));
else
  % take init model from input
  %--------------------------------
  model = init_model;
end

model.t = 0; 
model.exitflag = 0;

% main loop
%-----------------------------
while model.exitflag == 0 & options.tmax > model.t,
   
  model.t = model.t + 1;

  dW = (model.W1 - model.W2);
  norm_dW = norm( dW );
  
  projx = data.X'*dW;
  
  projx(inx1) = (projx(inx1) - model.W2'*dW)/norm_dW;
  projx(inx2) = (-projx(inx2) + model.W1'*dW)/norm_dW;
  
  [min_proj, min_inx] = min(projx);

  % bound for separating or eps-optimal separating hyperplane 
  if options.eps < 0, bound = norm_dW/2; else bound=norm_dW-options.eps/2; end
  
  if min_proj <= bound,
    xt=data.X(:,min_inx);

    % Updata - Kozinec's rule
    if data.y(min_inx) == 1,
      W1x = model.W1-xt;
      k = min(1, dW'*W1x/ (W1x'*W1x)); 
      model.W1 = model.W1 * (1 - k) + xt * k;
    else
      W2x = model.W2-xt;
      k = min(1, -dW'*W2x / (W2x'*W2x));
      model.W2 = model.W2 * (1 - k) + xt * k;
    end

    model.exitflag = 0;
  else
    model.exitflag = 1;
  end

  % print info
  if options.verb == 1 & mod(model.t,100) == 0,
     fprintf('iter %d: upper_bound = %f, lower_bound = %f, dif = %f\n', ...
       model.t, norm_dW/2, min_proj/2, (norm_dW-min_proj)/2 );  
  end
  
end

model.b = 0.5 * (model.W2'*model.W2 - model.W1'*model.W1);
model.W = model.W1 - model.W2;
model.margin = min_proj/2;
model.fun = 'linclass';

return;

