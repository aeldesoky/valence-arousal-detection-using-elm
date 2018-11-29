function model = eanders(distrib, options, init_model)
% EANDERS Epsilon-solution of the Generalized Andersson's task.
% 
% Synopsis:
%  model = eanders(distrib )
%  model = eanders(distrib, options)
%  model = eanders(distrib, options, init_model)
%
% Description:
%  This function is an implementation of the Schlesinger's iterative
%  algorithm which finds the epsilon-solution of the Generalized 
%  Anderson's task using the Kozinec's algorithm [SH10]. 
%
%  The goal of the GAT is find the binary linear classification
%  rule (g(x)=sgn(W'*x+b) with minimal probability of 
%  misclassification. The conditional probabilities are known to 
%  be Gaussians their paramaters belong to a given set of parameters. 
%  The true parameters are not known. The linear rule which 
%  guarantes the minimimal classification error for the worst 
%  possible case (the worst configuration of Gaussains) is 
%  sought for. 
%
% Input:
%  distrib [struct] Input set of labeld (1 or 2) Gassians:
%   .Mean [dim x ncomp] Mean veactors.
%   .Cov  [dim x dim x ncomp] Covariance matrices.
%   .y [1 x ncomp] Labels of Gaussian (1 or 2).
% 
%  options [struct] Determine stopping conditions:
%   .tmax [1x1] Maximal number of iterations.
%   .err [1x1] Desired classification error; must be 0<err<0.5; 
%     (default 0.05).
%
%  init_model [struct] Initial model:
%    W1, W2, t.
%
% Output:
%  model [struct] Binary linear classifier:
%   .W [dim x 1] Normal vector of the linear rule (hypeplane).
%   .b [1x1] Bias of the rule (shift from the origin).
% 
%   .t [1x1] Number of used iterations.
%   .exitflag [1x1] 1 ... solution with desired err was found.
%                   0 ... maximal number of iterations exceeded.
%                   -1 ... solution does not exist.
%   .W1, .W2 Auxciliary vectors; W=W1-W2.
%
% Example:
%  distrib = load('mars');
%  model = eanders(distrib,struct('err',0.06'));
%  figure; pandr( model, distrib );
%
% See also 
%  ANDRORIG, GANDERS, GGRADANDR, LINCLASS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 21-may-2004, VF
% 16-sep-2003, VF

if nargin < 2, options = []; else options=c2s(options); end
if ~isfield(options,'err'), options.err=0.05; end
if ~isfield(options,'tmax'), options.tmax=inf; end
if ~isfield(options,'zero_th'), options.zero_th=1e-6; end

% computes Mahalanobis distance correponding to the desired 
% misclassification error
desired_r = -icdf('norm',options.err,0,1);

% get dimension and number of distributions
[dim,ncomp]=size(distrib.Mean);

t = 0;
if nargin == 3,  
  t = init_model.t;
  W1 = init_model.W1;
  W2 = init_model.W2;
end

if t==0,
  t=1;
  W1 = mean(distrib.Mean(:,find( distrib.y==1)),2);
  W2 = mean(distrib.Mean(:,find( distrib.y==2)),2);
end

exitflag=0;
while exitflag == 0 & t < options.tmax,
   t=t+1;
  
   % compute f(x)=W'*x + b
   W=W1-W2;
   b=0.5*(W2'*W2 - W1'*W1);

   exitflag = 1;
   i=0;
   while exitflag ==1 & i < ncomp,
     i=i+1;

     mu_i = distrib.Mean(:,i);
     C_i = distrib.Cov(:,:,i);

     % denominator
     den = sqrt( W'*C_i*W );

     if den > options.zero_th,
       if distrib.y(i) == 1,  % 1st class
         r_i = ( W'*mu_i + b )/den;
         if desired_r >= r_i,  % stopping condition
            x0 = mu_i - ( desired_r/den )*C_i*W;       % overlapping point
            k = min([((W1-W2)'*(W1-x0))/( (W1-x0)'*(W1-x0) ), 1]);
            W1 = W1*(1-k) + x0*k;
            exitflag = 0;
         end
       elseif distrib.y(i)==2,
         r_i= -( W'*mu_i + b)/den;
         if desired_r >= r_i,
            x0 = mu_i + ( desired_r/den )*C_i*W;
            k = min([((W2-W1)'*(W2-x0))/( (W2-x0)'*(W2-x0)), 1]);
            W2 = W2*(1-k) + x0*k;
            exitflag = 0;
         end
       end 
     else % if den ~= 0,
       % solution does not exist - overlapping classes
       exitflag = -1;
     end

   end % while exitflag == 1,

end % while

% compute f(x)=W'*x + b
model.W=W1-W2;
model.b=0.5*(W2'*W2 - W1'*W1);
model.t = t;
model.exitflag = exitflag;
model.W1 = W1;
model.W2 = W2;
model.options = options;

[model.err,model.r]=andrerr(model,distrib);
if model.err < options.err, model.exitflag = 1; end
model.fun = 'linclass';

return;
