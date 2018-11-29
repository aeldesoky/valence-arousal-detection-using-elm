function varargout=mlsigmoid(arg1,arg2,arg3)
% MLSIGMOID Fitting a sigmoid function using ML estimation.
%
% Synopsis:
%  model = mlsigmoid(data,options)
%
% Description:
%  model = mlsigmoid(data,options) computes Maximum-Likelihood
%   estimation of parameters of sigmoid function [Platt99a]
%     p(y==1|x) = 1/(1+exp(A(1)*x+A(2))),
%
%   used to describe a posteriory probability of a hidden binary 
%   state y from {1,2}. The conditional probabilities p(x|y) are 
%   assumed  to be uni-variate Gaussian distribution. The training 
%   samples {(X(1),y(1)),...,(X(num_data),y(num_data))} assumed to 
%   be i.i.d. are given in data.X and data.y.
%
% Input:
%  data [struct] Input sample:
%   .X [1 x num_data] Values of discriminant function.
%   .y [1 x num_data] Corresponding class label (1 or 2).
%
%  options [struct] Control parameters:
%   .regul [1x1] If 1 then fitting is regularized to prevent 
%      overfitting (default 1).
%   .verb [1x1] If 1 then progress info is displayed (default 0).
%
% Output:
%  model.A [1x2] Parameters of sigmoid function.
%  model.logl [1x1] Value of the log-likelihood criterion.
%
% Example:
%  help demo_svmpout;
%
% See also 
%  SIGMOID.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 28-apr-2008, VF; fixed incorrect regularization for the positive labels
% 03-jun-2004, VF
% 11-oct-2003, VF
% 20-sep-2003, VF
% 08-may-2003, VF

if nargin > 2,
  % evaluates log-likelihood (objective function)
  [L,grad] = sigmoidlogl(arg1,arg2,arg3);
  varargout{1} = L;
  varargout{2} = grad;
  return;
end
  
% process inputs 
%-------------------------------------------------------
data = c2s(arg1);
if nargin == 1, options=[]; else options=c2s(arg2); end
if ~isfield(options,'verb'), options.verb=0; end
if ~isfield(options,'regul'), options.regul=1; end

% values of discriminant function
%-------------------------------------------------------
outs = data.X(:);

% targets 
%-------------------------------------------------------
N1=length(find(data.y==1));
N2=length(find(data.y==2));
if options.regul,
  % use regularization
  
  T1=(N1+1)/(N1+2);  % corrected 28. apr 2008; pointed out by Stijn Vanderlooy
  T2=1/(N2+2);
else
  T1=1; 
  T2=0;
end

targets=zeros(N1+N2,1);
targets(find(data.y==1))=T1;
targets(find(data.y==2))=T2;

% set options of Matlab optimizer
%-------------------------------------------------------
if options.verb, 
  opt=optimset('Display','on','GradObj','on');
else
  opt=optimset('Display','off','GradObj','on');
end

% run optimizer to maximize the log-likelihood
%------------------------------------------------------------
x0=[1 1];
[model.A,model.logl] = fminunc('mlsigmoid',x0,opt,targets,outs);

model.fun = 'sigmoid';
varargout{1} = model;

return;


%=======================================================
function [L,grad] = sigmoidlogl(A,targets,outs)
% SIGMOIDLOGL Returns log-likelihood of sigmoid model.
%
% [L,grad] = sigmoidlogl(A,targets,outs)
%
% Description:
%  It evaluates log-likelihood function
%   L = -sum( targets(i)*log(p_i)+(1-targets(i))*log(1-p_i)),
%
% where p_i = 1/(1+exp(A(1)*outs(i)+A(2))) is a sigmoid function.
%

tmp=exp(A(1)*outs+A(2));

p=1./(1+tmp); 

% prevents dividing by 0
inx=find(p==0); p(inx)=1e-12;
inx=find(p==1); p(inx)=1-1e-12;

L = - sum( targets.*log(p)+(1-targets).*log(1-p));

grad(1)=sum(targets.*outs+(outs.*tmp)./(1+tmp) - outs);
grad(2)=sum(targets + tmp./(1+tmp ) - 1);

return;
% EOF