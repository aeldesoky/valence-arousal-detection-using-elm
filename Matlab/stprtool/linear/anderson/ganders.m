function model = ganders( distrib, options, init_model )
% GANDERS Solves the Generalized Anderson's task.
%
% Synopsis:
%  model = ganders( distrib)
%  model = ganders( distrib, options)
%  model = ganders( distrib, options, init_model )
%
% Description:
%  This function is an implementation of the general framework 
%  to find the optimal solution of the Generalized Anderson's 
%  task  [SH10].
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
%  distrib [struct] Set of binary labeled Gaussians.
%   .Mean [dim x ncomp] Mean vectors.
%   .Cov [dim x dim x ncomp] Covariance matrices.
%   .y [1 x ncomp] Labels of the Gaussians (1 or 2).
% 
%  options [struct] Determines stopping conditions:
%   .tmax [1x1] Maximal number of iterations (default inf).
%   .eps [1x1] Minimal improvement of the optimized 
%     criterion (default 1e-6).
%   .mineps_tmax [1x1] Number of iterations of the one-dimensional 
%     numerical search (default 100).
%
%  init_model [struct] Initial model:
%    .W, .b, .t see below.
%
% Output:
%  model [struct] Binary linear classifer:
%   .W [dim x 1] Normal vector of the found hyperplane W'*x + b = 0.
%   .b [1x1] Bias of the hyperplane.
% 
%   .r [1x1] Mahalanobis distance for the cloasest Gaussian.
%   .err [1x1] Probability of misclassification.
%   .t [1x1] Number of iterations.
%   .exitflag [1x1] 0 ... maximal number of iterations was exceeded.
%                   1 ... solution was found.
%                  -1 ... solution (with err < 0.5) does not exist.
%
% Example:
%  distrib = load('mars');
%  model = ganders( distrib );
%  figure; pandr( model, distrib );
%
% See also 
%  ANDRORIG, EANDERS, GGRADANDR, ANDRERR, LINCLASS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 4-may-2004, VF
% 17-sep-2003, VF

if nargin < 2, options = []; else options=c2s(options); end
if ~isfield(options,'eps'), options.eps=1e-6; end
if ~isfield(options,'tmax'), options.tmax=inf; end
if ~isfield(options,'mineps_tmax'), options.mineps_tmax = 100; end

% get dimension and number of distributions
[dim,ncomp] = size(distrib.Mean);

% inicialization
exitflag = 0;
t = 0;

% add one constant coordinate
Mu = [distrib.Mean;ones(1,ncomp)];
Mu(:,find( distrib.y==2) ) = -Mu(:,find( distrib.y==2) );
C = zeros(dim+1,dim+1,ncomp);
C(1:dim,1:dim,:) = distrib.Cov;

if nargin == 3,
  W = [init_model.W; init_model.b];
  if isfield( init_model, 't' ), t = init_model.t; end
end

if t==0,
  [W,eflag] = optimal_hyperplane(Mu);

  if eflag <= 0, exitflag = -1; end
end
    
% find the minimal radius of all the ellipsoids
[minr,inx] = min_radius(W,Mu,C);
old_minr=minr;

% main cycle
while exitflag==0 & t < options.tmax,
  t = t + 1;
  
  % compute contact points 
  X0=zeros(dim+1,ncomp);
  for i=1:ncomp,
    X0(:,i)=Mu(:,i)-(minr/sqrt(W'*C(:,:,i)*W))*C(:,:,i)*W;
  end

  % find inprovig direction 
  [dW,eflag]=optimal_hyperplane( X0 );
   
  if eflag > 0,

    % find how much to move along the dw direction and move
    k = gat1dsearch(Mu,C,W,dW,options.mineps_tmax,0);

    % update
    W=W*(1-k)+dW*k;

    [minr,inx] = min_radius(W,Mu,C);

    if minr-old_minr < options.eps,
      exitflag = 1;
    end
    
  else
    exitflag = 1;
  end

  old_minr = minr;
end

% setup model
model.W = W(1:end-1);
model.b = W(end);
model.r = minr;
model.err = 1-cdf('norm',minr,0,1);
model.t = t;
model.exitflag = exitflag;
model.options = options;
model.fun = 'linclass';

return;


%----------------------------------------------------------
function [W,exitflag] = optimal_hyperplane(X)
% finds the optimal hyperplane passing through the origin

[dim,num_data]=size(X);

H=eye(dim);
f=zeros(dim,1);
b=-ones(num_data,1);
A=-X';

% quadratic programming
options=optimset('Display','off','Diagnostics','off','LargeScale','off');
[W,fval,exitflag]=quadprog(H,f,A,b,[],[],[],[],[],options);

return;

%----------------------------------------------------------
function [minr,inx] = min_radius(W,Mu,C);
% find radius of the miniaml ellipsoids

Radius = zeros(size(Mu,2),1);
for i = 1:size(Mu,2),
  den = sqrt(W'*C(:,:,i)*W);
  if den ~= 0,
    Radius(i) = W'*Mu(:,i)/sqrt(W'*C(:,:,i)*W);
  else
    Radius(i) = 0;
  end
end
[minr,inx]=min( Radius );

return;
