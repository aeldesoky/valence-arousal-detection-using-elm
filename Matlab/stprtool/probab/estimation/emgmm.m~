function model=emgmm(X,options,init_model)
% EMGMM Expectation-Maximization Algorithm for Gaussian mixture model.
% 
% Synopsis:
%  model = emgmm(X)
%  model = emgmm(X,options)
%  model = emgmm(X,options,init_model)
%
% Description:
%  This function implements the Expectation-Maximization algorithm 
%  (EM) [Schles68][DLR77] which computes the maximum-likelihood 
%  estimate of the paramaters of the Gaussian mixture model (GMM). 
%  The EM algorithm is an iterative procedure which monotonically 
%  increases log-likelihood of the current estimate until it reaches 
%  a local optimum. 
%
%  The number of components of the GMM is given in options.ncomp 
%  (default 2).
%
%  The following three stopping are condition used:
%   1. Improvement of the log-likelihood is less than given
%      threshold
%                logL(t+1)  - logL(t) < options.eps_logL
%   2. Change of the squared differences of a estimated posteriory 
%      probabilities is less than given threshold
%               ||alpha(t+1) - alpha(t)||^2 < options.eps_alpha
%   3. Number of iterations exceeds given threshold.
%               t >= options.tmax 
%
%  The type of estimated covariance matrices is optional:
%    options.cov_type = 'full'      full covariance matrix (default)
%    options.cov_type = 'diag'      diagonal covarinace matrix
%    cov_options.type = 'spherical' spherical covariance matrix
%
%  The initial model (estimate) is selected:
%    1. randomly (options.init = 'random') 
%    2. using C-means (options.init = 'cmeans')
%    3. using the user specified init_model.
%
% Input:
%  X [dim x num_data] Data sample.
%  
%  options [struct] Control paramaters:
%   .ncomp [1x1] Number of components of GMM (default 2).
%   .tmax [1x1] Maximal number of iterations (default inf).
%   .eps_logL [1x1] Minimal improvement in log-likelihood (default 0).
%   .eps_alpha [1x1] Minimal change of Alphas (default 0).
%   .cov_type [1x1] Type of estimated covarince matrices (see above).
%   .init [string] 'random' use random initial model (default);
%                  'cmeans' use K-means to find initial model.
%   .verb [1x1] If 1 then info is displayed (default 0).
% 
%  init_model [struct] Initial model:
%   .Mean [dim x ncomp] Mean vectors.
%   .Cov [dim x dim x ncomp] Covariance matrices.
%   .Priors [1 x ncomp] Weights of mixture components.
%   .Alpha [ncomp x num_data] (optional) Distribution of hidden state.
%   .t [1x1] (optional) Counter of iterations.
%
% Output:
%  model [struct] Estimated Gaussian mixture model:
%   .Mean [dim x ncomp] Mean vectors.
%   .Cov [dim x dim x ncomp] Covariance matrices.
%   .Prior [1 x ncomp] Weights of mixture components.
%   .t [1x1] Number iterations.
%   .options [struct] Copy of used options.
%   .exitflag [int] 0      ... maximal number of iterations was exceeded.
%                   1 or 2 ... EM has converged; indicates which stopping 
%                              was used (see above).
%  
% Example:
% Note: if EM algorithm does not converge run it again from different
% initial model.
%
% EM is used to estimate parameters of mixture of 2 Guassians:
%  true_model = struct('Mean',[-2 2],'Cov',[1 0.5],'Prior',[0.4 0.6]);
%  sample = gmmsamp(true_model, 100);
%  estimated_model = emgmm(sample.X,struct('ncomp',2,'verb',1));
%
%  figure; ppatterns(sample.X);
%  h1=pgmm(true_model,struct('color','r'));
%  h2=pgmm(estimated_model,struct('color','b'));
%  legend([h1(1) h2(1)],'Ground truth', 'ML estimation'); 
%  figure; hold on; xlabel('iterations'); ylabel('log-likelihood');
%  plot( estimated_model.logL );
%
% See also 
%  MLCGMM, MMGAUSS, PDFGMM, GMMSAMP.
%

% About: Statistical Patte7rn Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 26-may-2004, VF, initialization by K-means added
% 1-may-2004, VF
% 19-sep-2003, VF
% 16-mar-2003, VF


% processing input arguments 
% -----------------------------------------
if nargin < 2, options=[]; else options=c2s(options); end

if ~isfield( options, 'ncomp'), options.ncomp = 2; end
if ~isfield( options, 'tmax'), options.tmax = inf; end
if ~isfield( options, 'eps_alpha'), options.eps_alpha = 0; end
if ~isfield( options, 'eps_logL'), options.eps_logL = 0; end
if ~isfield( options, 'cov_type'), options.cov_type = 'full'; end
if ~isfield( options, 'init'), options.init = 'random'; end
if ~isfield( options, 'verb'), options.verb = 0; end

[dim,num_data] = size(X);

% setup initial model 
% ---------------------------------

if nargin == 3,

  % take model from input
  %-----------------------------

  model = init_model; 
  if ~isfield(model,'t'), model.t = 0; end
  if ~isfield(model,'Alpha'), 
     model.Alpha=-inf*ones(options.num_gauss,num_data);
  end
  if ~isfield(model,'logL'), model.logL=-inf; end
else
  
  % compute initial model
  %------------------------------------
  switch options.init,
    % random model
    case 'random' 
     % takes randomly first num_gauss trn. vectors as mean vectors
     inx = randperm(num_data);  
     inx=inx(1:options.ncomp);
     centers_X = X(:,inx);
  
    % K-means clustering
    case 'cmeans'
     tmp = cmeans( X, options.ncomp );
     centers_X = tmp.X;

    otherwise
     error('Unknown initialization method.');
  end

  knn = knnrule({'X',centers_X,'y',[1:options.ncomp]},1);
  y = knnclass(X,knn);

  % uses ML estimation of complete data
  model = mlcgmm( {'X',X,'y',y}, options.cov_type );

  model.Alpha = zeros(options.ncomp,num_data);
  for i = 1:options.ncomp,
    model.Alpha(i,find(y==i)) = 1;
  end
  model.logL= -inf;
  model.t = 1;
  model.options = options;
  model.fun = 'pdfgmm';
  
end


% Main loop of EM algorithm 
% -------------------------------------
model.exitflag = 0;
while model.exitflag == 0 & model.t < options.tmax,

  % counter of iterations
  model.t = model.t + 1;
  
  %----------------------------------------------------
  % E-Step
  % The distribution of hidden states is computed based
  % on the current estimate.
  %----------------------------------------------------

  newAlpha = (model.Prior(:)*ones(1,num_data)).*pdfgauss(X, model);
  newLogL = sum(log(sum(newAlpha,1)));  
  newAlpha = newAlpha./(ones(options.ncomp,1)*sum(newAlpha,1));

  %------------------------------------------------------
  % Stopping conditions.
  %------------------------------------------------------
  
  % 1) change in distribution of hidden state Alpha
  model.delta_alpha = sum(sum((model.Alpha - newAlpha).^2));
  
  % 2) change in log-Likelihood
  model.delta_logL = newLogL - model.logL(end);
  model.logL = [model.logL newLogL];
  
  if options.verb,
    fprintf('%d: logL=%f, delta_logL=%f, delta_alpha=%f\n',...
        model.t, model.logL(end), model.delta_logL, model.delta_alpha );
  end

  if options.eps_logL >= model.delta_logL,
    model.exitflag = 1;
  elseif options.eps_alpha >= model.delta_alpha,
    model.exitflag = 2;
  else

    model.Alpha = newAlpha;
  
    %----------------------------------------------------
    % M-Step
    % The new parameters maximizing expectation of 
    % log-likelihood are computed.
    %----------------------------------------------------

    tmp_model = melgmm(X,model.Alpha,options.cov_type);
    
    model.Mean = tmp_model.Mean;
    model.Cov = tmp_model.Cov;
    model.Prior = tmp_model.Prior;
  
  end
end % while main loop

return;


