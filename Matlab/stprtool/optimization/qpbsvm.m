function [x,fval,stat,Nabla] = qpbsvm(H,f,UB,x0,Nabla0,options)
% QPBSVM Solves QP task required for learning SVM without bias term.
%
% Synopsis:
%  [x,fval,stat,Nabla] = qpbsvm(H,f,UB)
%  [...] = qpbsvm(H,f,UB,x0,Nabla0)
%  [...] = qpbsvm(H,f,UB,[],[],options)
%  [...] = qpbsvm(H,f,UB,x0,Nabla0,options)
% 
% Description:
%  This function solves the following QP task:
%
%   min Q_P(x) = 0.5*x'*H*x + f'*x     s.t.  0 <= x <= UB
%    x
%
% Input:
%  H [n x n] Symmetric positive semidefinite matrix.
%  f [n x 1] Vector.
%  UB [1 x 1] Scalar.
%  x0 [n x 1] Initial solution (default zeros).
%  Nabla0 [n x 1] Nabla0 = H*x0 + f.
%  options [struct] 
%    .tmax [1 x 1] Maximal number of iterations.
%    .tolabs [1 x 1] Absolute tolerance stopping condition. 
%    .tolrel [1 x 1] Relative tolerance stopping condition. 
%    .verb [1 x 1] if > 0 then prints info every verb-th iterations.
%    .solver [string] 
%     'sca' ... Sequential Coordinatewise Algorithm (Generalize Gauss-Seidel); 
%     'scas' ... Greedy variant - udpate variable yielding the best improvement.
%  
% Output:
%  x [n x 1] Solution vector.
%  fval [1 x 1] Atteined value of the optimized QP criterion fval=Q_P(x);
%  stat [struct]
%   .exitflag [1 x 1] Indicates which stopping condition was used:
%      Q_P-Q_D <= tolabs           ->  exitflag = 1   Abs. tolerance.
%      Q_P-Q_D <= Q_P*tolrel       ->  exitflag = 2   Relative tolerance.
%      t >= tmax                   ->  exitflag = 0   Number of iterations.
%     where Q_P is value of the primal and Q_D of the dual criterion.
%   .t [1x1] Number of iterations.
%   .access [1x1] Access to elements of the matrix H.
%   .History [2x(t+1)] UB and LB with respect to number of iterations.
%  Nabla [dim x 1] Nabla = H*x+f.
%
% Example:
%  n=500; X = rand(n,n); H = X'*X; f = -10*rand(n,1); UB = 1; e = ones(n,1);
%  tic; [x1,fval1] = qpbsvm(H,f,UB); fval1, toc
%  tic; [x2,fval2] = quadprog(H,f,[],[],[],[],0*e,UB*e); fval2, toc
% 
% See also 
%


% Modifications:
% 20-nov-2006, VF

% cputimer
tic;

% options
if nargin < 3, error('Incorrect number of input arguments.'); end
if nargin < 4 || isempty(x0), x0 = zeros(length(f),1); end
if nargin < 5 || isempty(Nabla0), Nabla0 = H*x0+f; end
if nargin < 6, options = []; else options = c2s(options); end
if ~isfield(options,'tolabs'), options.tolabs = 0; end
if ~isfield(options,'tolrel'), options.tolrel = 1e-6; end
if ~isfield(options,'tolKKT'), options.tolKKT = 0; end
if ~isfield(options,'tmax'), options.tmax = inf; end
if ~isfield(options,'verb'), options.verb = 0; end
if ~isfield(options,'solver'), options.solver = 'sca'; end

% [x,exitflag,t,access,History] = qpbsvm_mex(H,f,UB,tmax,tolabs,tolrel,verb,x0,Nabla0)
[x,exitflag,t,access,History,Nabla] = qpbsvm_mex(H,f,UB,...
      options.solver,options.tmax,options.tolabs,options.tolrel,options.tolKKT,options.verb,x0,Nabla0);

stat.t = t;
stat.exitflag = exitflag;
stat.Q_P_History = History(1,:);
stat.Q_D_History = History(2,:);
stat.Q_P = History(1,end);
stat.Q_D = History(2,end);
stat.access = access;
fval = stat.Q_P;
stat.cputime = toc;

return;
% EOF
