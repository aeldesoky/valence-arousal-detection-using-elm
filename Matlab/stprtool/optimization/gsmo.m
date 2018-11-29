function [x,fval,stat,Nabla] = gsmo(H,f,a,b,LB,UB,x0,Nabla0,options)
% GSMO Generalized SMO algorithm for classifier design.
%
% Synopsis:
%  [x,fval,stat,Nabla] = gsmo(H,f,a,b,LB,UB)
%  [...] = gsmo(H,f,a,b,LB,UB,x0,Nabla0)
%  [...] = gsmo(H,f,a,b,LB,UB,[],[],options)
%  [...] = gsmo(H,f,a,b,LB,UB,x0,Nabla0,options)
% 
% Description:
%  This function implements the generalized SMO algorithm which solves 
%   the following QP task:
%
%   min Q_P(x) = 0.5*x'*H*x + f'*x  
%    x                                      
%
%   s.t.    a'*x = b 
%           LB(i) <= x(i) <= UB(i)   for all i=1:n
%
%  Reference:
%   S.S.Keerthi, E.G.Gilbert: Convergence of a generalized SMO algorithm for SVM 
%   classifier design. Machine Learning, Vol. 46, 2002, pp. 351-360.
%
% Input:
%  H [n x n] Symmetric positive semidefinite matrix.
%  f [n x 1] Vector.
%  a [n x 1] Vector which must not contain zero entries.
%  b [1 x 1] Scalar.
%  LB [n x 1] Lower bound; -inf is allowed.
%  UB [n x 1] Upper bound; inf is allowed.
%  x0 [n x 1] Initial solution.
%  Nabla0 [n x 1] Nabla0 = H*x0 + f.
%  options [struct] 
%    .tolKKT [1 x 1] Determines relaxed KKT conditions (default tolKKT=0.001);
%       it correspondes to $\tau$ in Keerthi's paper.
%    .verb [1 x 1] if > 0 then prints info every verb-th iterations (default 0)
%    .tmax [1 x 1] Maximal number of iterations (default inf).
%  
% Output:
%  x [n x 1] Solution vector.
%  fval [1 x 1] Atteined value of the optimized QP criterion fval=Q_P(x);
%  stat [struct]
%   .exitflag [1 x 1] Indicates which stopping condition was used:
%      relaxed KKT conditions satisfied  ->  exitflag = 1  
%      t >= tmax                         ->  exitflag = 0
%   .t [1x1] Number of iterations.
%   .access [1x1] Access to entries of the matrix H.
%  Nabla [dim x 1] Nabla = H*x+f.
%
% Example:
%  n=50; X = rand(n,n); H = X'*X; f = -10*rand(n,1); 
%  a=rand(n,1)+0.1; b=rand; tmp = rand(n,1); LB = tmp-1; UB = tmp+1;
%  tic; [x1,fval1] = gsmo(H,f,a,b,LB,UB); fval1, toc
%  tic; [x2,fval2] = quadprog(H,f,[],[],a',b,LB,UB); fval2, toc
% 
% See also 
%


% Modifications:
% 29-nov-2006, VF

% cputimer
tic;

% options
if nargin < 6, error('Minimal number of input arguments is six.'); end
if nargin < 7 || isempty(x0), 
  % find some feasible x0
  x0 = zeros(length(f),1); 
  xa = 0; i = 0;
  while x0'*a ~= b,
    i = i + 1;
    if i > length(a),
      error('Constraints cannot be satisfied.\n');
    end
    x0(i) = min(UB(i),max(LB(i),(b-xa)/a(i)));
    xa = xa + x0(i)*a(i);
  end
end
if nargin < 8 || isempty(Nabla0), Nabla0 = H*x0+f; end
if nargin < 9, options = []; else options = c2s(options); end
if ~isfield(options,'tolKKT'), options.tolKKT = 0.001; end
if ~isfield(options,'tmax'), options.tmax = inf; end
if ~isfield(options,'verb'), options.verb = 0; end

[x,exitflag,t,access,Nabla] = gsmo_mex(H,f,a,b,LB,UB,x0,Nabla0,...
                                         options.tmax,options.tolKKT,options.verb);

fval = 0.5*x'*(Nabla + f);

stat.t = t;
stat.exitflag = exitflag;
stat.access = access;
stat.cputime = toc;

return;
% EOF
