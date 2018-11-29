function [x,fval,stat] = qpssvm(H,f,b,I,x0,options)
% QPSSVM Solves QP task required for StructSVM learning.
%
% Synopsis:
%  [x,fval,stat] = qpssvm(H,f,b,I)
%  [x,fval,stat] = qpssvm(H,f,b,I,x0)
%  [x,fval,stat] = qpssvm(H,f,b,I,x0,options)
%  [x,fval,stat] = qpssvm(H,f,b,I,[],options)
% 
% Description:
%  This function solves the following QP task:
%
%   min 0.5*x'*H*x + f'*x
%    x
% subject to 
%   sum(x(find(I==k))) <= b   for all k=1:max(I)
%   x >= 0
%
% where I is a vector of indeices such that unique(I) = 1:max(I)
% and max(I) <= length(H).
%
% Input:
%  H [n x n] Symmetric positive semidefinite matrix.
%  f [n x 1] Vector.
%  b [1 x 1] Scalar.
%  I [n x 1] Vector of indices such that unique(I) = 1:max(I);
%  x0 [n x 1] Initial solution vector.
%  tmax [1 x 1] Maximal number of iterations.
%  tolabs [1 x 1] Absolute tolerance stopping condition. 
%  tolrel [1 x 1] Relative tolerance stopping condition. 
%  verb [1 x 1] if > 0 then prints info every verb-th iterations.
%
% Output:
%  x [n x 1] Solution vector.
%  fval [1 x 1] Value of optimized QP criterion.
%  exitflag [1 x 1] Indicates which stopping condition was used:
%    UB-LB <= tolabs           ->  exitflag = 1   Abs. tolerance.
%    UB-LB <= UB*tolrel        ->  exitflag = 2   Relative tolerance.
%    t >= tmax                 ->  exitflag = 0   Number of iterations.
%   where UB is value of the primal and LB of the dual criterion.
%    
%  t [1x1] Number of iterations.
%  access [1x1] Access to elements of the matrix H.
%  History [2x(t+1)] UB and LB with respect to number of iterations.
%
% Example:
%  I = [1,1,1,2,2,3,3,4,4,4,4,5,5]; 
%  Z = randn(length(I),100); H=Z*Z'; f=randn(length(I),1); b = 1;
%  [x1,fval1,stat1] = qpssvm(H,f,b,I)
%
%  A = zeros(max(I),length(I));
%  for i=1:max(I), A(i,find(I==i)) = 1; end
%  b = b*ones(max(I),1);
%  [x2,fval2] = quadprog(H,f,A,b,[],[],zeros(length(I),1))
%
% See also 
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2006, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 18-feb-2005, VF
% 17-feb-2005, VF

% cputimer
tic;

% options
if nargin < 4, error('Incorrect number of input arguments.'); end
if nargin < 5 || isempty(x0), x0 = zeros(length(f),1); end
if nargin < 6, options = []; else options = c2s(options); end
if ~isfield(options,'tolabs'), options.tolabs = 1e-6; end
if ~isfield(options,'tolrel'), options.tolrel = 1e-6; end
if ~isfield(options,'tmax'), options.tmax = inf; end
if ~isfield(options,'verb'), options.verb = 0; end
if ~isfield(options,'solver'), options.solver = 'sca'; end

[x,exitflag,t,access,History] = qpssvm_mex(H,f,b,uint16(I),x0,...
          options.tmax,options.tolabs,options.tolrel,options.verb);

stat.t = t;
stat.exitflag = exitflag;
stat.LB_History = History(1,:);
stat.UB_History = History(2,:);
stat.LB = History(1,end);
stat.UB = History(2,end);
stat.access = access;
fval = stat.UB;
stat.cputime = toc;

return;
% EOF