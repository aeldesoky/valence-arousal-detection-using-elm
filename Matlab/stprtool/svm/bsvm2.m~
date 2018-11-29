function model = bsvm2( data, options )
% BSVM2 Multi-class BSVM with L2-soft margin.
%
% Synopsis:
%  model = bsvm2( data, options ) 
%
% Description:
%  This function trains the multi-class SVM classifier based
%  on BSVM formulation (bias added to the objective function) and
%  L2-soft margin penalization of misclassifications.
%  The quadratic programming task is optimized by one of the
%  following algorithms:
%    mdm  ... Mitchell-Demyanov-Malozemov
%    imdm  ... Mitchell-Demyanov-Malozemov Improved 1.
%    iimdm  ... Mitchell-Demyanov-Malozemov Improved 2.
%    kozinec ... Kozinec algorithm.
%    keerthi ... NPA algorithm by Keerthi et al.
%    kowalczyk ... Based on Kowalczyk's maximal margin perceptron.
%
%  For more info refer to V.Franc: Optimization Algorithms for Kernel 
%  Methods. Research report. CTU-CMP-2005-22. CTU FEL Prague. 2005.
%  ftp://cmp.felk.cvut.cz/pub/cmp/articles/franc/Franc-PhD.pdf .
%
% Input:
%  data [struct] Training data:
%   .X [dim x num_data] Training vectors.
%   .y [1 x num_data] Labels (1,2,...,nclass).
%
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier. See 'help kernel'.
%   .arg [1 x nargs] Kernel argument(s).
%   .C [1x1] Regularization constant.
%   .solver [string] Solver to be used: 'mdm', 'imdm' (default), 'iimdm', 
%     'kozinec', 'kowalczyk','keerthi'.
%   .tmax [1x1] Maximal number of iterations (default inf).
%   .tolabs [1x1] Absolute tolerance stopping condition (default 0.0).
%   .tolrel [1x1] Relative tolerance stopping condition (default 0.001).
%   .thlb [1x1] Thereshold on the lower bound (default inf).
%   .cache [1x1] Number of columns of kernel matrix to be cached (default 1000).
%   .verb [1x1] If > 0 then some info is printed (default 0).
%
% Output:
%  model [struct] Multi-class SVM classifier:
%   .Alpha [nsv x nclass] Weights.
%   .b [nclass x 1] Biases.
%   .sv.X [dim x nsv] Support vectors.
%   .nsv [1x1] Number of support vectors.
%   .options [struct] Copy of input options.
%   .t [1x1] Number of iterations.
%   .UB [1x1] Upper bound on the optimal solution.
%   .LB [1x1] Lower bound on the optimal solution.
%   .History [2 x (t+1)] UB and LB with respect to t.
%   .trnerr [1x1] Training classification error.
%   .kercnt [1x1] Number of kernel evaluations.
%   .cputime [1x1] CPU time (measured by tic-toc).
%   .stat [struct] Statistics about optimization:
%     .access [1x1] Number of requested columns of matrix H.
%     .t [1x1] Number of iterations.
%     .UB [1x1] Upper bound on the optimal value of criterion.
%     .LB [1x1] Lower bound on the optimal value of criterion.
%     .LB_History [1x(t+1)] LB with respect to t.
%     .UB_History [1x(t+1)] UB with respect to t.
%     .NA [1x1] Number of non-zero elements in solution.
%
% Example:
%  data = load('pentagon');
%  options = struct('ker','rbf','arg',1,'C',10);
%  model = bsvm2(data,options )
%  figure; 
%  ppatterns(data); ppatterns(model.sv.X,'ok',12);
%  pboundary(model);
%
% See also 
%  SVMCLASS, OAASVM, OAOSVM, GMNP.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2005, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 09-sep-2005, VF
% 24-jan-2005, VF
% 29-nov-2004, VF
% 26-nov-2004, VF
% 16-Nov-2004, VF
% 31-may-2004, VF
% 23-jan-2003, VF

tic;

% process inputs 
%-------------------------------------------------------
data=c2s(data);
if nargin < 2, options=[]; else options=c2s(options); end
if ~isfield(options,'ker'), options.ker='linear'; end
if ~isfield(options,'arg'), options.arg=1; end
if ~isfield(options,'C'), options.C=inf; end
if ~isfield(options,'tmax'), options.tmax=inf; end
if ~isfield(options,'tolabs'), options.tolabs=0; end
if ~isfield(options,'tolrel'), options.tolrel=0.001; end
if ~isfield(options,'thlb'), options.thlb=inf; end
if ~isfield(options,'solver'), options.solver='imdm'; end
if ~isfield(options,'cache'), options.cache = 1000; end
if ~isfield(options,'verb'), options.verb=0; end

[dim,num_data]=size(data.X);
nclass = max(data.y);

% display info
%---------------------
if options.verb > 0,
  fprintf('Binary rules: %d\n', nclass);
  fprintf('Training data: %d\n', num_data);
  fprintf('Dimension: %d\n', dim);
  if isfield( options, 'ker'), fprintf('Kernel: %s\n', options.ker); end
  if isfield( options, 'arg'), fprintf('arg: %f\n', options.arg(1)); end
  if isfield( options, 'C'), fprintf('C: %f\n', options.C); end
  fprintf('Solver: %s\n', options.solver);
end

% call MEX implementation
[Alpha,b,exitflag,kercnt,access,trnerr,t,NA,UB,LB,History] = bsvm2_mex(...
    data.X,...
    data.y,...
    options.ker,...
    options.arg,...
    options.C,...
    options.solver,...
    options.tmax,...
    options.tolabs, ...
    options.tolrel,...
    options.thlb,...
    options.cache, ...
    options.verb );

% set up model
%-------------------------
sv_inx = find( sum(abs(Alpha),1) ~= 0 );
%sv_inx = find( sum(abs(Alpha),1) ~= inf );
Alpha = Alpha(:,sv_inx)';
for i = 1:size(Alpha,2),
  inx = find( data.y(sv_inx) ~= i);
  Alpha(inx,i) = -Alpha(inx,i);
end

model.Alpha = Alpha;
model.b = b;
model.sv.X = data.X(:,sv_inx);
model.sv.y = data.y(sv_inx);
model.sv.inx = sv_inx;
model.nsv = length(sv_inx);
model.options = options;
model.exitflag = exitflag;
model.trnerr = trnerr;
model.kercnt = kercnt;
model.stat.access = access;
model.stat.t = t;
model.stat.UB = UB;
model.stat.LB = LB;
model.stat.LB_History = History(1,:);
model.stat.UB_History = History(2,:);
model.stat.NA = NA;
model.cputime = toc;

if strcmpi('linear',options.ker) == 1,
  model.W = model.sv.X*model.Alpha;
  model.fun = 'linclass';
else
  model.fun = 'svmclass';
end

return;

% EOF
