function model = smo( data, options, init_model)
% SMO Sequential Minimal Optimization for binary SVM with L1-soft margin.
%
% Synopsis:
%  model = smo( data )
%  model = smo( data, options )
%  model = smo( data, options, init_model)
%
% Description:
%  This function is implementation of the Sequential Minimal 
%  Optimizer (SMO) [Platt98] to train the binary Support Vector 
%  Machines Classifier with L1-soft margin.
%           
% Input:
%  data [struct] Binary labeled training vectors:
%   .X [dim x num_data] Training vectors.
%   .y [a x num_data] Labels (1 or 2).
%
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier (default 'linear'); 
%     See 'help kernel'for more info.
%   .arg [1 x nargs] Kernel argument(s) (default 1).
%   .C Regularization constant (default C=inf). The constant C can 
%     be given as:
%      C [1x1] .. for all data.
%      C [1x2] .. for each class separately C=[C1,C2].
%      C [1xnum_data] .. for each training vector separately.
%   .eps [1x1] SMO paramater (default 0.001).
%   .tol [1x1] Tolerance of KKT-conditions (default 0.001).
%  
%  init_model [struct] Specifies initial model:
%   .Alpha [num_data x 1] Initial model. 
%   .b [1x1] Bias.
%  If not given then it is set to zero by default.
%
% Output:
%  model [struct] Binary SVM classifier:
%   .Alpha [nsv x 1] Weights (Lagrangians).
%   .b [1x1] Bias.
%   .sv.X [dim x nsv] Support vectors.
%   .nsv [1x1] Number of Support Vectors.
%   .kercnt [1x1] Number of kernel evaluations used by SMO.
%   .trnerr [1x1] Training classification error.
%   .margin [1x1] Margin of the found classifier.
%   .cputime [1x1] Used CPU time in seconds.
%   .options [struct] Copy of used options.
%
% Example:
%  trn = load('riply_trn');  
%  model = smo(trn,struct('ker','rbf','C',10,'arg',1));
%  figure; ppatterns(trn); psvm(model); 
%  tst = load('riply_tst');
%  ypred = svmclass( tst.X, model );
%  cerror( ypred, tst.y )
%
% See also 
%  SVMCLASS, SVMLIGHT, SVMQUADPROG.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 23-may-2004, VF
% 17-September-2001, V. Franc, created

% timer
tic;

% Input arguments 
%-------------------------------------------------------
if nargin < 2,  options = []; else options=c2s(options); end
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'C'), options.C = inf; end
if ~isfield(options,'eps'), options.eps = 0.001; end
if ~isfield(options,'tol'), options.tol = 0.001; end

[dim,num_data] = size(data.X);
if nargin < 3,
  init_model.Alpha = zeros(num_data,1); 
  init_model.b = 0;
end

% run optimizer
%----------------------------------------------------
[model.Alpha, model.b, model.nsv, model.kercnt, model.trnerr, model.margin]...
   = smo_mex(data.X, data.y, options.ker, options.arg, options.C, ...
     options.eps, options.tol, init_model.Alpha, init_model.b );

% set up output
%------------------------------------------------------
inx = find( model.Alpha );
model.sv.X = data.X(:,inx);
model.sv.y = data.y(inx);
model.sv.inx = inx;
model.Alpha = model.Alpha(inx);
model.Alpha(find(model.sv.y == 2)) = -model.Alpha(find(model.sv.y == 2));

% computes normal vector of the hypeprlane if linear kernel used
if strcmpi(options.ker,'linear'),
  model.W = model.sv.X*model.Alpha;
end

model.options = options;
model.fun = 'svmclass';

model.cputime = toc;

return;
