function model=svm1d(data,options)
% SVM1D Linear SVM for 1-dimensional input data.
%
% Synopsis:
%  model = svm1d( data )
%  model = svm1d( data, options )
%
% Description:
%  model = svm1d( data ) trains the linear SVM binary
%    classifier for the 1-dimensional training data.
%    The optimizer is based on a modification of the 
%    Sequential Minimal Optimizer (SMO) [Platt98]. 
%    The trainined classfier is defined as
%      q(x) = 1 if W*x + b >= 0
%           = 2 if W*x + b < 0
%
%  model = svm1d( data, options ) use to set up control
%    parameters for the SVM and the SMO algorithm.
%
% Input:
%  data [struct] Input 1-dimensional binary labeled training data:
%   .X [1 x num_data] Training numbers.
%   .y [1 x num_data] Labels (1 or 2).
%  
%  options [struct] Control parameters:
%   .C [1x1] SVM regularization constant (default C=inf). 
%   .eps [1x1] Tolerance of KKT-conditions (default eps=0.001).
%   .tol [1x1] Minimal change of variables (default tol=0.001).
%
% Output:
%  model [struct] Found SVM model:
%   .Alpha [nsv x 1] Weights.
%   .b [1x1] Bias of decision function.
%   .sv.X [1 x nsv] Support vectors.
%   .W [1x1] Explicit value of the normal vector (scalar).
%
%   .nsv [1x1] Number of Support Vectors.
%   .kercnt [1x1] Number of kernel evaluations (multiplications 
%     in this 1-d linear case) used by the SMO.
%   .trnerr [1x1] Training classification error.
%   .margin [1x1] Margin of found classifier.
%   .cputime [1x1] Used CPU time in seconds.
%   .options [struct] Copy of used options.
%
% See also 
%  SMO, SVMCLASS, KFD, KFDQP.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 17-may-2004, VF
% 14-may-2004, VF
% 15-july-2003, VF

% timer
tic;

% Process input arguments 
% --------------------------
[dim,num_data] = size(data.X);
if dim ~= 1,
  error('Inpu  data must be one-dimensional.');
end
  
if nargin < 2,  options = []; else options=c2s(options); end
if ~isfield(options,'C'), options.C = inf; end
if ~isfield(options,'eps'), options.eps = 0.001; end
if ~isfield(options,'tol'), options.tol = 0.001; end

% call MEX function
%---------------------------
[model.Alpha, model.b, model.nsv, model.kercnt, model.trnerr, model.margin]...
     = smo1d_mex(data.X, data.y, options.C, options.eps, options.tol);

% fill up the output structure
%---------------------------------
inx = find( model.Alpha );
model.sv.X = data.X(:,inx);
model.sv.y = data.y(inx);
model.sv.inx = inx;
model.Alpha = model.Alpha(inx);
model.Alpha( find(model.sv.y==2)) = -model.Alpha( find(model.sv.y==2 ));
model.W = model.sv.X*model.Alpha;
options.ker = 'linear';
options.arg = 1;
model.options = options;
model.fun = 'svmclass';
model.cputime = toc;

return;
% EOF
