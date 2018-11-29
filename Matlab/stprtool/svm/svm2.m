function model = svm2(data,options)
% SVM2 Learning of binary SVM classifier with L2-soft margin.
%
% Synopsis:
%  model = svm2(data)
%  model = svm2(data,options)
%
% Description:
%  This function learns binary Support Vector Machines
%  classifier with L2-soft margin. The corresponding quadratic 
%  programming task is solved by one of the following 
%  algorithms:
%    mdm  ... Mitchell-Demyanov-Malozemov (MDM) algorithm.
%    imdm ... Improved MDM algorithm (IMDM) (defaut).
%
%  For more info refer to V.Franc: Optimization Algorithms for Kernel 
%  Methods. Research report. CTU-CMP-2005-22. CTU FEL Prague. 2005.
%  ftp://cmp.felk.cvut.cz/pub/cmp/articles/franc/Franc-PhD.pdf .
%
% Input:
%  data [struct] Training data:
%   .X [dim x num_data] Training vectors.
%   .y [1 x num_data] Labels must equal 1 and/or 2.
%
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier. See 'help kernel'.
%   .arg [1 x nargs] Kernel argument(s).
%   .C [1x2] Regularization constants for class 1 and 2; 
%      if C is [1x1] then the same C is used for both classes.
%   .solver [string] Solver to be used: 'mdm', 'imdm' (default).
%   .tmax [1x1] Maximal number of iterations (default inf).
%   .tolabs [1x1] Absolute tolerance stopping condition (default 0.0).
%   .tolrel [1x1] Relative tolerance stopping condition (default 1e-3).
%   .thlb [1x1] Threshold on lower bound (default inf).
%   .cache [1x1] #of columns of kernel matrix to be cached (default 1000).
%   .verb [1x1] If > 0 then some info is displayed (default 0).
%
% Output:
%  model [struct] Binary SVM classifier:
%   .Alpha [nsv x 1] Weights of support vectors.
%   .b [1x1] Bias of decision function.
%   .sv.X [dim x nsv] Support vectors.
%   .sv.inx [1 x nsv] Indices of SVs (model.sv.X = data.X(:,inx)).
%   .nsv [int] Number of Support Vectors.
%   .kercnt [1x1] Number of kernel evaluations.
%   .trnerr [1x1] Classification error on training data.
%   .margin [1x1] Margin.
%   .options [struct] Copy of used options.
%   .cputime [1x1] Used CPU time in seconds (meassured by tic-toc).
%   .stat [struct] Statistics about optimization:
%     .access [1x1] Number of requested columns of matrix H.
%     .t [1x1] Number of iterations.
%     .UB [1x1] Upper bound on the optimal value of criterion. 
%     .LB [1x1] Lower bound on the optimal value of criterion. 
%     .LB_History [1x(t+1)] LB with respect to iteration.
%     .UB_History [1x(t+1)] UB with respect to iteration.
%     .NA [1x1] Number of non-zero entries in solution.
%
% Example:
%  data = load('riply_trn');
%  options = struct('ker','rbf','arg',1,'C',1);
%  model = svm2(data,options )
%  figure; ppatterns(data); psvm( model );
%
% See also
%  SVMCLASS, SVMLIGHT, SMO, GNPP.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2005, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 07-sep-2007, VF, now it is possible to use distinct reg. constants for both classes
% 09-sep-2005, VF
% 08-aug-2005, VF
% 24-jan-2005, VF
% 29-nov-2004, VF

% restart clock
tic;

if nargin < 2, options = []; else options = c2s(options); end
if ~isfield(options,'solver'), options.solver = 'imdm'; end
if ~isfield(options,'tolabs'), options.tolabs = 0; end
if ~isfield(options,'tolrel'), options.tolrel = 1e-3; end
if ~isfield(options,'thlb'), options.thlb = inf; end
if ~isfield(options,'tmax'), options.tmax = inf; end
if ~isfield(options,'C'), options.C = inf; end
if length(options.C)==1, options.C = [1 1]*options.C; end % C [1 x 2] 
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'cache'), options.cache = 1000; end
if ~isfield(options,'verb'), options.verb = 0; end

% call MEX implementation of QPC2 solver
[Alpha,b,exitflag,kercnt,access,errcnt,t,UB,LB,History] = svm2_mex(...
    data.X,...
    data.y,...
    options.ker,...
    options.arg,...
    options.C,...
    options.solver,...
    options.tmax,...
    options.tolabs, ...
    options.tolrel,...
    options.thlb, ...
    options.cache, ...
    options.verb );

% remove non-support vectors
inx = find(Alpha ~=0 );

% setup output model
model.Alpha = Alpha(inx);
model.b = b;
model.sv.X = data.X(:,inx);
model.sv.inx = inx;
model.sv.y = data.y(inx);
model.nsv = length(inx);
if strcmp( options.ker, 'linear'),
  model.W = model.sv.X * model.Alpha;
end
model.options = options;
model.kercnt = kercnt;
model.trnerr = errcnt/size(data.X,2);
model.errcnt = errcnt;
%model.margin = 1/sqrt(sum(abs(model.Alpha))-sum(model.Alpha.^2)...
%    /2/model.options.C);
model.exitflag = exitflag;
model.stat.access = access;
model.stat.t = t;
model.stat.UB = UB;
model.stat.LB = LB;
model.stat.LB_History = History(1,:);
model.stat.UB_History = History(2,:);
model.stat.NA = length(inx);
model.cputime = toc;
model.fun = 'svmclass';

return;
%EOF

