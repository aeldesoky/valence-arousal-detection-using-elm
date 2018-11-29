function model = svmquadprog(data,options)
% SVMQUADPROG SVM trained by Matlab Optimization Toolbox.
%
% Synopsis:
%  model = svmquadprog( data )
%  model = svmquadprog( data, options )
%
% Description:
%  This function trains binary Support Vector Machines classifer 
%  with L1 or L2-soft margin. The SVM quadratic programming task 
%  is solved by the 'quadprog.m' of the Matlab Optimization toolbox.
%
%  See 'help svmclass' to see how to classify data with found classifier.
%  
% Input:
%  data [struct] Binary labeled training data:
%   .X [dim x num_data] Vectors.
%   .y [1 x num_data] Training labels.
%
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier (default 'linear'). 
%      See 'help kernel' for more info.
%   .arg [1 x nargs] Kernel argument(s).
%   .C SVM regularization constant (default inf):
%     [1 x 1] .. the same for all training vectors.
%     [1 x 2] .. for each class separately C=[C1,C2],
%     [1 x num_data] .. each training vector separately.
%   .norm [1x1] 1 .. L1-soft margin penalization (default).
%               2 .. L2-soft margin penalization.
%
% Output:
%  model [struct] Binary SVM classifier:
%   .Alpha [nsv x 1] Weights.
%   .b [1x1] Bias of the decision function.
%   .sv.X [dim x nsv] Support vectors.
%   .nsv [1x1] Number of support vectors.
%   .kercnt [1x1] Number of used kernel evaluations.
%   .trnerr [1x1] Training classification error.
%   .margin [1x1] Margin of found classifier.
%   .cputime [1x1] Used CPU time in seconds.
%   .options [struct] Copy of used options.
%   .exitflag [1x1] Exitflag of the QUADPROG function. 
%     (if > 0 then it has converged to the solution).
%
% Example:
%  data = load('riply_trn');
%  options = struct('ker','rbf','arg',1,'C',10);
%  model = svmquadprog(data,options)
%  figure; ppatterns(data); psvm(model);
%
% See also 
%  SMO, SVMLIGHT, SVMCLASS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 31-may-2004, VF
% 16-may-2004, VF
% 17-Feb-2003, VF
% 28-Nov-2001, VF, used quadprog instead of qp
% 23-Occt-2001, VF
% 19-September-2001, V. Franc, renamed to svmmot.
% 8-July-2001, V.Franc, comments changed, bias mistake removed.
% 28-April-2001, V.Franc, flps counter added
% 10-April-2001, V. Franc, created

% timer
tic;

% input aruments
%-------------------------------------------
data=c2s(data);
[dim,num_data]=size(data.X);
if nargin < 2, options=[]; else options=c2s(options); end
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'C'), options.C = inf; end
if ~isfield(options,'norm'), options.norm = 1; end
if ~isfield(options,'mu'), options.mu = 1e-12; end
if ~isfield(options,'eps'), options.eps = 1e-12; end

% Set up QP task
%----------------------------

% labels {1,2} -> {1,-1}
y = data.y(:);
y(find(y==2)) = -1;

% compute kernel matrix
H = kernel(data.X,data.X,options.ker,options.arg).*(y*y');

% add small numbers to diagonal 
H = H + options.mu*eye(size(H));

Aeq = y';
beq = 0;

f = -ones(num_data,1);       % Alpha

LB = zeros(num_data,1);      % 0 <= Alpha
x0 = zeros(num_data,1);      % starting point

if options.norm==1,
  % L1-soft margin
  %---------------------
  if length(options.C) == 1,
    UB = options.C*ones(num_data,1);
  elseif length(options.C) == 2,
    UB=zeros(num_data,1);
    UB(find(data.y==1))=options.C(1);
    UB(find(data.y==2))=options.C(2);
  else
    UB=options.C(:);
  end
  vectorC=zeros(num_data,1);
else
  % L2-soft margin
  %---------------------
  UB=ones(num_data,1)*inf;
  vectorC = ones(num_data,1);
  if length(options.C) == 1,
    vectorC = vectorC*options.C;
  elseif length(options.C) == 2,
    inx1=find(data.y==1);inx2=find(data.y==2);
    vectorC(inx1)=options.C(1);
    vectorC(inx2)=options.C(2);
  else
    vectorC = options.C(:);
  end
  vectorC = 1./(2*vectorC);
  H = H + diag(vectorC);
end

% call optimization toolbox 
% ----------------------------------
qp_options = optimset('Display','off');
[Alpha,fval,exitflag] = quadprog(H, f, [],[],Aeq, beq, LB, UB, x0, qp_options);

inx_sv = find( Alpha > options.eps);

% compute bias
%--------------------------
% take boundary (f(x)=+/-1) support vectors 0 < Alpha < C
inx_bound = find( Alpha > options.eps & Alpha < (options.C - options.eps));

if length( inx_bound ) ~= 0,
  model.b = sum(y(inx_bound)-H(inx_bound,inx_sv)*...
     Alpha(inx_sv).*y(inx_bound))/length( inx_bound );
else
  disp('Bias cannot be determined.');
  model.b=0;
end

% compute margin
%------------------------
if options.norm == 1
  w2 = Alpha(inx_sv)'*H(inx_sv,inx_sv)*Alpha(inx_sv);
else
  w2 = Alpha(inx_sv)'*(H(inx_sv,inx_sv)-diag(vectorC(inx_sv)))*Alpha(inx_sv);
end

margin = 1/sqrt(w2);


% compute training classification error 
%-------------------------------------------
Alpha = Alpha.*y;
model.Alpha = Alpha( inx_sv );
model.sv.X = data.X(:,inx_sv );
model.sv.y = data.y(inx_sv );
model.sv.inx = inx_sv;
model.nsv = length( inx_sv );
model.margin = margin;
model.exitflag = exitflag;
model.options = options;
model.kercnt = num_data*(num_data+1)/2;
model.trnerr = cerror(data.y,svmclass(data.X, model));
model.fun = 'svmclass';

% used CPU time
model.cputime=toc;

return;
% EOF