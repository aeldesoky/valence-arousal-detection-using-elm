function model = kperceptr(data,options)
% KPERCEPTR Kernel Perceptron.
%
% Synopsis:
%  model = kperceptr(data)
%  model = kperceptr(data,options)
%
% Description:
%  This function is an implementation of the kernel version
%  of the Perceptron algorithm. The kernel perceptron search 
%  for the kernel binary classifier with zero emprical error.
%
% Input:
%  data [struct] Binary labeled training data:
%   .X [dim x num_data] Vectors.
%   .y [1 x num_data] Labels (1 or 2).
%  
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier (default 'linear').
%     See 'help kernel' for more info.
%   .arg [1 x nargs] Kernel argument.
%   .tmax [1x1] Maximal number of iterations (default inf).
%
% Output:
%  model [struct] Found kernel classifer:
%   .Alpha [nsv x 1] Multipliers of the training data.
%   .b [1x1] Bias of the decision rule.
%   .sv.X [dim x nsv] Training data with non-zero Alphas.
%   .exitflag [1x1] 1 ... Perceptron has converged.
%                   0 ... Maximal number of iterations exceeded.
%   .iter [1x1] Number of iterations.
%   .kercnt [1x1] Number of kernel evaluations.
%   .trnerr [1x1] Training classification error; Note: if exitflag==1 
%     then trnerr = 0.
%   .options [struct] Copy of options.
%   .cputime [real] Used cputime in seconds.
%  
%  If the linear kernel is used then model.W [dim x 1] contains 
%  normal vector of the separating hyperplane.
%  
% Example:
%  data = load('vltava');
%  model = kperceptr(data, struct('ker','poly','arg',2));
%  figure; ppatterns(data); pboundary(model);
% 
% See also SVMCLASS, SVM.
%

% Modifications: 
% 10-may-2004, VF
% 18-July-2003, VF
% 21-Nov-2001, V. Franc

tic;

% process inputs
%=======================================
if nargin < 2, options=[]; else options=c2s(options); end
if ~isfield( options, 'ker'), options.ker = 'linear'; end
if ~isfield( options, 'arg'), options.arg = 1; end
if ~isfield( options, 'tmax'), options.tmax = inf; end

[dim, num_data ] = size( data.X );

% change {1,2} --> {1,-1}
y = data.y;
y(find(y==2)) = -1;

% inicialize multiliers Alpha and bias
% =========================================
Alpha = zeros( num_data,1 );
dfce = zeros(num_data,1);  % cache for f(x_i)=<Phi(x_i),w>+b
b = 0;
exitflag = 0;
iter = 0;
kercnt = 0;

% Main loop
%====================================
while iter < options.tmax & exitflag == 0,
  
  iter = iter + 1;
  
  [min_dfce, inx ] = min( dfce.*y(:) );
  
  if min_dfce <=0 ,
    
    % Perceptron rule in terms of dot products
    old_alpha = Alpha(inx);
    old_b = b;
    Alpha(inx) = old_alpha + y(inx);
    b = old_b + y(inx);
    
    % Updates cache
    k_inx = kernel(data.X,data.X(:,inx), options.ker, options.arg );
    kercnt = kercnt + num_data;
    
    dfce = dfce + (Alpha(inx)-old_alpha)*k_inx + b - old_b;
    
  else
    exitflag = 1;
    
    % scales f(x) such that f(x)=1 for the patterns closest to
    % the separating hyperplane
    Alpha = Alpha/min_dfce;
    b = b/min_dfce;
  end
  
end

% fill up structure model
%=================================
inx = find( Alpha ~= 0);
Alpha = Alpha(inx);

model.Alpha = Alpha(:);
model.b = b;
model.sv.X = data.X(:,inx);
model.sv.y = data.y(inx);
model.sv.inx = inx;
model.nsv = length( inx);
model.options = options;
model.iter = iter;
model.kercnt = kercnt;
model.trnerr = length( find(dfce.*y(:) < 0))/num_data;
model.exitflag = exitflag;
model.fun = 'svmclass';

if strcmpi('linear',options.ker),
  model.W = model.sv.X * model.Alpha;
end

model.cputime = toc;

return;
