function model = minball(X,options)
% MINBALL Minimal enclosing ball in kernel feature space. 
%
% Synopsis:
%  model = minball(X)
%  model = minball(X,options)
%
% Description:
%  It computes center and radius of the minimal ball
%  enclosing data X mapped into a feature space induced 
%  by a given kernel. The problem leads to a special instance 
%  of the Quadratic Programming task which is solved by the 
%  GMNP solver (see 'help gmnp').
%
% Input:
%  X [dim x num_data] Input data.
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier (default 'linear'). See 'help kernel'.
%   .arg [1 x nargs] Kernel arguments.
%   .solver [string] Solver to be used (see 'help gmnp'); default 'imdm';
%   .C [1x1] Regularization constant (default []); If C > 0 it is equivalent 
%     to the Support Vector Data Description (or 1-class SVM) by Tax-Duin 
%     with quadratric penalization of overlapping data.
%
% Output:
%  model [struct] Center of the ball in the kernel feature space:
%   .sv.X [dim x nsv] Data determining the center.
%   .Alpha [nsv x 1] Data weights.
%   .r [1x1] Radius of the minimal enclosing ball.
%   .b [1x1] Squared norm of the center equal to Alpha'*K*Alpha.
%   .options [struct] Copy of used options.
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
%  options = struct('ker','rbf','arg',1);
%  model = minball(data.X,options);
%  [Ax,Ay] = meshgrid(linspace(-5,5,100),linspace(-5,5,100));
%  dist = kdist([Ax(:)';Ay(:)'],model);
%  figure; hold on; 
%  ppatterns(data.X); ppatterns(model.sv.X,'ro',12);
%  contour( Ax, Ay, reshape(dist,100,100),[model.r model.r]);
%
% See also 
%  KDIST.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2005, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 24-july-2008, VF: fixed problem with computing r; pointed out by Daewon Lee (MPI, Tuebingen)
% 09-nov-2006, VF, added C; requested by Hsiung, Chang
% 24-jan-2005, VF, Fast GMNP solver used.
% 25-aug-2004, VF, added model.fun = 'kdist' and .diag_add changed to .mu 
% 16-may-2004, VF
% 15-jun-2002, VF

% process input arguments
%-----------------------------------------
if nargin < 2, options = []; else options=c2s(options); end
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'solver'), options.solver = 'imdm'; end
if ~isfield(options,'C'), options.C = []; end


[dim,num_data] = size(X);

% set up QP problem
%-----------------------------------------
K = kernel( X, options.ker, options.arg );
f = -diag(K);

if isempty(options.C), H=2*K; else H=2*K+eye(size(K))/(2*options.C); end

% call GMNP solver
%----------------------------
[Alpha,fval,stat] = gmnp(H,f,options);

% take non-zero Alpha's
%---------------------
inx= find(Alpha > 0);
model.Alpha = Alpha(inx);

% compute radius
%---------------------
K = K(inx,inx);
model.b = model.Alpha'*K*model.Alpha;
% model.r = sum( sqrt( diag(K) - 2*K*model.Alpha + model.b ))/length(inx);
if isempty(options.C)
    model.r = sum( sqrt( diag(K) - 2*K*model.Alpha + model.b ))/length(inx);
else
    model.r = sum( sqrt( diag(K) - 2*K*model.Alpha + model.b - model.Alpha/(2*options.C) ))/length(inx);
end

% setup model
%---------------------
model.sv.X= X(:,inx);
model.sv.inx = inx;
model.nsv = length(inx);
model.options=options;
model.stat = stat;
model.fun = 'kdist';

return;
% EOF
