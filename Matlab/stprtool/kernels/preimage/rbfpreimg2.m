function z = rbfpreimg2(varargin)
% RBFPREIMG2 RBF pre-image problem by Gradient optimization.
%
% Synopsis:
%  z = rbfpreimg2(model)
%  z = rbfpreimg2(model,options)
%
%  Description:
%   z = rbfpreimg2(model) it uses gradient method to solve 
%     the pre-image problem for the Radial Basis Function (RBF) 
%     kernel. The function 'fminunc' of the Matlab Optimization 
%     toolbox is exploited for 1D search along the gradient 
%     direction.
%
%   z = rbfpreimg2(model,options) use to specify the control
%     parameters of the gradient optimization.
% 
% Input:
%  model [struct] Kernel expansion:
%   .Alpha [num_data x 1] Weight vector.
%   .sv.X [dim x num_data] Vectors determining the kernel expansion.
%   .options.arg [1x1] Argument of the RBF kernel (see 'help kernel').
%
%  options [struct] Control parameters:
%   .min_improvement [1x1] Minimal allowed improvement of the objective 
%     function in two consecutive steps (default 1e-3).
%  options.start_t [1x1] Starting value of the 1D search procedure 
%   (default 1e-3).
%
% Output:
%  z [dim x 1] Found preimage.
%
% See also 
%  RBFPREIMG, RBFPREIMG3, RSRBF, KPCAREC.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 10-jun-2004, VF
% 03-dec-2003, VF

% process input arguments
%----------------------------------
if nargin > 2,
  z=foo(varargin{1},varargin{2},varargin{3},varargin{4},varargin{5});
  return;
else
  model = varargin{1};
  if nargin < 2, options=[]; else options = c2s( varargin{2}); end
  if ~isfield(options,'min_improvement'), options.min_improvement = 1e-3; end
  if ~isfield(options,'start_t'), options.start_t = 1e-3; end
  if ~isfield(options,'attempts'), options.attempts = 10; end
end

[dim,num_sv]=size(model.sv.X);
ker = 'rbf';
arg = model.options.arg;
iXi = sum( model.sv.X.^2)';
s2 = arg^2;

% Selection of the starting point out of the model.sv.X.
% The point in which is the objective function minimal is taken.
% Minimum over 50 randomly drawn points is used.
%--------------------------------------------------------------

rand_inx = randperm( num_sv );
rand_inx = rand_inx(1:min([num_sv,50]));
Z = model.sv.X(:,rand_inx);

fval = kernel(Z,model.sv.X,ker,arg)*model.Alpha(:);
fval = -fval.^2;

[dummy, inx ] = min( fval );
z = Z(:,inx );

% Gradient descent optimization  
%--------------------------------------
change=inf;
opt=optimset('display','off','Diagnostics','off','LargeScale','off');
warning off MATLAB:divideByZero;
while change > options.min_improvement,
   
   % compute gradient
   dotp = kernel( model.sv.X,z,ker,arg ).*model.Alpha(:);
   dz = z*sum(dotp) - model.sv.X*dotp;

   % auxiciliary variables
   zXi = model.sv.X'*z;
   dzXi = model.sv.X'*dz;
   
   Ai = -(1/(2*s2)) * (iXi - 2*zXi + z'*z);
   Bi = -(1/s2) * (z'*dz - dzXi);
   C = -(1/(2*s2)) * (dz'*dz);
   
   % 1D-search to determine the size of the step in the gradient direction
   [t,fval] = fminunc('rbfpreimg2',options.start_t,opt,model.Alpha,Ai,Bi,C);
   
   old_z = z;
   z = z + dz*t;
      
   change = sum((z-old_z).^2);
end

return;

%---------------------------------
function f=foo(t,Alpha,Ai,Bi,C)
% 

f = Alpha(:)'*exp(Ai + Bi*t + C*t^2);
f = -f^2;

return;
% EOF


