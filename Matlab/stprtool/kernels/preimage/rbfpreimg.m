function x = rbfpreimg( model, options, init_point  )
% RBFPREIMG RBF pre-image by Schoelkopf's fixed-point algorithm.
%
% Synopsis:
%  x = rbfpreimg( model )
%  x = rbfpreimg( model, options  )
%  x = rbfpreimg( model, options, init_point  )
%
% Description:
%  x = rbfpreimg( model ) it is an implementation of the 
%   Schoelkopf's fixed-point algorithm to solve the pre-image 
%   problem for kernel expansion wiht RBF kernel [Schol98a]. 
%   The kernel expansion is given in the input structure model.
%
%  x = rbfpreimg( model, options  ) use structure options to 
%   set up control parameters: 
%    tmax ... maximal number of iterations.
%    eps ... minimal change in the norm of the optimized vector. 
%
%  x = rbfpreimg( model, options, init_point  ) use to set up
%    starting point of the optimization otherwise it is seleceted 
%    randomly.
%
% Input:
%  model [struct] Kernel expansion:
%   .Alpha [nsv x 1] Coefficients of the kernel expansion.
%   .sv.X [dim x nsv] Vectors of defining the expansion.
%   .options.ker [string] Must be equl to 'rbf'.
%   .options.arg [1x1] Argument of the RBF kernel.
% 
%  options [struct] Control parameters:
%   .tmax [1x1] Maximal number of iterations (default 1e6).
%   .eps [1x1] Minimal change of the optimized vector x.
%
%  init_point [dim x 1] Initial point of optimization.
%  
% Output:
%  x [dim x 1] Pre-image of the RBF kernel expansion.
%  
% See also 
%  RBFPREIMG2, RBFPREIMG3, RSRBF, KPCAREC.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 4-may-2004, VF
% 1-July-2003, VF
% 30-June-2003, VF

% input arguments
%-------------------------------
if nargin < 2, options = []; else options = c2s( options); end
if ~isfield(options,'tmax'), options.tmax = 1e6; end
if ~isfield(options,'eps'), options.eps = 1e-12; end
[dim,nsv] = size(model.sv.X);

% select initial point 
%-----------------------
if nargin == 3,
  x = init_point;
else
  x = model.sv.X*(2*rand(nsv,1)-ones(nsv,1));
end
old_x = x;

step=0;
change = inf;

% main loop
%----------------------
while step < options.tmax & change > options.eps,
 
   step = step + 1;
   
   k = kernel(model.sv.X, x, model.options.ker, model.options.arg );

   ka = k.*model.Alpha(:);

   den=sum(ka);
   if den > 0,
     x = model.sv.X*ka/den;
   else
     % reinit
     x = model.sv.X*(2*rand(nsv,1)-ones(nsv,1));
   end

   change = sum((old_x-x).^2);
   
   old_x = x;
end
 
return;
