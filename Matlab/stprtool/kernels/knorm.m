function kn = knorm(X,Alpha,ker,arg)
% KNORM Computes L2-norm in kernel space.
%
% Synopsis:
%  kn = knorm(X,Alpha,ker,arg)
%
% Description:
%  kn = knorm(X,Alpha,ker,arg) computes kn = sqrt(Alpha'*K*Alpha)
%   where K = kernel(X,ker,arg) is the kernel matrix.
%
% Input:
%   Alpha [num_data x 1] Expansion coefficients (weights).
%   X [dim x num_data] Data vectors.
%   ker [string] Kernel identifier (see "help kernel").
%   arg [...] Kernel argument.
%
% Ouput:
%   kn [1x1] Norm in the kernel space.
% 
% See also: 
%  KERNEL, KERNELPROJ.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2004, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 08-oct-2004, created

kn = sqrt( kernelproj_mex(X, Alpha(:), 0, X, ker, arg)*Alpha(:) );

return;

% EOF
