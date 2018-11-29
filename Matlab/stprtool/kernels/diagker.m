function diagK = diagker( X, ker, arg )
% DIAGKER Returns diagonal of kernel matrix of given data.
%
% Synopsis:
%  diagK = diagker( X, ker, arg )
%
% Description:
%  It returns the same vector as command
%    diagK = diag( kernel(X,ker,arg) )
% 
%  but it is more efficiently computed. See 'help kernel' 
%  to more info about implemented kernels and their arguments.
%
% Input:
%  X [dim x num_data] Input data.
%  ker [string] Kernel identifier.
%  arg [1 x nargs] Kernel argument.
%
% Output:
%  diagK [num_data x 1] Diagonal of kernel matrix.
%
% See also 
%  KERNEL.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 5-may-2004, VF

% MEX-File function. 
