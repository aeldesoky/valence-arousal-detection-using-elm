function K = kernel(varargin)
% KERNEL Evaluates kernel function.
%
% Synopsis:
%  K = kernel(X,ker,arg)
%  K = kernel(X1,X2,ker,arg)
%
% Description:
%  K = kernel( X, ker, arg ) returns kernel matrix K [n x n] 
%
%    K(i,j) = k(X(:,i),X(:,j))  for all i=1..n, j=1..n,
%
%   where k: a x b -> R is a kernel function given by 
%   identifier ker and argument arg:
%     
%   Identifier    Name           Definition
%   'linear'  ... linear kernel  k(a,b) = a'*b
%   'poly'    ... polynomial     k(a,b) = (a'*b+arg[2])^arg[1]
%   'rbf'     ... RBF (Gaussian) k(a,b) = exp(-0.5*||a-b||^2/arg[1]^2)
%   'sigmoid' ... Sigmoidal      k(a,b) = tanh(arg[1]*(a'*b)+arg[2])
%
%  K = kernel( X1, X2, ker, arg ) returns kernel matrix K [n1 x n2]
% 
%    K(i,j) = k(X1(:,i),X2(:,j))  for all i=1..n1, j=1..n2,
%
% Input:
%  X [dim x n] Single matrix of input vectors.
%  X1 [dim x n1], X2 [dim x n2] Pair of input matrices.
%  ker [string] Kernel identifier.
%  arg [1 x  narg] Kernel argument.
%
% Output:
%  K [n1 x n1] or K [n1 x n2] Kernel matrix.
% 
% Example:
%  X = rand(2,50);
%  K = kernel( X, 'rbf', 1);
%  figure; pcolor( K );
%
% See also:
%  DIAGKER, KERNELPROJ.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 19-sep-2004, VF
% 5-may-2004, VF

% MEX-File function.  

