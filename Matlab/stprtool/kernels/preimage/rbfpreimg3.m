function x = rbfpreimg3(model,nn)
% RBFPREIMG3 RBF pre-image problem by Kwok-Tsang's algorithm.
%
% Synopsis:
%  x = rbfpreimg3(model)
%  x = rbfpreimg3(model,nn)
%
% Description:
%  x = rbfpreimg3(model) is an implementation of the algorithm 
%    by [Kwok03] to solve the pre-image problem for kernel expansion 
%    with RBF kernel. The kernel expansion is given in the input 
%    structure model.
%
%  x = rbfpreimg3(model,nn) use to set number of nearest 
%    neighbours (default 10).
%
% Input:
%  model [struct] RBF kernel expansion:
%   .Alpha [nsv x 1] Weights.
%   .sv.X [dim x nsv] Vectors.
%   .options.arg [1x1] RBF kernel width.
%
%  nn [1x1] Number of nearest neighbours.
%  
% Output:
%  x [dim x 1] Pre-image of the RBF kernel expansion.
%  
% See also 
%  RBFPREIMG1, RBFPREIMG2, RSRBF, KPCAREC.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 17-may-2004, VF
% 21-Feb-2004, VF
% 17-February-2004, Petr Posik

[dim, num_data] = size(model.sv.X);

% default number of used nearest neighbours 
if nargin < 2, nn = min([num_data, 10]); end

K = kernel(model.sv.X, 'rbf', model.options.arg );

Const2 = model.Alpha(:)'*K*model.Alpha(:);
df2 = 1 + Const2 - 2*K*model.Alpha(:);
d2 = -2*model.options.arg^2 * log( 1 - 0.5*df2);

% select nn neighbours
[dummy, inx] = sort( df2 );
X = model.sv.X(:,inx(1:nn));
df2 = df2(inx(1:nn));
d2 = d2(inx(1:nn));

H = eye(nn,nn) - 1/nn * ones(nn,nn);

[U,L,V] = svd(X*H);
r = rank(L);

Z = L*V';

d02 = sum(Z.^2)';

z = -0.5*pinv(Z')*(d2-d02);
x = U*z + sum(X,2)/nn;

return;

