function dist = mahalan(X,Mean,Cov)
% MAHALAN Computes Mahalanobis distance.
%
% Synopsis:
%  dist = mahalan(X,Mean,Cov)
%
% Description:
%  It computes Mahalanobis distance between column vectors 
%  of matrix X and vector Mean with matrix Cov, i.e.,
%
%    dist(i) = (X(:,i)-Mean)'*inv(C)*(X(:,i)-Mean)
%
%  for all i=1:size(X,2).
%
% Input:
%  X [dim x num_data] Input data.
%  Mean [dim x 1] Vector.
%  Cov [dim x dim] Matrix.
%
% Output:
%  dist [1 x num_data] Mahalanobis distance.
%
% Example:
%  It plots isolines of Mahalanobis distance.
%
%  [Ax,Ay] = meshgrid(linspace(-5,5,100), linspace(-5,5,100));
%  dist = mahalan([Ax(:)';Ay(:)'],[0;0],[1 0.5; 0.5 1]);
%  figure; contour( Ax, Ay, reshape(dist,100,100));
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 28-apr-2004, VF

[dim, num_data] = size( X );

XC = X - repmat(Mean,1,num_data);
dist= sum((XC'*inv( Cov ).*XC')',1);

return;
