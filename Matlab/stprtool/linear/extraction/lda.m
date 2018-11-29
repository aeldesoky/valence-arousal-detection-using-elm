function model = lda(data,new_dim)
% LDA Linear Discriminant Analysis.
% 
% Synopsis:
%  model = lda(data)
%  model = lda(data,new_dim)
%
% Description:
%  This function is implementation of Linear Discriminant Analysis.
%  The goal is to train the linear transform which maximizes ratio 
%  between between-class and within-class scatter matrix of projected 
%  data.
%
% Input:
%  data [struct] Input labeled data:
%   .X [dim x num_data] Data sample.
%   .y [1 x num_data] Labels (1,2,...,nclass).
%
%  new_dim [1x1] Output data dimension (default new_dim = dim).
%
% Ouput:
%  model [struct] Linear projection:
%   .W [dim x new_dim] Projection matrix.
%   .b [new_dim x 1] Biases.
%
%   .mean_X [dim x 1] Mean value of data.
%   .Sw [dim x dim] Within-class scatter matrix.
%   .Sb [dim x dim] Between-class scatter matrix.
%   .eigval [dim x 1] Eigenvalues.
%
% Example:
%  in_data = load('iris');
%  model = lda( in_data, 2 );
%  out_data = linproj( in_data, model);
%  figure; ppatterns(out_data);
%
% See also 
%  LINPROJ, PCA.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 05-oct-2006, VF
% 25-may-2004, VF
% 3-may-2004, VF
% 20-may-2001, V.Franc, created

% process input arguments
%----------------------------------------
data=c2s(data);
[dim,num_data] = size(data.X);
nclass = max( data.y );

if nargin < 2, new_dim = dim; end

% compute within-class scatter matrix
%--------------------------------------
mean_X = mean( data.X, 2);
Sw=zeros(dim,dim);
Sb=zeros(dim,dim);

for i = 1:nclass,
  inx_i = find( data.y==i);
  X_i = data.X(:,inx_i);
  
  mean_Xi = mean(X_i,2);
  Sw = Sw + cov( X_i', 1);
  Sb = Sb + length(inx_i)*(mean_Xi-mean_X)*(mean_Xi-mean_X)';
end

% Compute projection matrix
%[U,D,V2]=svd( inv( Sw )*Sb ); ERROR !
%model.W = V(:,1:new_dim);

[V,D]=eig( inv( Sw )*Sb );
[D,inx] = sort(diag(D),1,'descend');

% take new_dim biggest eigenvectors
model.W = V(:,inx(1:new_dim));
model.eigval = diag(D);

% translation
model.b = -model.W'*mean_X;

model.Sw = Sw;
model.Sb = Sb;
model.mean_X = mean_X;

model.fun = 'linproj';

return;
% EOF