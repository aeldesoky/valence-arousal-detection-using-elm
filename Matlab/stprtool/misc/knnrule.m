function model=knnrule(data,K)
% KNNRULE Creates K-nearest neighbours classifier.
%
% Synopsis:
%  model=knnrule(data)
%  model=knnrule(data,K)
%
% Description:
%  It creates model of the K-nearest neighbour classifier.
%
% Input:
%  data.X [dim x num_data] Prototypes (training) data.
%  data.y [1 x num_data] Labels of training data.
%  K [1x1] Number of the nerest neighbours (default 1).
%
% Output:
%  model [struct] Model of K-NN classifier.
%   .X = data.X.
%   .y = data.y.
%   .K = K.
%   .num_data [1x1] number of prototypes.
%   .fun [string] Contains 'knnclass'.
%
% Example:
%  data=load('riply_trn');
%  model=knnrule(data,1);
%  figure; ppatterns(data); pboundary(model);
%
% See also 
%  KNNCLASS.
%

data=c2s(data);

if nargin <2, K=1; end

model=data;
model.fun='knnclass';
model.K=K;
model.num_data = size(data.X,2);

return;


