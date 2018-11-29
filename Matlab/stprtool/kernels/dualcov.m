function Z=dualcov(varargin)
% DUALCOV Dual representation of covariance matrix.
%
% Synopsis:
%  Z=dualcov(num_data)
%  Z=dualcov(labels, y)
%
% Description:
%  This function computes a matrix Z [num_data x num_data] which allows 
%  to express the sample covariance matrix of data sample X [dim x num_data] 
%  in terms of dot products. 
%
%  Z = dualcov(num_data) computes a matrix Z [num_data x num_data] such that 
%    cov(X',1) = X*Z*X'.
%
%  m = dualcov(labels,y) computes a matrix Z [length(y) x length(y)] such that
%    cov(X(:,find(labels==y))',1) = X*Z*X',
%
%  where labels [1 x num_data] is a vector of data labels and  y [1x1]  
%  is a label od class which covariance metrix is to be computed.
%
% Example:
%  Unlabeled data:
%   data = load('riply_trn');
%   ca = cov( data.X', 1)
%   cb = data.X*dualcov(size(data.X,2))*data.X'
%
%  Labeled data:
%   data = load('riply_trn');
%   ca1 = cov( data.X(:,find(data.y==1))',1)
%   cb1 = data.X*dualcov(data.y,1)*data.X'
%   ca2 = cov( data.X(:,find(data.y==2))',1)
%   cb2 = data.X*dualcov(data.y,2)*data.X'
%
% See also 
%  DUALMEAN.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 16-may-2004, VF
% 14-may-2004, VF
% 22-Jan-2003, VF
% 22-May-2001, V. Franc, created

if nargin == 2,
  % labeled data
  
  labels = varargin{1};
  y = varargin{2};
 
  num_data = length(labels);  
  inx_y=find( labels == y);
  n = length(inx_y);
  
  J=ones(n,1)/n; 
  Zy = ( eye(n,n)/n - J*J' ); 

  Z = zeros(num_data,num_data);
  Z(inx_y,inx_y) = Zy;
else
  % unlabeled data
  
  n = varargin{1};
  J=ones(n,1)/n; 
  Z = ( eye(n,n)/n - J*J' ); 
end

return;
% EOF