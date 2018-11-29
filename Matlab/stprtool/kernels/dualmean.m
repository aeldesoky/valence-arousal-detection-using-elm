function m=dualmean(varargin)
% DUALMEAN Computes dual representation of mean vector.
%
% Synopsis:
%  m = dualmean(num_data)
%  m = dualmean(labels,y)
%
% Description:
%  This function computes a vector m which allows to express the mean
%  vector of data sample X [dim x num_data] in terms of dot products. 
%
%  m = dualmean(num_data) computes a vector m [num_data x 1] such that 
%    mean(X,2) = X*m.
%
%  m = dualmean(labels,y) computes a vector m [length(y) x 1] such that
%    mean(X(:,find(labels==y)),2) = X*m,
%
%  where labels [1 x num_data] is a vector of data labels and  y [1x1]  
%  is a label od class which mean vector is to be computed.
% 
% Example:
%  Unlabeled data:
%   data = load('riply_trn');
%   ma = mean( data.X, 2)
%   mb = data.X*dualmean(size(data.X,2))
%
%  Labeled data:
%   data = load('riply_trn');
%   ma1 = mean( data.X(:,find(data.y==1)),2)
%   mb1 = data.X*dualmean(data.y,1)
%   ma2 = mean( data.X(:,find(data.y==2)),2)
%   mb2 = data.X*dualmean(data.y,2)
%
% See also 
%  DUALCOV.
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
  m = zeros(num_data,1);
  inx = find(labels==y);
  m(inx) = ones(length(inx),1)/length(inx);
else
  % unlabeled data
  num_data = varargin{1};
  m = ones(num_data,1)/num_data;
end

return;
% EOF