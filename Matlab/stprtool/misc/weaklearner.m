function model = weaklearner(data)
% WEAKLEARNER Produce classifier thresholding single feature.
%
% Synopsis:
%  model = weaklearner(data)
%
% Description:
%  This function produce a weak binary classifier which assigns
%  input vector x to classes [1,2] based on thresholding a single 
%  feature. The output is a model which defines the threshold 
%  and feature index such that the weighted error is minimized.
%  This weak learner can be used with the AdaBoost classifier
%  (see 'help adaboost') as a feature selection method.
%  
% Input:
%  data [struct] Training data:
%   .X [dim x num_data] Training vectors.
%   .y [1 x num_data] Binary labels (1 or 2).
%   .D [1 x num_data] Weights of training vectors (optional).
%    If not given then D is set to be uniform distribution.
% 
% Output:
%  model [struct] Binary linear classifier:
%   .W [dim x 1] Normal vector of hyperplane.
%   .b [1x1] Bias of the hyperplane.
%   .fun = 'linclass'.
%
% Example:
%  help adaboost
%
% See also: 
%  ADABOOST, ADACLASS.
% 

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2004, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 25-aug-2004, VF
% 11-aug-2004, VF

[dim,num_data] = size(data.X);

W = zeros(dim,1);
Errors = zeros(dim,1);

for i=1:dim,

  x = [-inf,sort(unique(data.X(i,:)))];
 
  err = []; f = [];
  for j=1:length(x)-1,
    
    f(j) = 0.5*(x(j)+x(j+1));

    y = ones(1,num_data);
    y(find(data.X(i,:)< f(j))) = 2;

    err(j) = sum((y(:)~=data.y(:)).*data.D(:));
    
  end

  [minerr1,inx1] = min(err);
  [minerr2,inx2] = min(1-err);
  if minerr1 < minerr2,
    W(i) = 1;
    Errors(i) = minerr1;
    b(i) = -f(inx1);
  else
    W(i) = - 1;
    Errors(i) = minerr2;
    b(i) = f(inx2);
  end
  
end

[dummy,inx] = min(Errors);

model.W = zeros(dim,1);
model.W(inx) = W(inx);
model.b = b(inx);
model.fun = 'linclass';

y = linclass(data.X,model);
err = sum((y(:)~=data.y(:)).*data.D(:));

return;

%EOF