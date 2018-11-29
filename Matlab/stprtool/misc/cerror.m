function error=cerror(y1,y2,label)
% CERROR Computes classification error.
%
% Synopsis:
%  error = cerror(ypred,ytrue)
%  error = cerror(ypred,ytrue,label) 
%
% Description:
%  error = cerror(ypred,ytrue) returns classification error, i.e.,
%      error=  find(ypred~=ytrue)/length(ytrue).
%
%  error = cerror(ypred,ytrue,label) considers only labels
%    find(ytrue==label), i.e., if ypred, ytrue are from {1,2} then
%
%   false_positives_rate = cerror(ypred,ytrue,2)
%   false_negatives_rate = cerror(ypred,ytrue,1)
%
% Input:
%  y1 [1 x n] Vector of integers (response of classifier).
%  y2 [1 x n] Vector of integers (ground truth).
%  label [int] Selected label.
%
% Output:
%  error [real] Error. 
%
% Example:
%  classifier  = [1,1,1,2]
%  groundtruth = [2,1,2,1]
%  error = cerror(classifier,groundtruth)
%  false_pos = cerror(classifier,groundtruth,2)
%  false_neg = cerror(classifier,groundtruth,1)
%
% See also 
%  ROC
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 08-aug-2005, VF
% 09-jun-2004, VF
% 14-Jan-2003, VF

y1=y1(:);y2=y2(:);

if nargin < 3,
  error=length(find((y1-y2)~=0))/length(y2);
else
  inx = find(y2==label);
  error = length( find(y1(inx)~=label) )/length(inx);
end

return;
