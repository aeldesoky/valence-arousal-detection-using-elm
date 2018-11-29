function [FP,FN]=roc(dfce,y)
% ROC computes Receiver Operating Characteristic (ROC) curves. 
%
% Synopsis:
%  [FP,FN]=roc(dfce,y)
%  
% Description:
%  It computes false positive rate FP and false negative rate FN
%  with rescpect to the shift of the bias of given binary decision 
%  function. The values of the decision function are given in dfce 
%  and y contains true labels (number 1 and/or 2). The vectors dfce 
%  and y must be of the same length. 
%  The bias is shifted from min(dfce) to max(dfce). 
%
% Input:
%  dfce [1 x num_data] Values of decision function returned by 
%   a classifier.
%  y [1 x num_data] True labels.
%
% Output:
%  FP [1 x num_data] False positive rate.
%  FN [1 x num_data] False negative rate.
%
% Example:
%  data = load('riply_trn');
%  model = fld(data);
%  [y_pred,dfce] = linclass(data.X,model);
%  [FP,FN] = roc(dfce,data.y);
%  figure; hold on; plot(FP,FN);
%  xlabel('false positives'); 
%  ylabel('false negatives');
%
% See also 
%  CERROR
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 26-aug-2005, VF
% 17-may-2004, VF
% 6-June-2003, VF
% 24-Feb-2003, VF


num_data=length(dfce);
n1=length(find(y==1));
n2=length(find(y==2));

[dfce,inx]=sort(dfce);
y = y(inx);

FP=zeros(1,num_data);
FN=zeros(1,num_data);

wrong1=0;
wrong2=n2;

for i=1:num_data,
  if y(i) == 1,
    wrong1=wrong1+1;
  else
    wrong2=wrong2-1;
  end
  
  FP(i)=wrong2/n2;
  FN(i)=wrong1/n1;
end

return;
