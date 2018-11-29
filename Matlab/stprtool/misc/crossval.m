function [itrn,itst]=crossval(num_data,num_folds)
% CROSSVAL Partitions data for cross-validation.
%
% Synopsis:
%  [itrn,itst] = crossval(num_data,num_folds)
%
% Description:
%  This function randomly partitions data into the training
%  and testing parts. The number of partitionings is determined
%  by the num_folds. If num_folds==1 then makes only one random
%  partitioning of data into training and testing in ratio 50:50.
% 
% Input:
%  num_data [1x1] number of data.
%  num_folds [1x1] number of folders.
%
% Output:
%  itrn{i} Indices of training data of i-th folder.
%  itst{i} Indices of testing data i-th folder.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 19-aug-2004, VF, comments repared
% 17-may-2004, VF
% 21-Jan-2003, VF
% 14-Jan-2003, VF
% 7-jan-2003, VF, created

% random partitioning 
inx=randperm(num_data);

itrn=cell(1,num_folds);
itst=cell(1,num_folds);

if num_folds == 1,
  half = fix(num_data/2);
  itrn{1}=inx(1:half);
  itst{1}=inx(half+1:end);
else

  num_column=ceil(num_data/num_folds);
  part=[1:num_data zeros(1,num_column*num_folds-num_data)];
  part=reshape(part,num_column,num_folds)';

  for i=1:num_folds,
    tst_inx=part(i,:);
    tst_inx=tst_inx(find(tst_inx));
    trn_inx=part(find([1:num_folds]~=i),:);
    trn_inx=trn_inx(find(trn_inx));

    itrn{i}=inx(trn_inx);
    itst{i}=inx(tst_inx);
  end
end

return;
% EOF
