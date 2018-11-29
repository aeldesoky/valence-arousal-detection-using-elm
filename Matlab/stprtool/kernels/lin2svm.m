function svm_model = lin2svm(kfe_model, lin_model)
% LIN2SVM Merges linear rule and kernel projection.
%
% Synopsis:
%  svm_model = lin2svm(kfe_model,lin_model)
%
% Description:
%  This function merges kernel feature extraction model
%  (data-type kernel projection) and linear classifier to 
%  create kernel (SVM) classifier.
%
% Input:
%  kfe_model [struct] Kernel data projection:
%   .Alpha [nsv x new_dim] Weight vector.
%   .b [new_dim x 1] Biases.
%   .oprions.ker [string] Kernel identifier (see 'help kernel').
%   .options.arg [1xnargs] Kernel arguments.
%
%  lin_model [struct] Linear classifier:
%   .W [dim x nfun] Weight vector(s).
%   .b [nfun x 1] Bias(es).
%  
% Output:
%  svm_model [struct] Kernel classifer:
%   .Alpha [nsv x nfun] Weight vector(s).
%   .b [nfun x 1] Bias(es).
%   .options [struct] Copy of kfe_model.options.
%
% Example:
%  data = load('riply_trn');
%  options = struct('ker','rbf','arg',1,'new_dim',10);
%  kpca_model = greedykpca(data.X,options);
%  proj_data = kernelproj(data,kpca_model);
%  lin_model = fld(proj_data);
%  kfd_model = lin2svm(kpca_model,lin_model);
%  figure; ppatterns(data); pboundary(kfd_model);
%
% See also 
%  LIN2QUAD, SVMCLASS, LINCLASS.
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 10-jun-2004, VF
% 02-Feb-2003, VF

svm_model.Alpha = kfe_model.Alpha*lin_model.W;
svm_model.b = lin_model.b+lin_model.W'*kfe_model.b;

svm_model.sv.X = kfe_model.sv.X;
svm_model.nsv = size(svm_model.sv.X,2);
svm_model.options = kfe_model.options;
svm_model.fun = 'svmclass';

return;
% EOF