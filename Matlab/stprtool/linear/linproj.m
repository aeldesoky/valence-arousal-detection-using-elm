function out1=linproj(arg1, model)
% LINPROJ Linear data projection.
%
% Synopsis:
%  Y = linproj(X, model)
%  out_data = linproj(in_data, model)
%
% Description:
%  Y = linproj(X, model) linearly projects data in X such that
%    Y = model.W'*X + model.b
%
%  out_data = linproj(in_data, model) projects in_data.X 
%    out_data.X = model.W'*in_data.X + model.b
%    out_data.y = in_data.y
%
% Input:
%  model [struct] linear projection:
%   .W [dim x ncomp] Projection matrix.
%   .b [ncomp x 1] Bias.
% 
% Example:
%  help pca;
%  help lda;
%
% See also 
%  PCA, LDA, KERNELPROJ.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 3-may-2004, VF
% 21-jan-03, VF
% 16-Jun-2002, VF

if isstruct(arg1),
  [dim,num_data]=size(arg1.X);

  out1.X = model.W'*arg1.X + model.b(:)*ones(1,num_data);
  out1.y = arg1.y;
  
else
  [dim,num_data]=size(arg1);

  out1 = model.W'*arg1 + model.b(:)*ones(1,num_data);
end

return;
