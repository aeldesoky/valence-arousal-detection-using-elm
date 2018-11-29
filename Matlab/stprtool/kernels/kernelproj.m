function out_data=kernelproj(in_data, model)
% KERNELPROJ Kernel projection.
%
% Synopsis:
%  Y = kernelproj(X, model)
%  out_data = kernelproj(in_data, model)
%
% Description:
%  Y = kernelproj(X, model) this function maps input vectors 
%    X [dim x num_data] onto vectors Y [new_dim x num_data]
%    using the kernel projection
%
%      Y(:,i) = Alpha' * kernel(sv.X, X(:,i), ker, arg) + b
%
%    where parameters of the projection are given in model:
%     .Alpha [nsv x new_dim] Multipliers.
%     .b [new_dim x 1] Bias.
%     .sv.X [dim x nsv] Vectors.
%     .options.ker [string] Kernel identifier.
%     .options.arg [1 x narg] Kernel argument.
%
%  out_data = kernelproj(in_data, model) assumes that in_data
%   is a structure containing vectors X and labels y. 
%   The output structute out_data is constructed as 
%
%     out_data.X = kernelproj(in_data.X, model)
%     out_data.y = in_data.y
%
% Example:
%  help kpca;
%  help gda;
%
% See also 
%  GDA, KPCA, LINPROJ, KERNEL.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 19-sep-2004, VF, core of the function rewritten to C
% 14-may-2004, VF
% 4-may-2004, VF

if isstruct(in_data)==1,
  out_data = in_data;

  if ~isempty(model.Alpha) & isfield(model, 'Alpha'),
    out_data.X = kernelproj_mex(in_data.X, model.Alpha, model.b, ...
        model.sv.X, model.options.ker, model.options.arg);
  else
    [dim,num_data]=size(in_data.X);
    out_data.X = model.b*ones(1,num_data);
  end
  
else
  if ~isempty(model.Alpha) & isfield(model,'Alpha'),
    
    out_data = kernelproj_mex(in_data, model.Alpha, model.b, ...
        model.sv.X, model.options.ker, model.options.arg);
  else
    [dim,num_data]=size(in_data);
    out_data = model.b*ones(1,num_data);
  end
end

return;
% EOF