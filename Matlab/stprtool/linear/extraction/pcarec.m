function Y = pcarec(X,model)
% PCAREC Computes reconstructed vector after PCA projection.
% 
% Synopsis:
%  Y = pcarec(X,model)
%
% Description:
%  The input vectorts X are projected onto Z using linear 
%  projection trained by the Principal Component Analysis (PCA). 
%  The vectors Y are computed from Z as a reconstruction of 
%  the original vectors X:
%
%       PCA     Reconstr
%    X  --->  Z   --->   Y
%
% Input:
%  X [dim x num_data] Input vectors.
%  model [struct] Linear projection trained by PCA. 
%
% Output:
%  Y [dim x num_data] Reconstructed vectors.
%
% See also 
%  LINPROJ, PCA, KPCAREC.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 25-may-2004, VF
% 5-may-2004, VF
% 22-apr-2004, VF
% 17-mar-2004, VF, created.

[dim,num_data] = size(X);

Y = model.W*linproj(X,model) + model.mean_X*ones(1,num_data);

return;  
% EOF  
