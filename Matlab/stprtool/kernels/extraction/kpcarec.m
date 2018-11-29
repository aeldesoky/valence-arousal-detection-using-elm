function Y = kpcarec(X,model)
% KPCAREC Reconstructs image after kernel PCA.
% 
% Synopsis:
%  Y = kpcarec(X,model)
%
% Description:
%  Input data X are projected using kernel projection trained
%  the by Kernel PCA [Mika99b]. The RBF kernel is assumed. This 
%  function computes the preimages Y from the input space 
%  corresponding to the projected data are.
%
%   X -> projection to -> preimage -> Y
%        kernel space     problem
%        by Kernel PCA
%
% Input:
%  X [dim x num_data] Input vectors.
%  model [struct] Kernel projection with RBF kernel;
%   see 'help kernelproj'. 
%
% Output:
%  Y [dim x num_data] Output data.
%
% See also  
%  KPCA, PCAREC.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 17-may-2004, VF
% 22-apr-2004, VF
% 17-mar-2004, VF, created.

[dim, num_data] = size(X);

fprintf('Projection data...');
Z = kernelproj(X, model );
fprintf('done.\n');

% allocate memory
Y = zeros(dim,num_data);
img = model;

fprintf('Computing preimages');
for i=1:num_data,
   fprintf('.');
  
   img.Alpha = model.Alpha*(Z(:,i) - model.b);

   Y(:,i) = rbfpreimg(img);       % Schoelkopf's fix-point algorithm
%   Y(:,i) = rbfpreimg2(img);  % Gradient method
%    Y(:,i) = rbfpreimg3(img,7);     % Kwok & Tsang

end
fprintf('done\n');

return;
% EOF