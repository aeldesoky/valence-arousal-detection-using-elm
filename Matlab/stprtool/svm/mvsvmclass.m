function [y,votes] = mvsvmclass(X,model)
% MVSVMCLASS Majority voting multi-class SVM classifier.
%
% Synopsis:
%  [y,votes] = mvsvmclass(X,model)
%
% Description:
%  [y,votes] = mvsvmclass(X,model) multi-class SVM classifier 
%    based on majority voting. The classifier involves nrule
%    binary rules each classifying into one of nclass labels.
%    The final decision is make for the class with majority 
%    votes.
%
% Input:
%  X [dim x num_data] Input vectors to be classified.
%
%  model [struct] Multi-class SVM majority voting classifier:
%   .Alpha [nsv x nrule] Weights.
%   .bin_y [2 x nrule] Translation between binary responses of
%     the discriminant functions and class labels.
%   .b [nrule x 1] Biases of discriminant functions.
%   .sv.X [dim x nsv] Support vectors.
%   .options.ker [string] Kernel identifier; see 'help kernel'.
%   .options.arg [1 x nargs] Kernel agrument(s).
%
% Output:
%  y [1 x num_data] Predicted labels.
%  votes [nclass x num_data] Number of votes for each class.
%
% Example:
%
% See also 
%  OAOSVM, SVMCLASS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications
%  11-Feb-2003, VF 
%  8-Feb-2003, VF 
%  3-Jun-2002, V.Franc

[dim,num_data] = size(X);
nclass = max( model.bin_y(:) );
nrule = size( model.Alpha, 2);

votes = zeros(nclass, num_data );

dfce = kernelproj( X, model );

for i=1:nrule,
  
  inx_pos = find( dfce(i,:) >= 0 );
  inx_neg = find( dfce(i,:) < 0 );

  votes( model.bin_y(1,i), inx_pos) = votes( model.bin_y(1,i), inx_pos) + 1;
  votes( model.bin_y(2,i), inx_neg) = votes( model.bin_y(2,i), inx_neg) + 1;

end

[dummy, y] = max( votes );

return;
% EOF
