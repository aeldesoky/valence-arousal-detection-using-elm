function y = knnclass(X,model)
% KNNCLASS k-Nearest Neighbours classifier.
%
% Synopsis:
%  y = knnclass(X,model)
%
% Description:
%  The input feature vectors X are classified using the K-NN
%  rule defined by the input model.
% 
% Input:
%  X [dim x num_data] Data to be classified.
%  model [struct] Model of K-NN classfier:
%   .X [dim x num_prototypes] Prototypes.
%   .y [1 x num_prototypes] Labels of prototypes.
%   .K [1x1] Number of used nearest-neighbours.
%
% Output:
%  y [1 x num_data] Classified labels of testing data.
%
% Example:
%  trn = load('riply_trn');
%  tst = load('riply_tst');
%  ypred = knnclass(tst.X,knnrule(trn,5));
%  cerror( ypred, tst.y )
%
% See also 
%  KNNRULE.
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 19-may-2003, VF
% 18-sep-2002, V.Franc

X=c2s(X);
model=c2s(model);

if ~isfield(model,'K'), model.K=1; end;

y = knnclass_mex(X,model.X,model.y, model.K);

return; 
% EOF

