function y=sigmoid(X,arg1,arg2)
% SIGMOID Evaluates sigmoid function.
%
% Synopsis:
%  y = sigmoid(X,model)
%  y = sigmois(X,A1,A2)
%
% Description:
%  y = sigmoid(X,model) returns 
%     y = 1/(1+exp(A(1)*X + A(2))
%
%  where A = model.A.
%
%  y = sigmois(X, A1, A2) allows A to be given as A = [A1 A2].
%
% Input:
%  X [1 x num_data] Inputs.
%  model.A [2 x 1] Sigmoid parameters.
%
% Output:
%  y [1 x num_data] Evaluated sigmoid function.
%
% See also 
%  MLSIGMOID.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 3-jun-2004, VF
% 8-may-2003, VF

if nargin == 3,
  y = 1./(1+exp(X*arg1 + arg2));
else
  model=c2s(arg1);
  y = 1./(1+exp(X*model.A(1)+model.A(2)));
end

return;
% EOF