function varargout=plane3(arg1,arg2)
% PLANE3 Plots plane in 3d.
%
% Synopsis:
%  h=plane3( W, b )
%  h=plane3( model )
%
% Description:
%  h=plane3( W, b ) plots the plane described as
%     W'*x + b = 0.
%
%  h=plane3( model ) structure model contains plane 
%   parametres W and b.
%
% Output:
%  h [1 x nobjects] Handles of used graphical objects.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 17-may-2004, VF
% 13-apr-2003, VF
% 11-aug-2003, VF

if isstruct(arg1),
  W = arg1.W;
  b = arg1.b;
  if nargin >= 2, col = arg2; else col = 'g'; end
else
  W = arg1;
  b = arg2;
  if nargin >= 3, col = arg3; else col = 'g'; end
end

ax=axis;

A=[ax(1) ax(3)];
B=[ax(2) ax(3)];
C=[ax(2) ax(4)];
D=[ax(1) ax(4)];

W=W(:);

A=[A, -(b + A*W(1:2))/W(3)];
B=[B, -(b + B*W(1:2))/W(3)];
C=[C, -(b + C*W(1:2))/W(3)];
D=[D, -(b + D*W(1:2))/W(3)];
  
X=[A(1) B(1) C(1) D(1)];
Y=[A(2) B(2) C(2) D(2)];
Z=[A(3) B(3) C(3) D(3)];

h=patch(X,Y,Z,col);
set(h,'faceAlpha',0.5);

if nargout >= 1, varargout{1} = h; end

return;
% EOF