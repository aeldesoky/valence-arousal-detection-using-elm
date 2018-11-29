function varargout=ellips(Center,Shape,radius,n)
% ELLIPS Creates ellipse.
%
% Synopsis:
%  [X,Y] = ellips(Center,Shape,radius,n)
%  [X,Y,Z] = ellips(Center,Shape,radius,n)
%
% Description:
%  This function interpolates ellipse by lines in 2d space
%  or by patches in 3d space respectivelly. The ellipsoid 
%  is described as
%         radius^2=(x-Center)'*Shape*(x-Center).
%   
%  The number of lines used for interpolation is given
%  by argument n in 2d case. In 3d case the argument
%  n has the same meaning in the Matlab function sphere(n).
%
% Input:
%  Center [2x1] or [3x1] Center of the ellipse.
%  Shape [2x2] or [3x3] Shape of the ellipse.
%  n [1x1] Density of interpolation (default 20).
%
% Example:
%
% 2d ellipse
%  [x,y] = ellips([1;1],[1 0.5;0.5 1],1);
%  figure; plot(x,y);
%
% 3d ellipsoid
%  [x,y,z] = ellips([1;1;1],[1 0 0;0 2 0; 0 0 3],1);
%  figure; mesh(x,y,z);
%

% About: Statistical Pattern Recognition Toolbox                                
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac                     
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>            
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>       
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>           

% Modifications:
% 30-apr-2004, VF
% 13-Feb-2003, VF
% 24-6-2000 V. Hlavac, comments polished.

Center=Center(:);

% input arguments processing
if nargin < 4, n=20; end

% dimension
dim = size(Shape,1);

if dim == 3,
   % 3D ellipsoid, creates the matrix X,Y,Z appropriate for using in
   % the functions mesh,surf, etc.

   [X,Y,Z] = sphere(n);

   X=radius*X;
   Y=radius*Y;
   Z=radius*Z;

   [A,p]=chol(Shape);
   if p ~= 0,
      error('Shape matrix must be positive definite');
   end
   A=inv(A);

   [ROWS,COLUMNS]=size(X);

   % transforms sphere to ellipse 
   for i=1:ROWS,
      P=[X(i,:);Y(i,:);Z(i,:)];
      Q=A*P;

      % if the translation is given then translate points
      Q=Q+repmat(Center(:),1,COLUMNS);

      X(i,:)=Q(1,:);
      Y(i,:)=Q(2,:);
      Z(i,:)=Q(3,:);
   end % for i=1:ROWS

   % return variables
   varargout{1}=X;
   varargout{2}=Y;
   varargout{3}=Z;

elseif dim == 2,

   from=0;
   to=2*pi;
   step=(to-from)/n;
   X=cos([from:step:to]);
   Y=sin([from:step:to]);

   X=radius*X;
   Y=radius*Y;

   [A,p]=chol(Shape);
   if p ~= 0,
      error('Shape matrix must be positive definite');
   end
   A=inv(A);

   P=[X;Y];
   Q=A*P;

   Q=Q+repmat(Center(:),1,n+1);

   X=Q(1,:);
   Y=Q(2,:);

   % return variables
   varargout{1}=X;
   varargout{2}=Y;

end
 
return;
%EOF
