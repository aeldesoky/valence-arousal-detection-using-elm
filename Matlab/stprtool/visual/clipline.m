function [x1,y1,x2,y2,inside] = clipline(W,b,window)
% CLIPLINE clips line into given window.
%
% Synopsis:
%  [x1,y1,x2,y2,inside] = clipline(W,b,window)
%
% Description:
%  This function returns 2d points (x1,y1) and (x2,y2) 
%  of a line segment given by line
%    W'*x + b = 0
%
%  clipped to the given window.
%
% Input:
%  W [2x1] Normal of line.
%  b [1x1] Line threshold.
%  window [4x1] Contains [left right top bottom].
%
% Output:
%  x1 [1x1], y1 [1x1] The first point of line segment.
%  x2 [1x1], y2 [1x1] The second point of line segment.
%  inside [1x1] 1 if line W'*x+b=0 intersects the window. 
%  
% Example:
%  figure; hold on;
%  axis([-1 1 -1 1]);
%  window = 0.5*axis;
%  [x1,y1,x2,y2]=clipline([-1;1],0,window);
%  plot([x1 x2],[y1 y2]);
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 30-apr-2004, VF


theta = -b;
alpha = W;
minx=window(1);
maxx=window(2);
miny=window(3);
maxy=window(4);

x=zeros(4,1);
y=zeros(4,1);

if alpha(1)==0,
   if alpha(2)~=0,
      x1=minx;
      y1=theta/alpha(2);
      x2=maxx;
      y2=y1;
      inside=1;
   else
      % if alpha == 0 then it means the bad input.
      x1=0;
      y1=0;
      x2=0;
      y2=0;
      inside=0;
   end
elseif alpha(2)==0,
   x1=theta/alpha(1);
   y1=miny;
   x2=x1;
   y2=maxy;
   inside=1;
else
   y(1)=maxy;
   x(1)=(theta-alpha(2)*y(1))/alpha(1);
   y(2)=miny;
   x(2)=(theta-alpha(2)*y(2))/alpha(1);

   x(3)=maxx;
   y(3)=(theta-alpha(1)*x(3))/alpha(2);
   x(4)=minx;
   y(4)=(theta-alpha(1)*x(4))/alpha(2);

   j=0;
   for i=1:4,
      if x(i) <= maxx & x(i) >= minx & y(i) <= maxy & y(i) >= miny,
         if j==0,
            j=j+1;
            x1=x(i);
            y1=y(i);
         elseif j==1,
            j=j+1;
            x2=x(i);
            y2=y(i);
         end
      end
   end

   if j<2,
      x1=0;
      y1=0;
      x2=0;
      y2=0;
      inside=0;
   else
      inside=1;
   end
end % elseif alpha(2)==0
 
return; 
