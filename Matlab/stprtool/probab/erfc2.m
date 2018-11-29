function y=erfc2(x,m,v)
% ERFC2 Normal cumulative distribution function.
%
% Synopsis:
%  y=erfc2(x,m,v)
% 
% Description:
%  It computes normal cumulative distribution function 
%
%  erfc2(x) = 1/sqrt(2*pi)*v) * 
%    integral from x to inf of exp(-(t-m)^2/(2*v^2) dt
%
% Input:
%  x [1x1] Real number.
%  m [1x1] Mean value.
%  v [1x1] Variance.
%
% Output:
%  y [1x1] Error function.
%
% See also 
%  ERFC.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 28-apr-2004, VF

x2=(x-m)/(sqrt(2)*v);
y=erfc(x2)/2;

return;
