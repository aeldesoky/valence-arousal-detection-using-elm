function c=marker_color(i)
% MARKER_COLOR Returns marker color.
%
% Synopsis:
%  c=marker_color(i)
%
% See also MARKER_TYPE
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
%  7-jan-2003, VF, created

%COLORS=['r','b','y','m','c','k','g'];
COLORS=['b','r','g','k','m','c','y'];

c=COLORS(mod(i-1,length(COLORS))+1);

return;

