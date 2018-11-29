function m=marker_type(i)
% MARKER_TYPE Returns marker type.
%
% Synopsis:
%  m=marker_type(i)
%
% See also MARKER_COLOR. 
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
%  7-jan-2003, VF, created

MARKERS=['x','o','*','.','^','s','d','p','h','+'];

m=MARKERS(mod(i-1,length(MARKERS))+1);

return;
