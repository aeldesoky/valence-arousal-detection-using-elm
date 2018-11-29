function [h,m,s] = sectohms(sec)
% SECTOHMS Converts seconds to HOUR:MIN:SEC format.
%
% Synopsis:
%  [h,m,s] = sectohms(sec)
%

% Modifications:
% 2-jul-2007, VF, created

    
h = floor(sec/60/60);

sec = sec - h*60*60;

m = floor(sec/60);

sec = sec - m*60;

s = sec;

return;