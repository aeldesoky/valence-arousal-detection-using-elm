function X = dsamp(Px,num_data)
% DSAMP Generates samples from discrete distribution.
%
% Synopsis:
%  X = dsamp(Px,num_data)
%
% Input:
%  Px [dim x 1] Discrete probability distribution; it must
%   satisfy sum(Px) = 1 and min(Px) >= 0.
%  num_data [1 x 1] Number of samples to be generated.
% 
% Example:
%  Px = [0.2 0.3 0.1 0.4];
%  X = dsamp(Px,1000);
%  Px_emp = hist(X,4)/1000
%
% See also
%  GSAMP, GMMSAMP

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2005, Written by Vojtech Franc and Vaclav Hlavac
% Czech Technical University Prague
% Faculty of Electrical Engineering
% Center for Machine Perception

% Modifications:
% 26-mar-2005, VF

unif_X = rand(1,num_data);
X = zeros(1,num_data);

cum_Px = [0; tril(ones(length(Px)))*Px(:)];
cum_Px(end) = inf;
for i=1:length(Px),
  inx = find( unif_X >= cum_Px(i) & unif_X < cum_Px(i+1) );
  X(inx) = i;
end

return;
% EOF