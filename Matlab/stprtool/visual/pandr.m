function varargout = pandr(model,distrib)
% PANDR Visualizes solution of the Generalized Anderson's task.
%
% Synopsis:
%  h = pandr(model)
%
% Description:
%  It vizualizes solution of the Generalized Anderson's task 
%  for bivariate input Gaussians.
%  
%  The input of the task are two sets of Gaussians which 
%  describe the first and second class. The Gaussians are denoted as 
%  the ellipses (shape -> covariance, center -> mean). 
%  The output of the task is the linear classifier denoted as a line 
%  separating the 2D feature space.
%
% Input:
%  model [struct] Linear classifier:
%   .W [2 x 1] Normal vector of the separating hyperplane.
%   .b [real] Bias of the hyperplane.
%
%  distrib [struct] Set of binary labeled Gaussians:
%   .Mean [2 x ncomp] Mean vectors.
%   .Cov [2 x 2 x ncomp] Covariance matrices.
%   .y [1 x ncomp] Labels 1 or 2. 
%
% Output:
%  h [1 x nobjects] Handles of used graphics objects.
%  
% Example:
% 

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 4-may-2004, VF
% 24-feb-2003, VF
% 30-sep-2002, VF

[err,r,inx] = andrerr( model, distrib );

[dim, ncomp ] = size( distrib.Mean );
for i=1:ncomp,
  p(i) = exp(-0.5*r^2)/(2*pi*sqrt(det(distrib.Cov(:,:,i))));
end

h1 = pgauss( distrib, {'p',p});
h2 = pline( model );

if nargout > 0,
  varargout{1} = [h1 h2];
end

return;
