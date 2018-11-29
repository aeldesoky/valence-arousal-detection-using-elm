function varargout = pgmm( model, options )
% PGMM Vizualizes Gaussian mixture model.
%
% Synopsis:
%  pgmm( model );
%  pgmm( model, options );
%  h = pgmm( ... );
%
% Description:
%  It vizualizes univariate (dim=1) or bivariate (dim=2) Gaussin mixture 
%  model (GMM). In the univariate case it also displays mixture components. 
%  It returns handles of used graphics objects.
%
%  In the case of bivariate GMM trhere are two options of visualization:
%    countours of p.d.f.  ... options.visual = 'contour' (default)
%    surface of p.d.f.    ... options.visual = 'surf'
% 
% Input:
%  model.Mean [dim x ncomp] Mean values.
%  model.Cov [dim x ncomp] Covariances.
%  model.Prior [dim x ncomp] Mixture weights.
%
%  options.comp [1x1] If 1 (default) then it plots also mixture components.
%  options.visual [string] If equal to 'contour' then contour function is 
%    used if 'surf' then surf functions is used (see above).
%  options.adj_axes [1x1] If 1 (default) then axes are set to display 
%   whole mixture otherwise unchanged.
%  options.color [string] Color of GMM plot in univariate case (default 'b').
%
% Output:
%  h [1 x nobjects] Handles of used graphics object.
%
% Example:
%
% Univariate case:
%  model1 = c2s({'Mean',[-3 0 3],'Cov',[0.5 1 0.8],'Prior',[0.4 0.3 0.3]});
%  figure; pgmm(model1);
%
% Bivariate case:
%  model2.Mean(:,1) = [-1;-1]; model2.Cov(:,:,1) = [1,0.5;0.5,1];
%  model2.Mean(:,2) = [1;1]; model2.Cov(:,:,2) = [1,-0.5;-0.5,1];
%  model2.Prior = [0.4 0.6];
%  figure; pgmm(model2);
%  figure; pgmm(model2,{'visual','surf'});
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 2-may-2004, VF
% 29-apr2004, VF
% 8-mar-2004, VF

if nargin >= 2, options=c2s(options); else options = []; end
if ~isfield(options,'color'), options.color = 'b'; end
if ~isfield(options,'comp_color'), options.comp_color = options.color; end
if ~isfield(options,'comp'), options.comp = 1; end
if ~isfield(options,'adj_axes'), options.adj_axes = 1; end
if ~isfield(options,'visual'), options.visual = 'contour'; end

[dim,ncomp] = size( model.Mean);
a = axis;
old_hold = ishold;
hold on;

% univariate variances can be given as a vector
if size(model.Cov,1) ~= size(model.Cov,2), 
  model.Cov = reshape(model.Cov,1,1,ncomp); 
end

% univariate case
if dim == 1,
  
  a = axis;
  if options.adj_axes == 1,
    min_x = min(a(1),min(model.Mean)-sqrt(max(model.Cov))*3);
    max_x = max(a(2),max(model.Mean)+sqrt(max(model.Cov))*3);
  else
    min_x = a(1); 
    max_x = a(2);
  end
  
  x = linspace(min_x,max_x,200);

  h = plot(x, pdfgmm(x, model), options.color );

  if options.comp == 1,
    for i=1:ncomp,
      if length(options.comp_color) < i,
        comp_color = options.comp_color(1);
      else
        comp_color = options.comp_color(i);
      end
      h(i+1) = plot(x, ...
        model.Prior(i)*pdfgauss(x, model.Mean(:,i),model.Cov(:,:,i)), ...
        ['--' comp_color]);
    end
  end

% bivariate case
elseif dim ==2,

  a = axis;
  if options.adj_axes == 1,
    tmp=[];
    for i=1:ncomp,
      tmp(i) = max(eig(model.Cov(:,:,i)));
    end
    margin = sqrt(max(tmp))*3;
    min_x = min(a(1),min(model.Mean(1,:))-margin);
    max_x = max(a(2),max(model.Mean(1,:))+margin);
    min_y = min(a(3),min(model.Mean(2,:))-margin);
    max_y = max(a(4),max(model.Mean(2,:))+margin);
  else
    min_x = a(1);
    max_x = a(2);
    min_y = a(3);
    max_y = a(4);
  end
  
  [Ax,Ay] = meshgrid(linspace(min_x,max_x,50), linspace(min_y,max_y,50));
  y = pdfgmm([Ax(:)';Ay(:)'],model );
  
  switch options.visual
    case 'contour'
      h = contour( Ax, Ay, reshape(y,50,50)); 
    case 'surf'
      h = surf( Ax, Ay, reshape(y,50,50)); 
      shading interp;
 end
  
else
  error('GMM must be univariate or bivariate.');
end

% return handles if required
if nargout > 0, varargout{1} = h; end

if ~old_hold, hold off; end

return;
