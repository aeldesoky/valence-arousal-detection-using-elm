function varargout=pgauss(model,options)
% PGAUSS Vizualizes set of bivariate Gaussians.
%
% Synopsis:
%  pgauss(model)
%  pgauss(model,options)
%  h = pgauss(...)
%
% Description:
%  pgauss(model) visualizes a set of bivariate Gaussians as
%   isolines (ellipse) with equal probability density functions.
%   The Gaussians are given by mean vectors model.Mean [2xncomp]
%   and covariance matrices model.Cov [2x2xncomp]. If labels
%   model.y [1xncomp] are given then the Gaussians are distinguished
%   by colors correspoding to labels.
% 
%  pgauss(model,options) structure options controls visualization;
%   If options.fill equals 1 then Ellipses are filled otherwise only
%   contours are plotted. The isolines to be drawn are given by 
%   values of probability distribution function in field 
%   options.p [1xncomp]. If length(option.p)==1 then isolines for
%   all Gaussians are drawn for the same value.
%  
%  h = pgauss(...) returns handles of used graphics objects.
%    
% Input:
%  model [struct] Parameters of Gaussian distributions:
%   .Mean [2 x ncomp] Mean vectors of ncomp Gaussians.
%   .Cov [2 x 2 x ncomp] Covariance matrices.
%   .y [1 x ncomp] (optional) Labels of Gaussians used to distingush 
%     them by colors. If y is not given then y = 1:ncomp is used.
%  
%  options.p [1 x ncomp] Value of p.d.f on the draw isolines.
%   If not given then p is computed to make non-overlapping isolines.
%  options.fill [int] If 1 then ellipses are filled (default 0).
%
% Output:
%  h [1 x nobjects] Handles of used graphics objects.
%
% Example:
%  data = load('riply_trn');
%  model = mlcgmm( data );
%  figure; hold on;
%  ppatterns(data);
%  pgauss( model );
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 23-aug-2004, VF, uses model.y to color plots in 1D case
% 30-apr-2004, VF

[dim,ncomp]=size(model.Mean);

if nargin < 2, options=[]; else options=c2s(options); end
if ~isfield( options, 'fill'), options.fill= 0; end
if ~isfield( options, 'interp'), options.interp = 30; end
if ~isfield(model,'y'), model.y = 1:ncomp; end

oldhold=ishold;
hold on;
h=[];

if dim == 1,
 
  % univariate variances can be given as a vector
  if size(model.Cov,1) ~= size(model.Cov,2), 
    model.Cov = reshape(model.Cov,1,1,ncomp); 
  end

  a = axis;
  x = linspace(a(1),a(2),options.interp*3);
  for i=1:ncomp,
    y = pdfgauss(x, model.Mean(i),model.Cov(:,:,i));
    h = [h, plot(x,y,marker_color(model.y(i)))];
  end
  
else

 % computes isolines values automatically
 if ~isfield( options, 'p'), 
  
  minr=inf;
  
  if ncomp > 1,
  
    for i=1:ncomp-1,
      for j=i+1:ncomp,
        distr.Mean=[model.Mean(:,i),model.Mean(:,j)];
        distr.Cov(:,:,1)=model.Cov(:,:,i);
        distr.Cov(:,:,2)=model.Cov(:,:,j);
        fld = androrig( distr );
      
        if minr > fld.r, 
          minr = fld.r; 
          minCov = 0.5*(model.Cov(:,:,j)+model.Cov(:,:,i));
        end
      end
    end

    for i=1:ncomp,
      options.p(i) = exp(-0.5*(minr*0.95)^2)/...
                     (2*pi*sqrt(det(model.Cov(:,:,i))));
    end
  
  else
    minr = 1;
    options.p = exp(-0.5*(minr*0.95)^2)/...
                     (2*pi*sqrt(det(model.Cov)));
  end

 elseif length(options.p) == 1 , 
  options.p = options.p*ones(ncomp,1); 
 end

 for i=1:ncomp,

  r = sqrt( -2*log(options.p(i)*2*pi*sqrt(det(model.Cov(:,:,i)))) );
  
  [x,y] = ellips(model.Mean(:,i),inv(model.Cov(:,:,i)),...
      r,options.interp );
  
  if options.fill,
    h=[h,fill(x,y,marker_color(model.y(i)))];
    h=[h,plot(model.Mean(1,i),model.Mean(2,i),'k+')];
    h=[h,text(model.Mean(1,i),model.Mean(2,i),['gauss ' num2str(i)])];
    set(h(end),'HorizontalAlignment','center', ...
      'VerticalAlignment','bottom','fontsize',12 );
    T = [[options.p(i);length(x)] [x;y]];
    tmp1=clabel(T,'fontsize',12);
    h = [h,tmp1(:)'];
  else 
    h=[h,plot(x,y,marker_color(model.y(i)))];
    h=[h,plot(model.Mean(1,i),model.Mean(2,i),['+' marker_color(model.y(i))])];
    h=[h,text(model.Mean(1,i),model.Mean(2,i),['gauss ' num2str(i)])];
    set(h(end),'HorizontalAlignment','center', ...
      'VerticalAlignment','bottom','fontsize',12,...
      'color',marker_color(model.y(i)));
    T = [[options.p(i);length(x)] [x;y]];
    tmp1=clabel(T,'fontsize',12);
    set(tmp1,'Color',marker_color(model.y(i)));
    h = [h,tmp1(:)'];
  end

  drawnow;
 end
end

if ~oldhold, hold off; end

if nargout > 0, varargout{1} = h; end

return;
