function varargout=psvm(model,options)
% PSVM Plots decision boundary of binary SVM classifier.
%
% Synopsis:
%  h = psvm(...)
%  psvm(model)
%  psvm(model,options)
%
% Description:
%  This function samples the Support Vector Machiones (SVM) decision 
%  function f(x) in 2D feature space and interpolates isoline 
%  width f(x)=0. The isolines f(x)=+1 and f(x)=-1 are plotted as well. 
%
% Input:
%  model [struct] Model of binary SVM classifier:
%   .Alpha [1 x nsv] Weights of training data.
%   .b [real] Bias of decision function.
%   .sv.X [dim x nsv] Support vectors.
%   .options.ker [string] Kernel function identifier.
%      See 'help kernel' for more info.
%   .options.arg [1 x nargs] Kernel argument(s).
%
% options [struct] Controls apperance:
%  .background [1x1] If 1 then backgroud is colored according to 
%    the value of decision function (default 0).
%  .sv [1x1] If 1 then the support vectors are marked (default 1).
%  .sv_size [1x1] Marker size of the support vectors.
%  .margin [1x1] If 1 then margin is displayed (default 1).
%  .gridx [1x1] Sampling in x-axis (default 25).
%  .gridy [1x1] Sampling in y-axis (default 25).
%  .color [int] Color of decision boundary (default 'k').
%
% Output:
%  h [struct] Handles of used graphical objects.
%
% Example:
%  data = load('riply_trn');  
%  model = smo( data, struct('ker','rbf','arg',1,'C',10) );
%  figure;  ppatterns(data);  
%  psvm( model, struct('background',1) );
%
% See also 
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 25-may-2004, VF
% 10-may-2004, VF
% 5-oct-2003, VF, returns handles
% 14-Jan-2003, VF
% 21-oct-2001, V.Franc
% 16-april-2001, V. Franc, created

% Process input arguments 
%------------------------------------------------
if nargin < 2, options=[]; else options=c2s(options); end
if ~isfield(options,'background'), options.background = 0; end
if ~isfield(options,'sv'), options.sv = 1; end
if ~isfield(options,'margin'), options.margin = 1; end
if ~isfield(options,'gridx'), options.gridx = 25; end
if ~isfield(options,'gridy'), options.gridy = 25; end
if ~isfield(options,'sv_size'), options.sv_size = 12; end
if ~isfield(options,'color'), options.color = 'k'; end

% variable for handles
h = [];

% get axis
a = axis;
old_hold = ishold;
hold on;

% plot Support Vectos
if options.sv,
  h.sv = ppatterns(model.sv.X,'ok' ,options.sv_size);
end

% limits of current figure
xmin=a(1);
xmax=a(2);
ymin=a(3);
ymax=a(4);
  
% makes grid 
[X,Y] = meshgrid(xmin:(xmax-xmin)/options.gridx:xmax,...
                 ymin:(ymax-ymin)/options.gridy:ymax);

% generate samples
tst_data=[reshape(X',1,prod(size(X)));reshape(Y',1,prod(size(Y)))];

% classify points
[pred_labels, dec_fun] = svmclass(tst_data,model);

% compute color limits
l=(-min(dec_fun)+max(dec_fun))/2;

% reshape dec_fun
Z = reshape(dec_fun,size(X,1),size(X,2))';

% colors background 
if options.background,
  h.background = pcolor(X,Y,Z);
end

% smooth shading
shading interp;

% plots decision boundary
[dummy,h.boundary] = contour(X,Y,Z,[0,0],options.color);

% plots margins
if options.margin,
   [dummy,h.margin_plus] = contour(X,Y,Z,[1,1],[options.color,'--']);
   [dummy,h.margin_minus] = contour(X,Y,Z,[-1,-1],[options.color,'--']);
end

% set color limits and colormap
if options.background,
  set(h.background, 'LineStyle','none' );
  set(gca,'Clim',[-l l]);
  
  % creates colormap 
  g=gray(64);
  cmp=[g(33:end,:);flipud(g(33:end,:))];
  cmp(1:32,1)=cmp(1:32,1)/2;
  cmp(1:32,3)=cmp(1:32,3)/2;
  cmp(33:end,3)=cmp(33:end,3)/2;
  colormap(cmp)
end

if ~old_hold, hold off; end

if nargout >= 1, varargout{1} = h; end

return;
% EOF