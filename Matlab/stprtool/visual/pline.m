function varargout=pline(arg1,arg2,arg3)
% PLINE Plots line in 2D.
%
% Synopsis:
%  h=pline(W,b)
%  h=pline(W,b,line_style)
%  h=pline(model)
%  h=pline(model,options)
%
% Description:
%  h=pline(W,b) plots the line in 2D space described implicitely as
%   W'*x + b = 0 ,
%  where W, x are vectors [2x1] and b is scalar or explicitely as
%   y = W*x + b ,
%  where W, x and b are scalars.
%
%  h=pline(W,b,line_style) defines parameter line_style of plot 
%   function (default 'k-').
%  
%  h=pline(model) parameters of the line are given in structure 
%   model with fields model.W and model.b.
%
%  h=pline(model,options) argument options controls apperance
%   of the plotted line; options.win [left right top bottom] 
%   determines window to which the line is plotted and 
%   options.line_style is described above.
%
% Output:
%  h [1x1] handle of plotted line.
%
% Example:
%
% Plot horizontal and vertical axes with dashed line:
%  figure; hold on; axis([-1 1 -1 1]);
%  pline(inf,0,'--'); 
%  pline(0,0,'--');
%
% Plot Fisher linear discriminat for Riply's data set:
%  data = load('riply_trn');
%  model = lfld( mlcgmm(data));
%  figure; ppatterns(data);
%  pline(model);
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 29-apr-2004, VF
% 13-july-2003, VF
% 20-jan-2003, VF
% 8-jan-2003, VF, A new coat.

oldhold=ishold;
hold on;

if nargin >= 2 & isstruct(arg2)==0,
   model.W = arg1;
   model.b = arg2;
   if nargin > 2, options.line_style = arg3; else options = []; end
else
  model=c2s(arg1);
  if nargin < 2, options = []; else options=c2s(arg2); end
end

if ~isfield(options,'win'), options.win = axis; end
if ~isfield(options,'line_style'), options.line_style = 'k-'; end

if length(model.W)==1,
   model.W = [model.W; -1];
end

[x1,y1,x2,y2,in]=clipline(model.W,model.b,options.win);

h = plot([x1,x2],[y1,y2],options.line_style);
if nargout>0, varargout{1}= h; end

if ~oldhold, hold off; end

return;

