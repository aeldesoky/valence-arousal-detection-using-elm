function varargout=pboundary(model,options)
% PBOUNDARY Plots decision boundary of given classifier in 2D.
%
% Synopsis:
%  h = pboundary(model)
%  h = pboundary(model,options)
%
% Description:
%  This function plots decision boundary of given classifier in 
%  2-dimensional feature space. The classification function
%  must be specified in the field model.eval. The decision
%  bounary is interpolated from the response of the classifier
%    y = feval( model.fun, X, model).
%
% Input:
%  model [struct] Model of classifier.
%   .fun [string] Classification function.
%  
%  options [struct] Controls visualization:
%   .gridx [1x1] Sampling density in x-axis (default 200).
%   .gridy [1x1] Sampling density in y-axis (default 200).
%   .line_style [string] Used line-style to plot decision boundary.
%   .fill [1x1] If 1 then the class regions are filled. 
%  
% Output:
%  h [1 x nobjects] Handles of used graphics objects.
%
% Example:
%  data = load('riply_trn');
%  figure; 
%  ppatterns(data);
%  pboundary( knnrule(data,1) );
%
% See also 
%  PPATTERNS, PLINE.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 1-may-2004, VF
% 19-may-2003, VF

% process inputs
if nargin < 1, error('Not enough input arguments.'); end
model=c2s(model);

if nargin < 2, options=[]; else options=c2s(options); end
if ~isfield(options,'fill'), options.fill=0; end
if ~isfield(options,'gridx'), options.gridx = 200; end
if ~isfield(options,'gridy'), options.gridy = 200; end
if ~isfield(options,'line_style'), options.line_style = 'k'; end

% get hold option
old_hold = ishold;
hold on;

%
V = axis;
dx = (V(2)-V(1))/options.gridx;
dy = (V(4)-V(3))/options.gridy;

[X,Y] = meshgrid(V(1):dx:V(2),V(3):dy:V(4));

% make testing points
tst_data=[reshape(X',1,prod(size(X)));reshape(Y',1,prod(size(Y)))];

% classify points
D = feval( model.fun, tst_data, model );

h = plot_boundary( D, V(1):dx:V(2), V(3):dy:V(4), ...
    options.fill, options.line_style );

if ~old_hold, hold off; end

% return handles if required
if nargout > 0, varargout{1} = h; end

return;

%-------------------------------------------------------
function h = plot_boundary( L, X_pos, Y_pos, fill_regions, linestyle )
% Plots decision boudary. 
%

dx=X_pos(2)-X_pos(1);
dy=Y_pos(2)-Y_pos(1);
m = length( X_pos );
n = length( Y_pos );

Z = NaN*ones( m+2, n+2 );
num_classes = max(L(:));

% mask=fspecial('gauss',[5 5],1) ; % fspecial is from images toolbox 
mask = [0.0030    0.0133    0.0219    0.0133    0.0030;
        0.0133    0.0596    0.0983    0.0596    0.0133;
        0.0219    0.0983    0.1621    0.0983    0.0219;
        0.0133    0.0596    0.0983    0.0596    0.0133;
        0.0030    0.0133    0.0219    0.0133    0.0030];

h = [];      
for i = 1:num_classes,
  
  A=L;
  A(find(L==i))=1;
  A(find(L~=i))=-1;
    
  A = reshape( A', m, n );
  A = filter2(mask, A);
   
  Z(2:end-1,2:end-1) = A;  
  
  [cc,tmp_h] = contour([X_pos(1)-dx,X_pos,X_pos(end)+dx],...
                   [Y_pos(1)-dy,Y_pos,Y_pos(end)+dy],Z',[-0 0],linestyle);  
  h = [h tmp_h(:)'];
  
  if fill_regions,
   while ~isempty(cc)
     len = cc(2,1);
     tmp_h = fill(cc(1,2:len+1),cc(2,2:len+1),marker_color(i));
     h = [h tmp_h(:)'];
     cc(:,1:len+1) = [];
   end  
  end
end

return;

