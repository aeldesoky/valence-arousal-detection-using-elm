function varargout=ppatterns(data,arg1,arg2)
% PPATTERNS Plots pattern as points in feature space.
% 
% Synopsis:
%  ppatterns(data,marker_size)
%  ppatterns(data,'num')
%  ppatterns(X,marker,marker_size)
%  ppatterns(X,y)
%  ppatterns(X,y,marker_size)
%  ppatterns(X,y,'num')
%
% Description:
%  ppatterns(data,marker_size) plots data.X as points
%   distinguished by marker and its color according to 
%   given labels data.y. The marker size can be prescribed.
%
%  ppatterns(data,'num') plots data.X in distinguished 
%   by numbers and colors according to given labels data.y. 
%   The marker size can be determined by argument marker_size.
%
%  ppatterns(X,marker,marker_size) plots data X. Marker type
%   can be determined by argument marker. The marker size can 
%   be determined by argument marker_size.
%
%  ppatterns(X,y,...) instead of structure data, which contains
%   items X and y these can enter the function directly.
%
%  If dimension of input data is greater than 3 then
%  only first 3 dimensions are assumed and data are plotted 
%  in 3D space.
%
% Output:
%  H [struct] Handles of used graphical objects.
%
% Example:
%  data = load('riply_trn');  
%  figure; ppatterns(data);
%  figure; ppatterns(data,'num');
%  figure; ppatterns(data.X,'xk',10);
%  
% See also 
%  PLINE.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 25-may-2004, VF
% 11-mar-2004, VF, 
% 5-oct-2003, VF, returns handles
% 12-feb-2003, VF, 1D, 3D added
% 7-jan-2003, VF, created

oldhold = ishold;
hold on;

% takes care of the case when X,y is used insted of structure data
% ppatterns(X,y,...) -> ppatterns(data,...)
if nargin > 1 & ~isstruct(data) & size(data,2)==length(arg1) & ~isstr(arg1),
  data.X=data;
  data.y=arg1;
  if nargin >= 3, H=ppatterns(data,arg2); else H=ppatterns(data); end 
  if nargout >= 1, varargout{1}=H; end
  return;
end

% ppatterns(data) or ppatterns(data,marker_size) 
if isstruct(data) == 1 & (nargin < 2 | isstr(arg1)==0),
  
  if nargin < 2, marker_size = 6; else marker_size = arg1; end
  
  H = [];
  for i = min(data.y):max(data.y),
   
    inx = find(data.y == i);
    if ~isempty(inx),
      
      if size(data.X,1)==1,
        h = plot(data.X(1,inx),zeros(1,length(inx)),marker_type(i));
      elseif size(data.X,1)==2,
        h = plot(data.X(1,inx),data.X(2,inx),marker_type(i));
      else
        h = plot3(data.X(1,inx),data.X(2,inx),data.X(3,inx),marker_type(i));
      end
               
      set(h,'Color',marker_color(i));
      set(h,'MarkerSize',marker_size);
      H = [H, h];
    end
  end
% ppatterns(data,marker)
elseif isstruct(data) == 1 & nargin == 2 & isstr(arg1)==1 & strcmpi(arg1,'num'),

  marker_size = 12;
  H_Points = [];
  H_Num = [];
  for i = min(data.y):max(data.y),
    inx = find(data.y==i);
    if ~isempty(inx),
      if size(data.X,1)==1,
        h = plot(data.X(1,inx),zeros(1,length(inx)),'o');
      elseif size(data.X,1)==2,
        h = plot(data.X(1,inx),data.X(2,inx),'o');
      else
        h = plot3(data.X(1,inx),data.X(2,inx),data.X(3,inx),'o');
      end
       
      set(h,'Color',marker_color(i));
      set(h,'MarkerSize',marker_size);
      
      H_Points = [H_Points, h ];
         
      if size(data.X,1)==1,
        h = text(data.X(1,inx),zeros(1,length(inx)),num2str(i));
      elseif size(data.X,1)==2,
        h = text(data.X(1,inx),data.X(2,inx),num2str(i));
      else
        h = text(data.X(1,inx),data.X(2,inx),data.X(3,inx),num2str(i));
      end

      set(h,'HorizontalAlignment','center');
      set(h,'VerticalAlignment','middle');
      set(h,'Color',marker_color(i));
      set(h,'FontSize',marker_size-2);
      
      H_Num = [H_Num, h(:)'];
    end
  end
  H = [H_Points, H_Num];
else
  if nargin < 2, marker = 'kx'; else marker = arg1; end
  if nargin < 3, marker_size = 6; else marker_size = arg2; end
  if size(data,1)==1,
    h = plot(data(1,:),zeros(1,size(data,2)),marker);
  elseif size(data,1)==2,
    h = plot(data(1,:),data(2,:),marker);
  else
    h = plot3(data(1,:),data(2,:),data(3,:),marker);
  end
  
  H = h;
  
  set(h,'MarkerSize',marker_size);
end

if oldhold,
  hold on;
else
  hold off;
end

if nargout>=1, varargout{1} = H; end

return;


