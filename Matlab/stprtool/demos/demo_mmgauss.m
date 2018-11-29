function demo_mmgauss(action,hfigure,varargin)
% DEMO_MMGAUSS Demo on minimax estimation for Gaussian.
%
% Synopsis:
%  demo_mmgauss
% 
% Description:
%  demo_mmgauss demonstrates the minimax estimation algorithm 
%  [SH10] for bivariate Gaussian distribution. The training data 
%  is supposed to contain samples which well describing the 
%  probability distribution function (pdf), i.e., which have 
%  high value of pdf. The samples do not have to be i.i.d. in 
%  contrast to the ML estimation.
%  
%  The estimated model is visualized as an ellipsoid:
%  shape is influenced by the covariance matrix and the center
%  corresponds to the mean vector.
%  The lower (red) and upper (blue) bound on the optimal value 
%  of the optimized minimax criterion is displayed at the bottom
%  part of the window.
%
% Control:
%  Epsilon     - Stopping condition. The algorithm stops if the 
%                difference between lower and the upper bound
%                is less then the epsilon.
%               
%  Iterations  - Number of iterations after which the model 
%                is re-displayed.
%
%  FIG2EPS     - Exports figure to the PostScript file.
%  Load data   - Loads input data sample from file.
%  Create data - Invokes program for creating data sample.
%  Reset       - Resets the demo.
%  Play        - Runs the algorithm.
%  Stop        - Stops the running algorithm.
%  Step        - Performs one iteration of the algorithm.
%  Info        - Invokes the info box.
%  Close       - Closes the program.
%
% See also MMGAUSS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 2-may-2004, VF
% 19-sep-2003, VF
% 3-mar-2003, VF
% 11-june-2001, V.Franc, comments added.
% 24. 6.00 V. Hlavac, comments polished.

% Used functions: PPOINTS, ELLIPS

BORDER=0.25;           % space between window limits and the points
CENTERSIZE=10;         % size of center point
LINE_WIDTH=1;
AXIST_ADD=10;
DATA_IDENT='Finite sets, Enumeration';   % file identifier

if nargin < 1,
   action = 'initialize';
end

% what action is required ?
switch lower(action)

case 'initialize'
   % == Initialize user interface control and figure window ================

   % == Figure =============================================================
   left=0.2;
   width=0.6;
   bottom=0.1;
   height=0.8;
   hfigure=figure('Name','Minimax learning', ...
      'Visible','off',...
      'NumberTitle','off', ...
      'Units','normalized', ...
      'Position',[left bottom width height],...
      'tag','Demo_mmgauss',...
      'doublebuffer','on',...
      'backingstore','off');

   % == Axes ===============================================================
   % axes with prob.
   left=0.1;
   width=0.65;
   bottom=0.1;
   height=0.28;
   haprob=axes(...
       'Units','normalized', ...
      'NextPlot','add',...
      'Position',[left bottom width height]);
   title('blue - log p(x), red - \Sigma \alpha(x) log p(x)',...
      'Parent',haprob,...
      'VerticalAlignment','bottom',...
      'Units','normalized',...
      'HorizontalAlignment','left',...
      'Position',[0 1 0]);
   htxsteps=xlabel('step number');

   % points
   height=0.45;
   bottom=0.5;
%      'XTick',[],'YTick',[], ...
   haset=axes(...
      'Units','normalized', ...
      'NextPlot','add',...
      'Position',[left bottom width height]);
   ylabel('feature y');
   xlabel('feature x');

    % == Buttons ===========================================================
   % -- Export to EPS ---------
   width=0.1;
   left=0.75-width;
   bottom=0.95;
   height=0.04;
   hbtclose = uicontrol(...
    'Units','Normalized', ...
      'Callback','fig2eps(gcf)',...
        'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','FIG2EPS');
   %----------------------------------

   % Close button
   left=0.8;
   bottom=0.05;
   height=0.05;
   width=0.15;
   hbtclose = uicontrol(...
    'Units','Normalized', ...
      'Callback','close(gcf)',...
        'ListboxTop',0, ...
        'Position',[left bottom width height], ...
        'String','Close');

   % Info button: call stanard info box
   bottom=bottom+1.5*height;
   hbtinfo = uicontrol(...
    'Units','Normalized', ...
      'Callback','demo_mmgauss(''info'',gcf)',...
        'ListboxTop',0, ...
        'Position',[left bottom width height], ...
        'String','Info');

   % Step button: perform one adaptation step
   bottom=bottom+1.5*height;
    hbtstep = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Step', ...
      'Interruptible','off',...
      'Callback','demo_mmgauss(''step'',gcf)');

   % Stop button: stop process of adaptation
   bottom=bottom+height;
   hbtstop = uicontrol(...
    'Units','Normalized', ...
        'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Stop', ...
      'Callback','set(gcbo,''UserData'',1)',...
      'Enable','off');

   % Play button: start up adaptation
   bottom=bottom+height;
   hbtplay = uicontrol(...
    'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Play', ...
      'Callback','demo_mmgauss(''play'',gcf)');

   % Reset button: set up t = 0
   bottom=bottom+height;
    hbtreset = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Reset', ...
      'Callback','demo_mmgauss(''reset'',gcf)');

   % Create data
   bottom=bottom+1.5*height;
    hbtcreat = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Create data', ...
      'Callback','demo_mmgauss(''creatdata'',gcf)');

   % Load data
   bottom=bottom+1*height;
    hbtload = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Load data', ...
      'Callback','demo_mmgauss(''getfile'',gcf)');

   % == Edit line ==========================================================

   % epsilon
   bottom=0.95-height;
   htxeps=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'String','epsilon');
   bottom=bottom-height;
   hedeps = uicontrol(...
    'Units','normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'Style','edit',...
      'String','0.1');

   % Iterations
   bottom=bottom-1.5*height;
   htxiter=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'String','Iterations');
   bottom=bottom-height;
   hediter = uicontrol(...
    'Units','normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'Style','edit',...
      'String','1');

   % == Texts ===========================================================

   htitle1=title('No data loaded',...
      'Parent',haset,...
      'VerticalAlignment','bottom',...
      'HorizontalAlignment','left',...
      'Units','normalized',...
      'Position',[0 1 0]);

   %=====================================================================
   % Store handlers
   handlers=struct(...
      'ellipse',struct('handler',-1,'center',-1,'mi',[],'sigma',[],'t',0,'N',[]),...
      'plot1',struct('handler',-1,'topps',[],'axist',0,'time',[]),...
      'plot2',struct('handler',-1,'minps',[]),...
      'title1',htitle1,...
      'btstep',hbtstep,...
      'btstop',hbtstop,...
      'btclose',hbtclose,...
      'btplay',hbtplay,...
      'btreset',hbtreset,...
      'btinfo',hbtinfo,...
      'btload',hbtload,...
      'btcreat',hbtcreat,...
      'txsteps',htxsteps,...
      'txeps',htxeps,...
      'txiter',htxiter,...
      'aset',haset,...
      'aprob',haprob,...
      'editer',hediter,...
      'edeps',hedeps);
   set(hfigure,'UserData',handlers)

   % Reset
   demo_mmgauss('reset',hfigure);

   % Put figure on desktop
   set(hfigure,'Visible','on');
   drawnow;


case 'play'
   % == One step learning ==============================================
   h=get(hfigure,'UserData');

   % get data set
   sets=get(h.aset,'UserData');

   % are data sets loaded ?
   if isempty(sets)==1,
      return;
   end

   % disable button
   set([h.editer,h.edeps,h.btstep,h.btclose,h.btplay,...
        h.btreset,h.btinfo,h.btload,h.btcreat,h.txeps,h.txiter],...
      'Enable','off');
   % enable stop button
   set(h.btstop,'Enable','on');

   % get # of iter and epsilon
   iter=str2num(get(h.editer,'String'));
   epsilon=str2num(get(h.edeps,'String'));

   % set stop button
   set(h.btstop,'UserData',0);

   % start point for plot
   if h.ellipse.t==0 & iter > 1,
%%      [mi,sigma,solution,minp,topp]=mmln(sets.X,epsilon,1,0);
      options.eps = epsilon;
      options.tmax = 1;
      model=mmgauss(sets.X,options);
      mi = model.Mean;
      sigma = model.Cov;
      minp = model.lower_bound(end);
      topp = model.upper_bound(end);
      solution = model.exitflag;

      h.plot1.time=[1];
      h.plot2.minps=[minp];
      h.plot1.topps=[topp];
   end


   % Play - adaptation process
   play=1;
   while play==1 & get(h.btstop,'UserData')==0,

      % perform one learning step
%      [mi,sigma,solution,minp,topp,h.ellipse.N,h.ellipse.t]=...
%         mmln(sets.X,epsilon,iter,h.ellipse.t,h.ellipse.N);

   options.tmax = iter+h.ellipse.t;
   options.eps = epsilon;
   if h.ellipse.t == 0,
     model=mmgauss(sets.X,options);
   else
     init_model.t = h.ellipse.t;
     init_model.Alpha = h.ellipse.N;
     model=mmgauss(sets.X,options,init_model);
   end

   mi= model.Mean;
   sigma = model.Cov;
   solution = model.exitflag;
   minp = model.lower_bound(end);
   topp = model.upper_bound(end);
   h.ellipse.N = model.Alpha;
   h.ellipse.t = model.t;
     
     
      text=sprintf('iteration t=%d ',h.ellipse.t);
      if solution==1,
         text=strcat(text,', algorithm has converged.');
         play=0;
      else
         % add new value to plot
         h.plot1.time=[h.plot1.time,h.ellipse.t];
         h.plot2.minps=[h.plot2.minps,minp];
         h.plot1.topps=[h.plot1.topps,topp];

         % is axis to be changed ?
         if h.ellipse.t > h.plot1.axist,
            h.plot1.axist=h.ellipse.t+iter*AXIST_ADD;
            set(h.aprob,'XLim',[1 h.plot1.axist]);
         end
         set(h.plot2.handler,'XData',h.plot1.time,'YData',h.plot2.minps,'Visible','on');
         set(h.plot1.handler,'XData',h.plot1.time,'YData',h.plot1.topps,'Visible','on');

         % erase old ellipse, compute and plot new one
         set(h.ellipse.center,'XData',mi(1),'YData',mi(2),'Visible','on');
         r=sqrt(max(mahalan(sets.X,mi,sigma)));
%%%         [x,y]=ellipse(inv(sigma),50,r,mi);
         [x,y]=ellips(mi,inv(sigma),r,40);
         set(h.ellipse.handler,'XData',x,'YData',y,'Visible','on');
       end

      % comment
      set(h.txsteps,'String',text);

      % store data
      set(hfigure,'UserData',h);

      % flush it on desktop
      drawnow;
   end

   % disable button
   set([h.edeps,h.editer,h.btstep,h.btclose,h.btplay,...
        h.btreset,h.btinfo,h.btload,h.btcreat,h.txeps,h.txiter],...
      'Enable','on');
   % enable stop button
   set(h.btstop,'Enable','off');


case 'step'
   % == One step learning ==============================================
   h=get(hfigure,'UserData');

   % get data set
   sets=get(h.aset,'UserData');

   % are data sets loaded ?
   if isempty(sets)==1,
      return;
   end

   % get # of iter and epsilon
   iter=str2num(get(h.editer,'String'));
   epsilon=str2num(get(h.edeps,'String'));

   % start point for plot
   if h.ellipse.t==0 & iter > 1,
%%      [mi,sigma,solution,minp,topp]=mmln(sets.X,epsilon,1,0);
      options.eps = epsilon;
      options.tmax = 1;
      model=mmgauss(sets.X,options);
      mi = model.Mean;
      sigma = model.Cov;
      minp = model.lower_bound(end);
      topp = model.upper_bound(end);
      solution = model.exitflag;
      
      
      h.plot1.time=[1];
      h.plot2.minps=[minp];
      h.plot1.topps=[topp];
   end

   % perform one learning step

%%   [mi,sigma,solution,minp,topp,h.ellipse.N,h.ellipse.t]=...
%%      mmln(sets.X,epsilon,iter,h.ellipse.t,h.ellipse.N);

   options.tmax = iter+h.ellipse.t;
   options.eps = epsilon;
   if h.ellipse.t == 0,
     model=mmgauss(sets.X,options);
   else
     init_model.t = h.ellipse.t;
     init_model.Alpha = h.ellipse.N;
     model=mmgauss(sets.X,options,init_model);
   end
   mi= model.Mean;
   sigma = model.Cov;
   solution = model.exitflag;
   minp = model.lower_bound(end);
   topp = model.upper_bound(end);
   h.ellipse.N = model.Alpha;
   h.ellipse.t = model.t;

   text=sprintf('iteration t=%d ',h.ellipse.t);
   if solution==1,
      text=strcat(text,',solution was found.');
   else
      % add new value to plot
      h.plot1.time=[h.plot1.time,h.ellipse.t];
      h.plot2.minps=[h.plot2.minps,minp];
      h.plot1.topps=[h.plot1.topps,topp];

      % is axis to be changed ?
      if h.ellipse.t > h.plot1.axist,
         h.plot1.axist=h.ellipse.t+iter*AXIST_ADD;
         set(h.aprob,'XLim',[1 h.plot1.axist]);
      end
      set(h.plot2.handler,'XData',h.plot1.time,'YData',h.plot2.minps,'Visible','on');
      set(h.plot1.handler,'XData',h.plot1.time,'YData',h.plot1.topps,'Visible','on');

      % erase old ellipse, compute and plot new one
      set(h.ellipse.center,'XData',mi(1),'YData',mi(2),'Visible','on');
      r=sqrt(max(mahalan(sets.X,mi,sigma)));
%%%      [x,y]=ellipse(inv(sigma),50,r,mi);
      [x,y]=ellips(mi,inv(sigma),r,40);
      set(h.ellipse.handler,'XData',x,'YData',y,'Visible','on');
   end

   % comment
   set(h.txsteps,'String',text);

   % flush it on desktop
   drawnow;

   set(hfigure,'UserData',h);


case 'redraw'
   % == Redraw contents of axes ====================================

   h=get(hfigure,'UserData');                   % uicontrol handlers

   % get sets with points
   sets=get(h.aset,'UserData');
   if isempty(sets)==1,
      return;
   end

   % clears axes
   axes(h.aset);
   set(get(h.aset,'Children'),'EraseMode','normal');
   %%%   cla;
   clrchild(h.aset);
   h.ellipse.handler=plot([0],[0],...
      'Parent',h.aset,'LineWidth',LINE_WIDTH,...
      'EraseMode','xor','Color','r','Visible','off');
   h.ellipse.center=line(0,0,'Marker','+',...
      'EraseMode','xor','Color','r','MarkerSize',CENTERSIZE,'Visible','off');

   win=cmpwin(min(sets.X'),max(sets.X'),BORDER,BORDER);
%%%   axis(win);
   setaxis(h.aset,win);
   axes(h.aset);

   % plot mixture
%%   pplot(sets.X,sets.I);
%%   ppoints(sets.X,sets.I);
   ppatterns(sets.X);

   set(hfigure,'UserData',h);

case 'getfile'
   % == Invoke standard open file dialog ====================================
   % Opens file and checks if contains appropriate data, if yes than loads data.

   h=get(hfigure,'UserData');

   % change path to directory
%%   wres=what('minimax');
%%   cd(wres.path);

   [name,path]=uigetfile('*.mat','Open file');
   if name~=0,
      file.pathname=strcat(path,name);
      file.path=path;
      file.name=name;
%      if checkdat(file.pathname,DATA_IDENT,2,[0])==1,
      if check2ddata(file.pathname),
         set(h.btload,'UserData',file);
         demo_mmgauss('loadsets',hfigure);
      else
         errordlg('This file does not contain required data.','Bad file','modal');
      end
   end


case 'loadsets'
   % == Load sets ==================================================================
   % Get file name from the pop up menu according to menu pointer.

   h=get(hfigure,'UserData');

   % Get file name with sets
   file=get(h.btload,'UserData');

   % Load sets
   sets=load(file.pathname);
   sets.I = sets.y;
   sets.K = size(sets.X,2);
   sets.N = 2;

   % store loaded sets
   set(h.aset,'UserData',sets);

   % call reset
   demo_mmgauss('reset',hfigure);

   % call redraw
   demo_mmgauss('redraw',hfigure);

   drawnow;


case 'reset'
   % == Reset adaptation process, set up t=0 ================

   h=get(hfigure,'UserData');                     % get handlers

   % get file
   file=get(h.btload,'UserData');

   % get data set
   sets=get(h.aset,'UserData');

   % zeroes parameters of the separation line
   h.ellipse.mi=[];
   h.ellipse.sigma=[];
   h.ellipse.t=0;
   h.ellipse.N=[];
   h.plot1.topps=[];
   h.plot2.minps=[];
   h.plot1.time=[];
   h.plot1.axist=0;

   % comment
   text=sprintf('iteration t=0 ');
   set(h.txsteps,'String',text);

   % hide ellipse
   if h.ellipse.handler==-1,
      h.ellipse.handler=plot([0],[0],...
         'Parent',h.aset,'LineWidth',LINE_WIDTH,...
         'EraseMode','xor','Color','k','Visible','off');
      h.ellipse.center=line(0,0,'Marker','x',...
         'EraseMode','xor','Color','k','MarkerSize',CENTERSIZE,'Visible','off');
   else
      set(h.ellipse.handler,'Visible','off');
      set(h.ellipse.center,'Visible','off');
   end % if h.ellipse.handler~=-1

   % clear axes prob.
   axes(h.aprob);
   cla;
   axis auto;
   h.plot1.handler=plot([0],[0],'b','Parent',h.aprob,...
      'EraseMode','background','Visible','off');
   h.plot2.handler=plot([0],[0],'r','Parent',h.aprob,...
      'EraseMode','background','Visible','off');

   % set up handlers and flush queue with graph. objects
   set(hfigure,'UserData',h);

   % creat comment
   axes(h.aset);
   if isempty(sets)==0,
      set(h.title1,'String',sprintf('Data file: %s, %d vectors',...
          file.name,sum(sets.K))); 
   else
      set(h.title1,'String','No data loaded');

      pos=get(h.aset,'Position');
      fsize=min(pos(3),pos(4))/8;
      axis([-1 1 -1 1]);
      builtin('text',0,0,'Press ''Load data'' button.',...
         'HorizontalAlignment','center',...
         'FontUnits','normalized',...
         'Clipping','on',...
         'FontSize',fsize);
      builtin('text',0,-fsize*2,...
         'Load sample data from ../toolboxroot/data/mm\_samples/ ',...
         'HorizontalAlignment','center',...
         'FontUnits','normalized',...
         'Clipping','on',...
         'FontSize',fsize*0.65);
   end

   drawnow;

case 'creatdata'
   % == Invoke data set creator ============================================
%%   creatset('finite',1,'demo_mmgauss','created',hfigure);
   createdata('finite',1,'demo_mmgauss','created',hfigure);

case 'created'
   % == Load new created data set ===========================================

   % get handler and make this figure active
   figure(hfigure);
   h=get(hfigure,'UserData');

   % get file name
   path=varargin{1};
   name=varargin{2};
   pathname=strcat(path,name);

%%   if checkdat(pathname,DATA_IDENT,2,[0])==1,
   if check2ddata(pathname),
      file.pathname=pathname;
      file.path=path;
      file.name=name;
      set(h.btload,'UserData',file);
      demo_mmgauss('loadsets',hfigure);
   else
      errordlg('This file does not contain required data.','Bad file','modal');
   end

case 'info'
   % == Call standard Matlab`s info box =========================================
   helpwin(mfilename);
end % switch

%%%%%%%%%%%%%%%%%%%%
function []=clrchild(handle)
% function []=clraxis(handle)
%
% CLRCHILD clears children of an object with the given handle.
%
% See also DELETE.
%
% Statistical Pattern Recognition Toolbox, Vojtech Franc, Vaclav Hlavac
% (c) Czech Technical University Prague, http://cmp.felk.cvut.cz
% Written Vojtech Franc (diploma thesis) 02.01.2000
% Modifications
% 24. 6.00 V. Hlavac, comments polished.

delete(get(handle,'Children'));

return;

function [win]=cmpwin(mins,maxs,xborder,yborder)
%
%  [win]=cmpwin(mins,maxs,xborder,yborder)
%
% CMPWIN computes appropriate size of the axes.
%

% Statistical Pattern Recognition Toolbox, Vojtech Franc, Vaclav Hlavac
% (c) Czech Technical University Prague, http://cmp.felk.cvut.cz
% Written Vojtech Franc (diploma thesis) 02.01.2000
% Modifications
% 24. 6.00 V. Hlavac, comments polished.

dx=max( (maxs(1)-mins(1)), 1 )*xborder;
dy=max( (maxs(2)-mins(2)), 1 )*yborder;

%x1=floor(mins(1)-dx);
%x2=ceil(maxs(1)+dx);
%y1=floor(mins(2)-dy);
%y2=ceil(maxs(2)+dx);
x1=(mins(1)-dx);
x2=(maxs(1)+dx);
y1=(mins(2)-dy);
y2=(maxs(2)+dx);

win=[x1 x2 y1 y2];

return;

function []=setaxis(handle,rect)
% function []=setaxis(handle,rect)
%
% SETAXIS sets scaling for the x- and y-axes
%   on the plot with a given handle.
%
% See also AXIS.
%
% Statistical Pattern Recognition Toolbox, Vojtech Franc, Vaclav Hlavac
% (c) Czech Technical University Prague, http://cmp.felk.cvut.cz
% Written Vojtech Franc (diploma thesis) 02.01.2000
% Modifications

set(handle,'XLim',rect(1:2));
set(handle,'YLim',rect(3:4));

if size(rect,2)>=6,
   set(handle,'ZLim',rect(5:6));
end

return;
