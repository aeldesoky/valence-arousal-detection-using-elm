function result = demo_linclass(action,hfigure,varargin)
% DEMO_LINCLASS Demo on the algorithms learning linear classifiers.
%
% Synopsis:
%  demo_linclass
%
% Description:
%  DEMO_LINCLASS demonstrates use of the algorithms which find 
%  linear decision hyperplane between two (dichotomy)
%  vector sets. The demo requires 2D input training data.
%
%  The program vizualizes found hyperplane in the current 
%  algorithm step. The missclassified vector used by the
%  demonstrated iterative algorithms for update is vizualized
%  as well. Text description of the found solution is
%  printed at the bottom part of window.
%
%  Following algorithms can be tested:
%
%  Perceptron  - Perceptron learning rule (see 'help perceptron').
%  Kozinec     - Kozinec's algorithm (see 'help ekozinec', eps=-1). 
%  e-Kozinec   - Kozinec's algorithm finding eps-optimal hyperplane
%                (see 'help ekozinec', eps > 0).
%  Linear SVM   - Linear Supprot Vector Machines for separable data
%                (see 'help smo', C=inf, ker='linear').
%
% Control:
%  Algorithm  - Dselects algorithm for testing.
%  Epsilon    - Input parameter of 'ekozinec' algorithm 
%               (see 'help ekozinec').
%  Iterations - Number of iterations in one step.
%  Animation  - Enables/dissables animation - smooth changeover 
%               between two algorithm states.
%
%  FIG2EPS     - Exports screen to the PostScript file.
%  Load data   - Loads input training data from file.
%  Create data - Invokes program for creating training data.
%  Reset       - Sets the algorithm to the initial state.
%  Play        - Runs the algorithm.
%  Stop        - Stops the running algorithm.
%  Step        - Perform one step of the algorithm.
%  Info        - Invoke the info box.
%  Close       - Close the program.
%
% See also 
%  PERCEPTRON, EKOZINEC, SVM.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 19-sep-2003, VF
% 17-Feb-2003, VF
% 24. 6.00 V. Hlavac, comments polished.
% 11-dec-2000 V. Franc, a little increasing of code readibility
% 15-dec-2000

LINE_WIDTH=1;        % width of separation line
ANIM_STEPS=10;       % number of steps during the line animation
BORDER=0.2;          % minimal space between axis and points
%DATA_IDENT='Finite sets, Enumeration';  % file identifier
ALGOS=['Perceptron';'Kozinec   ';'e-Kozinec ';  'LinearSVM '];
WRONGPOINT_SIZE = 13;

% if number of arguments is less then 1, that means first call of this
% function. Every other calls set up at least argument action
if nargin < 1,
   action = 'initialize';
end

% what action is required ?
switch lower(action)

case 'initialize'
   % == Initialize user interface control and figure =======

   % == Figure =============================================
   left=0.1;
   width=0.8;
   bottom=0.1;
   height=0.8;
   hfigure=figure('Name','Linear discriminant function', ...
      'Visible','off',...
    'NumberTitle','off', ...
      'Units','normalized', ...
      'Position',[left bottom width height],...
      'Units','normalized', ...
      'tag','Demo_Linclass',...
      'RendererMode','manual');

   % == Axes =========================================
   left=0.1;
   width=0.65;
   bottom=0.35;
   height=0.60;
   haxes1=axes(...
       'Units','normalized', ...
      'Box','on', ...
      'DrawMode','fast',...
      'UserData',[],...
      'Position',[left bottom width height]);
   xlabel('feature x_1');
   ylabel('feature x_2');

   % == Comment window =================================
   % Comment Window frame
   bottom=0.05;
   height=0.2;
   uicontrol( ...
        'Style','frame', ...
        'Units','normalized', ...
        'Position',[left bottom width height], ...
        'BackgroundColor',[0.5 0.5 0.5]);

   % Text label
   uicontrol( ...
        'Style','text', ...
        'Units','normalized', ...
        'Position',[left height-0.01 width 0.05], ...
        'BackgroundColor',[0.5 0.5 0.5], ...
        'ForegroundColor',[1 1 1], ...
        'String','Comment Window');

   % Edit window
   border=0.01;
   hconsole=uicontrol( ...
        'Style','edit', ...
        'HorizontalAlignment','left', ...
        'Units','normalized', ...
        'Max',10, ...
        'BackgroundColor',[1 1 1], ...
        'Position',[left+border bottom width-2*border height-0.05], ...
        'Enable','inactive',...
        'String','');


    % == Buttons ===========================================
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
      'Callback','demo_linclass(''info'',gcf)',...
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
      'Callback','demo_linclass(''step'',gcf)');

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
      'Callback','demo_linclass(''play'',gcf)');

   % Reset button: set up t = 0
   bottom=bottom+height;
    hbtreset = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Reset', ...
      'Callback','demo_linclass(''reset'',gcf)');

   % Creat data
   bottom=bottom+1.5*height;
    hbtcreat = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Create data', ...
      'Callback','demo_linclass(''creatdata'',gcf)');

   % Load data
   bottom=bottom+1*height;
    hbtload = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Load data', ...
      'Callback','demo_linclass(''getfile'',gcf)');

   % == Popup menus ======================================
   % Pop up menu for the selection between algorithms
    % title
   bottom=0.95-height;
   htxalgo=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width height], ...
      'String','Algorithm');
   % popup menu
   bottom=bottom-height;
   hpualgo=uicontrol( ...
      'Style','popup', ...
      'Units','normalized', ...
      'CallBack','demo_linclass(''epshandler'',gcf)',...
      'Position',[left bottom width height], ...
      'String',ALGOS);


   % == Edit line ========================================
   % epsilon
   bottom=bottom-1.2*height;
   htxeps=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'Enable','off',...
      'String','epsilon');
   bottom=bottom-height;
   hedeps = uicontrol(...
    'Units','normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'Style','edit',...
      'Enable','off',...
      'CallBack','demo_linclass(''epshandler'',gcf)',...
      'String','1e-2');

   % # of iterations
   bottom=bottom-1.1*height;
   htxiter=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'String','Iterations');

   bottom=bottom-0.9*height;
   hediter = uicontrol(...
    'Units','normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'Style','edit',...
      'CallBack','demo_linclass(''iterhandler'',gcf)',...
      'String','1');

   % == Check boxes ==============================================
   % Make check box to determine if a line will be drawn in one 
   % step or smooth plot.
   bottom=bottom-height*1.2;
    hxbanim = uicontrol(...
    'Style','checkbox', ...
    'Units','normalized', ...
    'ListboxTop',0, ...
    'Position',[left bottom width height], ...
    'String','Animation');

    % ============================================================
   % Store handlers
   handlers=struct(...
      'line',struct('handler',-1,'alpha',0,'theta',0,'t',0),...
      'btstep',hbtstep,...
      'btstop',hbtstop,...
      'btclose',hbtclose,...
      'btplay',hbtplay,...
      'btreset',hbtreset,...
      'btinfo',hbtinfo,...
      'btload',hbtload,...
      'btcreat',hbtcreat,...
      'pualgo',hpualgo,...
      'console',hconsole,...
      'editer',hediter,...
      'edeps',hedeps,...
      'txeps',htxeps,...
      'axes1',haxes1,...
      'xbanim',hxbanim);
   set(hfigure,'UserData',handlers)

   % Reset
   demo_linclass('reset',hfigure);

   % Put figure on desktop
   set(hfigure,'Visible','on');
   drawnow;

case 'iterhandler'
   % == Handler for edit line Iterations ===============
   h=get(hfigure,'UserData');

   iter=round(str2num(get(h.editer,'String')));
   if isempty(iter) | iter < 1, iter=1; end
   set(h.editer,'String',num2str(iter));


case 'epshandler'
   % == Handler for edit line Epsilon =======================
   h=get(hfigure,'UserData');

   % if algorithm e-Kozinec is selected then ...
   if get(h.pualgo,'Value')==3,
      set(h.edeps,'Enable','on');
      set(h.txeps,'Enable','on');

      epsil=str2num(get(h.edeps,'String'));
      if epsil < 0,
         epsil=1;
         set(h.edeps,'String',num2str(epsil));
      end
   else
      set(h.edeps,'Enable','off');
      set(h.txeps,'Enable','off');
   end


case 'creatdata'
   % == Invoke data set creator ================================
   createdata('finite',2,'demo_linclass','created',hfigure);


case 'created'
   % == Load new created data set =============================

   % get handler and make this figure active
   figure(hfigure);
   h=get(hfigure,'UserData');

   % get file name
   path=varargin{1};
   name=varargin{2};
   pathname=strcat(path,name);

%   if checkdat(pathname,DATA_IDENT,2,[0 0])==1,
   if check2ddata(pathname),
      file.pathname=pathname;
      file.path=path;
      file.name=name;
      set(h.btload,'UserData',file);
      demo_linclass('loadsets',hfigure);
      demo_linclass('reset',hfigure);
   else
      errordlg('This file does not contain required data.','Bad file','modal');
   end


case 'getfile'
   % == Invoke standard open file dialog ===================
   % Opens file and checks if contains appropriate data, 
   % if yes loads data.

   h=get(hfigure,'UserData');

   [name,path]=uigetfile('*.mat','Open file');
   if name~=0,
      file.pathname=strcat(path,name);
      file.path=path;
      file.name=name;
      if check2ddata(file.pathname),
%      if checkdat(file.pathname,DATA_IDENT,2,[0 0])==1,
         set(h.btload,'UserData',file);
         demo_linclass('loadsets',hfigure);
         demo_linclass('reset',hfigure);
      else
         errordlg('This file does not contain required data.','Bad file','modal');
      end
   end


case 'redraw'
   % == Redraw points in axes ==================================
   h=get(hfigure,'UserData');           % uicontrol handlers

   % get point sets
   sets=get(h.axes1,'UserData');
   if isempty(sets)==1,
      return;
   end

   % clears axes
   set(get(h.axes1,'Children'),'EraseMode','normal');
   clrchild(h.axes1);

   h.line.handler=line('EraseMode','xor','Color','k','Visible','off','Parent',h.axes1);
   set(hfigure,'UserData',h);                   % uicontrol handlers

%%%   pplot(sets.X,sets.I);
%%   ppoints(sets.X,sets.I);
   ppatterns(sets);
   demo_linclass('drawline',hfigure,h.line.theta,h.line.alpha);

   drawnow;


case 'loadsets'
   % == Load data sets ========================================
   % Get file name from the pop up menu according to menu pointer. 
   % Than clear axes,load new file and appear the points from the file.

   h=get(hfigure,'UserData');

   % Clear axes
   clrchild(h.axes1);

   set(h.axes1, ...
      'Box','on', ...
      'DrawMode','fast' );
   xlabel('feature x_1');
   ylabel('feature x_2');

   % No line
   h.line.handler=-1;
   set(hfigure,'UserData',h);

   % Get file name with sets
   file=get(h.btload,'UserData');

   % Load sets
   sets=load(file.pathname);
   %%%
   sets.I=sets.y;
   sets.N=2;
   for ii=1:max(sets.y),
     sets.K(ii)=length(find(sets.y==ii));
   end

   % store loaded sets
   set(h.axes1,'UserData',sets);

   % set axes according to current point set
   win=cmpwin(min(sets.X'),max(sets.X'),BORDER,BORDER);
   setaxis(h.axes1,win);

   axes(h.axes1);

   % plots points
%%   ppoints(sets.X,sets.I);
   ppatterns(sets);

   drawnow;
   
case 'play'
   % == Start up the adaptation process ============================

   % Get handle to data.
   h=get(hfigure,'UserData');               

   if h.line.handler==-1,
      return;
   end

   % Check if data are loaded.
   sets=get(h.axes1,'UserData');
   if isempty(sets)==1,
      return;
   end

   % Disable and enable buttons.
   set([h.btinfo h.btstep h.btclose h.btplay h.btreset h.btload h.btcreat ...
       h.pualgo h.editer],'Enable','off');
   set(h.btstop,'Enable','on');

   set(h.btstop,'UserData',0);

   h.stop = 0;
   set(hfigure,'UserData',h);
   
   % Play - adaptation process
   while h.stop==0 & get(h.btstop,'UserData')==0,
      demo_linclass('step',hfigure);
      h=get(hfigure,'UserData');                  
   end

   % Enable and dissable buttons.
   set([h.btinfo h.btstep h.btclose h.btplay h.btreset h.pualgo ...
      h.editer h.btload h.btcreat],'Enable','on');
   set(h.btstop,'Enable','off');


case 'step'
   % == Perform one adaptation step ======================================
   h=get(hfigure,'UserData');        % get handlers we will need...

   if h.line.handler==-1,
      return;
   end

   % get sets
   sets=get(h.axes1,'UserData');

   % no data set loaded
   if isempty(sets)==1,
      return;
   end

   [alpha,theta,solution,t]=exec(hfigure);

   if mod(h.stepcnt,2)==1 & h.line.t >0, 
     if solution ~= 1,
        demo_linclass('drawline',hfigure,theta,alpha);
     else
        text=makeinfo(t,alpha,theta,solution);
        set(h.console,'String',text );
     end 
   else
      if get(h.xbanim,'Value')==0,
        demo_linclass('drawline',hfigure,theta,alpha);
      else
        demo_linclass('animline',hfigure,theta,alpha);
      end
   
      if solution==1 | solution==-1,
        h.stop=1;
      end
      h.line.alpha = alpha;
      h.line.theta = theta;
      h.line.t = t;
   
      if solution==0 | solution ==1,
         % appear time and line parameters
         text=makeinfo(t,alpha,theta,solution);
      elseif solution==-1,
        text=sprintf('Solution does not exist.\n');
      end 
      set(h.console,'String',text );
   end


   %  store new solution
   h.stepcnt=h.stepcnt+1;
   set(hfigure,'UserData',h);


case 'animline'
   % == Smooth transition of line from old to new position ===============

   h=get(hfigure,'UserData');                     % get handlers

   % old position of line is...
   alpha2=h.line.alpha;
   theta2=h.line.theta;
   t2=h.line.t;

   % New position get from input arguments
   theta1=varargin{1};
   alpha1=varargin{2};

   if t2~=0,
      % move line
    step=1/ANIM_STEPS;
       for k=0:step:1,
       alpha=(1-k)*alpha2+k*alpha1;      % smooth transition of alpha
        theta=(1-k)*theta2+k*theta1;      % --//--                      theta

          demo_linclass('drawline',hfigure,theta,alpha);
      end
   else
      % it is first step
      demo_linclass('drawline',hfigure,theta1,alpha1); % first step
   end % if t2~=0


case 'reset'
   % == Reset adaptation process, set up zero step ================

   h=get(hfigure,'UserData');                     % get handlers

   % get data set
   sets=get(h.axes1,'UserData');

   % get file
   file=get(h.btload,'UserData');

   % zeroize parameters of the separation line
   h.line.t=0;
   h.line.theta=0;
   h.line.alpha=[0;0];
   
   h.stepcnt=0;

   if h.line.handler==-1,
      % create line
      axes(h.axes1);
      
      h.line.handler=...
        line('EraseMode','xor','Color','k','Visible','off','Parent',h.axes1);

      h.badpoint.handler=line('EraseMode','xor','Color','k','Visible','off',...
       'Parent',h.axes1,...
       'Marker','o',...
       'MarkerSize',WRONGPOINT_SIZE);
      
      
      drawnow;
   else
      % change parameters of line
      set(h.line.handler,'Visible','off');
      set(h.badpoint.handler,'Visible','off');
   end % if hline==-1

   % set up handlers and flush queue with graph. objects
   set(hfigure,'UserData',h);

   % create comment
   if isempty(sets)==0,
      consoletext=sprintf('Step t=0\nNo separation line');
      titletext=sprintf('Data file: %s, %d vectors',file.name,sum(sets.K));
   else
      consoletext=sprintf(['No data loaded.\nPress Load data button.\n',...
            'Load sample data from ../toolboxroot/data/binary_separable']);
      titletext='';

      pos=get(h.axes1,'Position');
      fsize=min(pos(3),pos(4))/10;
      setaxis(h.axes1,[-1 1 -1 1]);
      axis([-1 1 -1 1]);
      
      builtin('text',0,0,'Press ''Load data'' button.',...
         'HorizontalAlignment','center',...
         'FontUnits','normalized',...
         'Clipping','on',...
         'FontSize',fsize);
   end

   % show comment
   set(h.console,'String',consoletext );

   % print title
   pos=get(h.axes1,'Position');
   fsize=(1-pos(2)-pos(4))*1;
   title(titletext,...
      'Parent',h.axes1,...
      'VerticalAlignment','bottom',...
      'HorizontalAlignment','left',...
      'FontUnits','normalized',...
      'Units','normalized',...
      'Position',[0 1 0],...
      'FontSize',fsize);


case 'drawline'
   % == Draw separation line ============================

   h=get(hfigure,'UserData');              % get handlers

   % get new line position from input arguments
   theta=varargin{1};
   alpha=varargin{2};

   if mod(h.stepcnt,2)==1 & h.line.t >0, 
     set(h.badpoint.handler,'Visible','on',...
         'XData',alpha(1),'YData',alpha(2));
   else
     set(h.badpoint.handler,'Visible','off');
 
     % Cut off line along axes
     [x1,y1,x2,y2,in]=cliplin1(alpha,theta,getaxis(h.axes1));

     % erase old line
     set(h.line.handler,'Visible','off');

      % draw new line if is in the axes
      if in==1,
         set(h.line.handler,...
           'XData',[x1 x2],...
           'YData',[y1 y2],...
           'LineWidth',LINE_WIDTH,...
           'Visible','on');
      end
   end

    % flush draw queue
   drawnow;
 

case 'info'
   % == Invokes standard Matlab`s info box ==========================
   helpwin(mfilename);
end



%========================================
function [text]=makeinfo(t,alpha,theta,solution)
% assembles text description about current solution state

if solution==1,
   txline{1}=sprintf('Solution was found after t=%d step(s).',t);
else
   txline{1}=sprintf('Step t=%d',t);
end
txline{2}=sprintf('Linear decision function:');
txline{3}=sprintf('%f x_1 + %f x_2 + %f = 0',alpha(1),alpha(2),-theta);

text='';
for i=1:3,
   text=strvcat(text,txline{i});
end

return;


%===========================================
function [alpha,theta,solution,tplus1]=exec(hfigure);

h=get(hfigure,'UserData');                           

if h.line.handler==-1,
   return;
end

% get sets
sets=get(h.axes1,'UserData');

% no data set loaded
if isempty(sets)==1,
   return;
end

% get parameters
t=h.line.t;
alpha=h.line.alpha;
theta=h.line.theta;

iter=max(1,round(str2num(get(h.editer,'String'))));
epsil=str2num(get(h.edeps,'String'));

if mod(h.stepcnt,2) ~= 0 & t > 0,
  iter = -1;
end

% perform one adaptation step
switch get(h.pualgo,'Value')
 case 1
%
%[alpha,theta,solution,tplus1]=perceptr(sets.X,sets.I,iter,t,alpha,theta);
  if iter==-1, tmp_options.tmax = -1; else tmp_options.tmax=t+iter; end
   init_model.t=t;
   init_model.W=alpha;
   init_model.b=-theta;
   tmp_model = perceptron(sets,tmp_options,init_model);
   
   if iter==-1,
     if tmp_model.last_update, alpha=sets.X(:,tmp_model.last_update); 
     else alpha=tmp_model.W; end
   else
     alpha = tmp_model.W;
   end
   solution=tmp_model.exitflag;
   tplus1 = tmp_model.t;
   theta = -tmp_model.b;

 case 2
%%
%%  [alpha,theta,solution,tplus1]=kozinec(sets.X,sets.I,iter,t,alpha,theta);
  if iter==-1, tmp_options.tmax = -1; else tmp_options.tmax=t+iter; end
   tmp_options.eps=-1;
   if t~=0, 
    init_model.t=t;
    init_model.W=alpha;
    init_model.b=-theta;
    tmp_model = ekozinec(sets,tmp_options, init_model);
   else
     tmp_model = ekozinec(sets,tmp_options);
   end
   
   if iter==-1,
     if tmp_model.last_update, alpha=sets.X(:,tmp_model.last_update); 
     else alpha=tmp_model.W; end
   else
     alpha = tmp_model.W;
   end
   solution=tmp_model.exitflag;
   tplus1 = tmp_model.t;
   theta = -tmp_model.b;

 case 3
%    [alpha,theta,solution,tplus1]=ekozinec(sets.X,sets.I,epsil,iter,t,...
%                                           alpha,theta);      
  if iter==-1, tmp_options.tmax = -1; else tmp_options.tmax=t+iter; end
   if t~=0, 
    init_model.t=t;
    init_model.W=alpha;
    init_model.b=-theta;
    tmp_options.eps=epsil;
    tmp_model = ekozinec(sets,tmp_options,init_model);
   else
    tmp_options.eps=epsil;
    tmp_model = ekozinec(sets,tmp_options);
   end
   
   if iter==-1,
     if tmp_model.last_update, alpha=sets.X(:,tmp_model.last_update); 
     else alpha=tmp_model.W; end
   else
     alpha = tmp_model.W;
   end
   solution=tmp_model.exitflag;
   tplus1 = tmp_model.t;
   theta = -tmp_model.b;

 case 4
    tmp_options.ker='linear';
    tmp_options.C = 1e12;
%    [alpha,theta,solution]=svmmot(sets.X,sets.I);
    tmp_model=smo(sets,tmp_options);
%    if tmp_model.exitflag <=0, solution=-1; else solution=1; end
    if tmp_model.trnerr > 0, solution = -1; else solution=1; end
    theta=-tmp_model.b;
    alpha=tmp_model.sv.X*tmp_model.Alpha(:);
    tplus1=1;
 case 5
%    [alpha,theta,solution,tplus1]=psum(sets.X,sets.I,10,iter,t,alpha,theta);
% case 6
%    [alpha,theta,solution,tplus1]=psumv(sets.X,sets.I,iter,t,alpha,theta);
% case 7
%    [alpha,theta,solution,tplus1]=navara1(sets.X,sets.I,iter,t,alpha,theta);
% case 8
%    [alpha,theta,solution,tplus1]=navarah1(sets.X,sets.I,iter,t,alpha,theta);
% case 9
%    [alpha,theta,solution]=simplex(sets.X,sets.I);
%    if solution==0, solution=-1; end
%    tplus1=1;

end

return


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

%========================================                                       
function [rect]=getaxis(handle)                                                 
% function [rect]=getaxis(handle)                                               
%                                                                               
% GETAXIS returns a row vector containing the scaling for                       
%   the plot with a given handle.                                               
%                                                                               
% See also AXIS.                                                                
%                                                                               
% Statistical Pattern Recognition Toolbox, Vojtech Franc, Vaclav Hlavac         
% (c) Czech Technical University Prague, http://cmp.felk.cvut.cz                
% Written Vojtech Franc (diploma thesis) 02.01.2000                             
% Modifications                                                                 
% 24. 6.00 V. Hlavac, comments polished.                                        
                                                                                
rect=[get(handle,'XLim'),get(handle,'YLim'),get(handle,'ZLim')];                
                                                                                
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


function [x1,y1,x2,y2,in]=cliplin1(alpha,theta,window)                          
% [x1,y1,x2,y2,in]=cliplin1(alpha,theta,window)                                 
%                                                                               
% CLIPLIN1 clips the line given by the equation alpha*x=theta along             
%   the window. It returns two points on the border of the window.              
%   If the line is in the window then the argument is equal to 1                
%   else it returns 0.                                                          
%                                                                               
% See also CLIPLIN2.                                                            
%                                                                               
% Statistical Pattern Recognition Toolbox, Vojtech Franc, Vaclav Hlavac         
% (c) Czech Technical University Prague, http://cmp.felk.cvut.cz                
% Written Vojtech Franc (diploma thesis) 20.10.1999, 23.12.1999                 
% Modifications                                                                 
% 24. 6.00 V. Hlavac, comments polished.                                        
                                                                                
minx=window(1);                                                                 
maxx=window(2);                                                                 
miny=window(3);                                                                 
maxy=window(4);                                                                 
                                                                                
x=zeros(4,1);                                                                   
y=zeros(4,1);    
if alpha(1)==0,                                                                 
   if alpha(2)~=0,                                                              
      x1=minx;                                                                  
      y1=theta/alpha(2);                                                        
      x2=maxx;                                                                  
      y2=y1;                                                                    
      in=1;                                                                     
   else                                                                         
      % if alpha == 0 then it means the bad input.                              
      x1=0;                                                                     
      y1=0;                                                                     
      x2=0;                                                                     
      y2=0;                                                                     
      in=0;                                                                     
   end                                                                          
elseif alpha(2)==0,                                                             
   x1=theta/alpha(1);                                                           
   y1=miny;                                                                     
   x2=x1;                                                                       
   y2=maxy;                                                                     
   in=1;                                                                        
 else      
     y(1)=maxy;                                                                   
   x(1)=(theta-alpha(2)*y(1))/alpha(1);                                         
   y(2)=miny;                                                                   
   x(2)=(theta-alpha(2)*y(2))/alpha(1);                                         
                                                                                
   x(3)=maxx;                                                                   
   y(3)=(theta-alpha(1)*x(3))/alpha(2);                                         
   x(4)=minx;                                                                   
   y(4)=(theta-alpha(1)*x(4))/alpha(2);                                         
                                                                                
   j=0;                                                                         
   for i=1:4,                                                                   
      if x(i) <= maxx & x(i) >= minx & y(i) <= maxy & y(i) >= miny,             
         if j==0,                                                               
            j=j+1;                                                              
            x1=x(i);                                                            
            y1=y(i);                                                            
         elseif j==1,                                                           
            j=j+1;                                                              
            x2=x(i);                                                            
            y2=y(i);                                                            
         end                                                                    
      end                                                                       
   end                                                                          
                                                                                
   if j<2,                                                                      
      x1=0;                                                                     
      y1=0;                                                                     
      x2=0;                                                                     
      y2=0;                                                                     
      in=0;                                                                     
   else                                                                         
      in=1;                                                                     
   end                                                                          
end % elseif alpha(2)==0                                                        
     

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

function model=ekozinec(data,options,init_model)
% EKOZINEC Kozinec's algorithm for eps-optimal separating hyperplane.
%
% Synopsis:
%  model = ekozinec(data)
%  model = ekozinec(data,options)
%  model = ekozinec(data,options,init_model)
%
% Description:
%  This function is implementation of the Kozinec's algorithm
%  [Koz73] with eps-optimality stopping condition ([SH10], Chapter 5).
%  It transforms the binary-class problem to the single-class problem.
%  The nearest point from convex hull of transformed data which
%  determines the optimal separating hyperplane is found by the 
%  Kozinec's algorithm.
%  
%  model = ekozinec(data) finds eps-optimal separating hyperplane for
%   linearly separable data.
%
%  model = ekozinec(data,options) specifies stopping conditions of
%   the algorithm in structure options:
%    .eps [1x1] ... controls how close is the found solution to
%        the optimal hyperplane in terms of Euclidien norm of the
%        normal vector (default 0.01).
%    .tmax [1x1]... maximal number of iterations.
%  
%   If tmax==-1 then it only returns index (model.last_update)
%   of data vector which should be used by the algorithm for updating
%   the linear rule in the next iteration.
%  
%  model = ekozinec(data,options,init_model) specifies initial model
%   which must contain:
%    .W [dim x 1] ... normal vector.
%    .b [1x1] ... bias of hyperplane.
%    .t [1x1] (optional) ... iteration number.
%
% Input:
%  data [struct] Labeled (binary) training data. 
%   .X [dim x num_data] Input vectors.
%   .y [1 x num_data] Labels (1 or 2).
%
%  options [struct] 
%   .eps [1x1] Precision of found hyperplane (see above).
%   .tmax [1x1] Maximal number of iterations (default tmax=inf).
%     If tmax==-1 then it does not perform any iteration but returns only 
%     index of the point which should be used to update linear rule.
%  
%  init_model [struct] Initial model; must contain items
%    .W, .b and .t (see below).
%
% Output:
%  model [struct] Found linear classifier:
%   .W [dim x 1] Normal vector of hyperplane.
%   .b [1x1] Bias of hyperplane.
%  
%   .exitflag [1x1] 1 ... eps-optimality condition satisfied
%                   0 ... number of iterations exceeded tmax.
%   .t [int] Number of iterations.
%   .last_update [d x 1] Index of the last point used for update.
%
% Example:
%  data = genlsdata( 2, 50, 1);
%  model = ekozinec(data)
%  figure; ppatterns(data); pline(model); 
%
% See also 
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 3-may-2004, VF
% 17-sep-2003, VF
% 16-Feb-2003, VF
% 21-apr-2001, V.Franc, created

% get data dimensions
[dim,num_data] = size(data.X);

% Process input arguments 
% ------------------------------------
if nargin < 2, options = []; else options = c2s( options ); end

if ~isfield(options,'tmax'), options.tmax = inf; end
if ~isfield(options,'eps'), options.eps = 0.01; end

if nargin < 3, 
  % make inital model
  model.b = 0;
  if data.y(1)==1, W = [data.X(:,1);1]; else W = -[data.X(:,1);1]; end
else
  % get inital model from input
  model = init_model;
  W = [model.W; model.b];
end
if ~isfield( model,'t'), model.t = 0; end

model.exitflag = 0;
model.last_update = 0;

% Add one constant coordinates to the data and swap                             
% points from the second class along the origin.
% ---------------------------------------------------- 
data.X = [data.X; ones(1,num_data)];
dim=dim+1;
inx = find(data.y==2);
data.X(:,inx) = -data.X(:,inx);

if options.tmax == -1,
  % return index of point which should be used to update linear rule
  %----------------------------------------------------------------------
  norm_W = norm(W);
  [min_proj,min_inx]=min( data.X'*W/norm_W );

  % bound for separating or eps-optimal separating hyperplane 
  if options.eps < 0, bound = norm_W/2; else bound = norm_W - options.eps; end

  model.last_update = min_inx;
  if min_proj <= bound, model.exitflag = 0; else model.exitflag = 1; end
else 

  % main loop
  %----------------
  while model.exitflag == 0 & options.tmax > model.t,

    model.t = model.t + 1;

    norm_W = norm(W);
    [min_proj,min_inx]=min( data.X'*W/norm_W );

    % bound for separating or eps-optimal separating hyperplane 
    if options.eps < 0, bound = 0; else bound=norm_W - options.eps; end
    
    if min_proj <= bound,
      model.last_update = min_inx;
      xt=data.X(:,min_inx);
    
      k=min(1,(W-xt)'*W/((W-xt)'*(W-xt)) );
      W = W*(1-k) + xt*k;
      
      model.exitflag = 0;
    else
      model.exitflag = 1;
    end
  end
  
end

model.W = W(1:dim-1);
model.b = W(dim);
model.fun = 'linclass';

return;

