function result = demo_anderson(action,hfigure,varargin)
% DEMO_ANDERSON Demo on Generalized Anderson's task.
%
% Synopsis:
%  demo_anderson
%
% Description:
%  This demo demonstrates the algorithms which solve 
%  the Generalized Anderson`s Task (GAT) [SH10]. The GAT is an 
%  instance of the non-Bayesian task of decision under 
%  non-random intervention. 
% 
%  The goal of is to find a binary linear classification
%  rule (g(x)=sgn(W'*x+b) (line in 2D) with minimal probability of
%  misclassification. The conditional probabilities are known to
%  be Gaussians their paramaters belong to a given set of 
%  parameters. The true parameters are not known. The linear rule 
%  which guarantes the minimimal classification error for the worst
%  possible case (the worst configuration of Gaussains) is
%  sought for.
%  
%  The found solution (hyperplane, line in 2D) is vizualized 
%  as well as the input Gaussians which describe input classes.
%
%  Following algorithms can be tested:
%     
%  Eps-solution - Finds epsilon-solution of the GAT in finite number
%              of iterations if such solution exist. The epsilon means
%              desired classification error.
%  Original  - Original Anderson-Bahadur's algorithm defined for 
%              two Gaussians only (each class one Gaussian).
%  Optimal   - Implementation of general algorithm propsed by Schlesinger.
%              It finds the optimal solution.
%  Gradient  - Fast and simple implementation which uses the generalized
%              gradient descent optimization.
%
% Control:
%  Algorithm  - select algorithm for testing.
%  Parameter  - parameters for the selected algorithm.
%  Iterations - number of iterations in one step.
%  Animation  - enable/dissable animation.
%
%  FIG2EPS     - export screen to the PostScript file.
%  Load data   - load input point sets from file.
%  Create data - call interactive program for creating sets of Gaussians.
%  Reset       - set the tested algorithm to the initial state.
%  Play        - run the tested algorithm.
%  Stop        - stop the running algorithm.
%  Step        - perform only one step.
%  Info        - display the info box.
%  Close       - close the program.
%
% See also 
%  EANDERS, ANDRORIG, GGRADANDR, GANDERS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 17-sep-2003, VF
% 11-June-2001, V.Franc, comments added.
% 24. 6.00 V. Hlavac, comments polished.

% constants
ALGONAMES=['Eps-solution ';...
           'Original     ';...
           'Gradient     ';...
           'Optimal      '];
PREC_TITLE=['Max error (0-50) [%]';...    % e-Optimal solution
            'd(lambda,ni) (0,inf)';...    % Original Anderson`s solution
            'd( min r ) (0,inf)  ';...    % Gradient method
            'd( min r ) (0,inf)  '];      % General solution 
DEF_PRECISION=[10,1e-3,1e-3,1e-3];  % default values of the precision of algo. 1,2,...
BORDER=0.5;
DATA_IDENT='Infinite sets, Normal distributions';       % M-file identifier
PLOT_FCE='pandr2d';     % outlined ellipsoids
%PLOT_FCE='pandr2df';     % outlined ellipsoids

% if number of arguments is less then 1, that means first call of this function. Every
% other calls set up at least argument action
if nargin < 1,
   action = 'initialize';
end

% what action is required ?
switch lower(action)

case 'initialize'
   % == Initialize user interface control and figure window ================

   % == Figure =============================================================
   left=0.1;
   width=0.8;
   bottom=0.1;
   height=0.8;
   hfigure=figure('Name','Anderson`s task', ...
      'Visible','off',...
      'Units','normalized', ...
       'NumberTitle','off', ...
      'Position',[left bottom width height],...
      'tag','Demo_Anderson',...
      'Units','normalized', ...
      'RendererMode','manual');

   % == Axes ===============================================================
   left=0.1;
   width=0.65;

   % axes for showing sets
   bottom=0.35;
   height=0.60;
   haxes1=axes(...
       'Units','normalized', ...
      'Box','on', ...
      'DrawMode','fast',...
      'NextPlot','add',...
      'Layer','top',...
      'UserData',[],...
      'Position',[left bottom width height]);
   xlabel('feature x');
   ylabel('feature y');

   % == Comment window =====================================================

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
      'Callback','demo_anderson(''info'',gcf)',...
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
      'Callback','demo_anderson(''step'',gcf)');

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
      'Callback','demo_anderson(''play'',gcf)');

   % Reset button: set up t = 0
   bottom=bottom+height;
    hbtreset = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Reset', ...
      'Callback','demo_anderson(''reset'',gcf)');

   % Create data
   bottom=bottom+1.5*height;
    hbtcreat = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Create data', ...
      'Callback','demo_anderson(''creatdata'',gcf)');

   % Load data
   bottom=bottom+1*height;
    hbtload = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Load data', ...
      'Callback','demo_anderson(''getfile'',gcf)');

   % == Check boxes ===============================================================

   % Make chack box to determine if a line will be drawn in one step or smoothly.
   bottom=bottom+height*1.2;
    hxbanim = uicontrol(...
    'Style','checkbox', ...
       'Units','normalized', ...
    'ListboxTop',0, ...
       'Position',[left bottom width height], ...
    'String','Animation');


   % == Popup menus ==========================================================

   % Pop up menu for the selection between algorithms
    % title
   bottom=0.95-height;
   htxalgo=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'String','Algorithm');
   % popup menu
   bottom=bottom-0.9*height;
   hpualgo=uicontrol( ...
      'Style','popup', ...
      'Units','normalized', ...
      'CallBack','demo_anderson(''algohandler'',gcf)',...
      'Position',[left bottom width height], ...
      'UserData',1,...
      'String',ALGONAMES);


   % == Edit lines ================================================================

   bottom=0.95-3.5*height;
   % Precision of solution
   htxprec=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'String',PREC_TITLE(1,:));

   bottom=bottom-height;
   hedprec = uicontrol(...
    'Units','normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'Style','edit',...
      'String',num2str(DEF_PRECISION(1)) );

   % # of iterations
   bottom=bottom-1.5*height;
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
      'String',1);

    % ==============================================================================

   % Store handlers
   handlers=struct(...
      'line',struct('handler',-1,'alpha',0,'alpha1',0,'alpha2',0,'lambda',0,'theta',0,'t',0),...
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
      'edprec',hedprec,...
      'editer',hediter,...
      'txprec',htxprec,...
      'axes1',haxes1,...
      'xbanim',hxbanim);

   set(hfigure,'UserData',handlers);

   % Reset adaptation, t=0
   demo_anderson('reset',hfigure);

   % Put figure on desktop
   set(hfigure,'Visible','on');

   drawnow;


case 'creatdata'
   % == Invoke data set creator ============================================
%   creatset('normal',2,'demo_anderson','created',hfigure);

   createdata('gauss',2,'demo_anderson','created',hfigure);

case 'created'
   % == Load new created data set ===========================================

   % get handler and make this figure active
   figure(hfigure);
   h=get(hfigure,'UserData');

   % get file name
   path=varargin{1};
   name=varargin{2};
   pathname=strcat(path,name);

%   if checkdat(pathname,DATA_IDENT,2,[0 0])==1,
  if check2dgauss(pathname),
      file.pathname=pathname;
      file.path=path;
      file.name=name;
      set(h.btload,'UserData',file);
      demo_anderson('loadsets',hfigure);
      demo_anderson('reset',hfigure);
   else
      errordlg('This file does not contain required data.','Bad file','modal');
   end


case 'getfile'
   % == Invoke standard open file dialog ====================================
   % Opens file and checks if contains apropriate data, if yes loads data.

   h=get(hfigure,'UserData');

   % change path to directory
%%   wres=what('anderson');
%%   cd(wres.path);

   [name,path]=uigetfile('*.mat','Open file');
   if name~=0,
      file.pathname=strcat(path,name);
      file.path=path;
      file.name=name;
%      if checkdat(file.pathname,DATA_IDENT,2,[0 0])==1,
      if check2dgauss(file.pathname),
         set(h.btload,'UserData',file);
         demo_anderson('loadsets',hfigure);
         demo_anderson('reset',hfigure);
      else
         errordlg('This file does not contain required data.','Bad file','modal');
      end
   end


case 'algohandler'
   % == Handler for check box 'Algorithm' =======================================
   % If new algorithm is selected then prepare data for it.

   h=get(hfigure,'UserData');

   if get(h.pualgo,'UserData') ~= get(h.pualgo,'Value'),
      set(h.pualgo,'UserData',get(h.pualgo,'Value'));

      set(h.edprec,'String',num2str(DEF_PRECISION(get(h.pualgo,'Value'))) );
      set(h.txprec,'String',PREC_TITLE(get(h.pualgo,'Value'),:));

      demo_anderson('loadsets',hfigure);
      demo_anderson('reset',hfigure);
   end


case 'loadsets'
   % == Load sets ==============================================
   % Get given file name and load the data set from him.

   h=get(hfigure,'UserData');                   % uicontrol handlers

   % Get file name
   file=get(h.btload,'UserData');
   if isempty(file)==0,

      % Load sets
      sets=load(file.pathname);

      sets.MI = sets.Mean;
      sets.SIGMA = reshape(sets.Cov,2,2*size(sets.Mean,2));
      sets.I = sets.y;
      sets.K = zeros(1,max(sets.y));
      sets.N = size(sets.Mean,1);
      for i=1:max(sets.y),
        sets.K(i) = length(find(sets.y==i));
      end
      
      % algorithm 2 (Original Anderson`s) solution find solution 
      % for two distributions only
      if get(h.pualgo,'Value')==2,
         % get only one distribution from each class
         class1=0;
         class2=0;
         NI=[];
         NMI=[];
         NSIGMA=[];
         NK=[1 1];
         N = 2;
         i=0;
         while i<sets.K | class1==0 |class2==0,
            i=i+1;
            if sets.I(i)==1 & class1==0,
               class1=1;
               NI=[NI,1];
               NMI=[NMI,sets.MI(:,i)];
               NSIGMA=[NSIGMA,sets.SIGMA(:,(i-1)*2+1:i*2)];
            elseif sets.I(i)==2 & class2==0,
               class2=1;
               NI=[NI,2];
               NMI=[NMI,sets.MI(:,i)];
               NSIGMA=[NSIGMA,sets.SIGMA(:,(i-1)*2+1:i*2)];
            end
         end % while

         % replace old values
         sets.MI=NMI;
         sets.SIGMA=NSIGMA;
         sets.I=NI;
         sets.K=NK;
      end % if get(h.pualgo,....
   else
      % No set is loaded.
      sets=[];
   end

   % store sets
   set(h.axes1,'UserData',sets);


case 'play'
   % == Start up the adaptation process =======================================
   % Perform the adaptation step by step until the solution is found or stop
   % button is pushed down.

   h=get(hfigure,'UserData');                      % get handlers

   % get sets
   sets=get(h.axes1,'UserData');

   % if some set is loaded then perform on step
   if isempty(sets)==0,

      % Disable buttons everything axcept
      set([h.btinfo h.btstep h.btclose h.btplay h.btreset h.btload h.pualgo ...
            h.btcreat h.editer h.edprec h.xbanim],'Enable','off');

      % Only stop button can be pushed down
      set(h.btstop,'Enable','on');

      % Stop button was not pushed down
      set(h.btstop,'UserData',0);
      play=1;

      % get arguments from dialog
      anim=get(h.xbanim,'Value');

      % Play - adaptation process
      while play==1 & get(h.btstop,'UserData')==0,

         % get arguments from dialog
         oldtheta=h.line.theta;
         oldalpha=h.line.alpha;

         % call algorithm
         [h,text,play,solution]=callalgo(h,sets);

         % appear result
         set(h.console,'String',text );
         drawnow;

         if play~=-1,
            % plot result
            if h.line.handler==-1,
               axes(h.axes1);
               h.line.handler=feval(PLOT_FCE,sets.MI,sets.SIGMA,sets.I,...
                  h.line.alpha,h.line.theta );
            else
               feval(PLOT_FCE,sets.MI,sets.SIGMA,sets.I,h.line.alpha,h.line.theta,...
                  h.line.handler,anim,oldalpha,oldtheta);
            end
         end

         % hands on control to MATLAB
         drawnow;
      end % while play == 1 & get(...

      %  store new solution
      set(hfigure,'UserData',h);

      % enable these buttons
      set([h.btinfo h.btstep h.btclose h.btplay h.btreset h.btload ...
         h.btcreat h.pualgo h.editer h.edprec h.xbanim],'Enable','on');

      % disable stop button
      set(h.btstop,'Enable','off');

   else % isempty(sets)==0,
      % write down description
      text=sprintf(['No data loaded.\nPress Load data or Create data button.\n', ...
         'Load sample data from ../toolboxroot/data/anderson_task/' ]  );
      set(h.console,'String',text );
   end   


case 'step'
   % == Perform one adaptation step ================================================
   h=get(hfigure,'UserData');                     % get handlers we will need...

   % get sets
   sets=get(h.axes1,'UserData');

   % if some set is loaded then perform on step
   if isempty(sets)==0,

      % get arguments from dialog
      anim=get(h.xbanim,'Value');
      oldtheta=h.line.theta;
      oldalpha=h.line.alpha;

      % call algorithm
      [h,text,play,solution]=callalgo(h,sets);

      % appear result
      set(h.console,'String',text );
      drawnow;

      if play~=-1,
         % plot result
         if h.line.handler==-1,
            axes(h.axes1);
            h.line.handler=feval(PLOT_FCE,sets.MI,sets.SIGMA,sets.I,...
               h.line.alpha,h.line.theta );
         else
            feval(PLOT_FCE,sets.MI,sets.SIGMA,sets.I,h.line.alpha,...
               h.line.theta,h.line.handler,anim,oldalpha,oldtheta);
         end
      end

      drawnow;

      % store data
      set(hfigure,'UserData',h);

   else % isempty(sets)==0,
      % write down description
      text=sprintf(['No data loaded.\nPress Load data or Create data button.\n', ...
         'Load sample data from ../toolboxroot/data/anderson_task/' ]  );
      set(h.console,'String',text );
   end

case 'reset'
   % == Reset adaptation process ==================================
   % Sets t=0 and redraws axes.

   h=get(hfigure,'UserData');                     % get handlers

   % zeroize all parameters of the solution
   h.line.t=0;
   h.line.theta=0;
   h.line.alpha=[0;0];
   h.line.alpha1=[0;0];
   h.line.alpha2=[0;0];
   h.line.lambda=0;
   set(hfigure,'UserData',h);

   % No line
   h.line.handler=-1;
   set(hfigure,'UserData',h);
   %%%   cla;
   clrchild(h.axes1);
   %%%   win=axis;
   win=getaxis(h.axes1);
   %%%   axis([0 1 0 1]);
   setaxis(h.axes1,[0 1 0 1]);
   %%%   axis(win);
   setaxis(h.axes1,win);

   % Redraw points
   sets=get(h.axes1,'UserData');

   % if some points are loaded than appear it
   if isempty(sets)==0,

      % set axes according to current points MI

      if sum(sets.K) < 3,
         win=cmpwin(min(sets.MI'),max(sets.MI'),BORDER*2,BORDER*2);
      else
         win=cmpwin(min(sets.MI'),max(sets.MI'),BORDER,BORDER);
      end

      %%% axis(win);
      setaxis(h.axes1,win);

%%      pplot(sets.MI,sets.I);
%%      ppoints(sets.MI,sets.I);
      ppatterns(sets.MI,sets.I);

      % comment window text
      consoletext=sprintf('Press Step button or Play button.\n');

      file=get(h.btload,'UserData');
      titletext=sprintf('File: %s, # of distributions K = %d',file.name,sum(sets.K));

   else % if isempty(sets)==0,
      % comment window text
%      consoletext=sprintf(['No data loaded.\n' ...
%             'Press Load button or Create data button.\n']);
%      consoletext=sprintf('No data loaded.\nPress Load data button.\n');
      consoletext=sprintf(['No data loaded.\nPress Load data or Create data button.\n', ...
         'Load sample data from ../toolboxroot/data/anderson_task/' ]  );
      titletext='';

      pos=get(h.axes1,'Position');
      fsize=min(pos(3),pos(4))/10;

      %%%      axis([-1 1 -1 1]);
      setaxis(h.axes1,[-1 1 -1 1]);
      builtin('text',0,0,'Press ''Load data'' button.',...
         'HorizontalAlignment','center',...
         'FontUnits','normalized',...
         'Clipping','on',...
         'FontSize',fsize);
   end

   % print comment
   set(h.console,'String',consoletext );

   % print title
   pos=get(h.axes1,'Position');
   fsize=(1-pos(2)-pos(4))*1;
   title(titletext,...
      'VerticalAlignment','bottom',...
      'HorizontalAlignment','left',...
      'FontUnits','normalized',...
      'Units','normalized',...
      'Position',[0 1 0],...
      'FontSize',fsize);


case 'info'
   % == Call standard Matlab info box =========================================
   helpwin(mfilename);


end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [h,text,play,solution]=callalgo(h,sets)

% get arguments from dialog
precision=str2num(get(h.edprec,'String'));
iter=max([1,fix(str2num(get(h.editer,'String')))]);

% get parameters
t=h.line.t;
alpha=h.line.alpha;
alpha1=h.line.alpha1;
alpha2=h.line.alpha2;
lambda=h.line.lambda;
theta=h.line.theta;

% perform algorithm
switch get(h.pualgo,'Value')
case 4
   % General solution

   distr = sets;
   options.eps = precision;
   options.tmax = t+iter;
   init_model.W=alpha;
   init_model.b=-theta;
   init_model.t = t;
  
%   [nalpha,ntheta,solution,minr,nt,maxerr]=...
%      ganders(sets.MI,sets.SIGMA,sets.I,iter,precision,t,alpha,theta);
   model = ganders( distr, options, init_model);

   nalpha = model.W;
   ntheta = -model.b;
   solution = model.exitflag;
   minr = model.r;
   nt = model.t;
   maxerr = model.err;
   
  text=sprintf(['Interation(s) t=%d\nLinear rule q(x)=sgn([%f , %f]*x %+f)\n'...
                 'Minimal r = %.8f\nClassification error = %.4f%%'],...
      nt,nalpha(1),nalpha(2),-ntheta,minr,maxerr*100);
case 1

   distr = sets;
   options.err = precision/100;
   options.tmax = t+iter;
   init_model.W1=alpha1;
   init_model.W2=alpha2;
   init_model.t = t;

   % e-Optimal solution
   model = eanders( distr, options, init_model );

   nalpha = model.W;
   ntheta = -model.b;
   alpha1 = model.W1;
   alpha2 = model.W2;
   solution = model.exitflag;
   minr = model.r;
   nt = model.t;
   maxerr = model.err;
  
%   [nalpha,ntheta,solution,nt,alpha1,alpha2]=...
%       eanders(sets.MI,sets.SIGMA,sets.I,iter,precision/100,t,alpha1,alpha2);
%    if sum(nalpha)==0,
%      solution=-1;
%      nalpha=alpha;
%      ntheta=theta;
%   end
% text=sprintf('Step t=%d\nLine [%f , %f]*x=%f',nt,nalpha(1),nalpha(2),ntheta);
  text=sprintf(['Interation(s) t=%d\nLinear rule q(x)=sgn([%f , %f]*x %+f)\n'...
                 'Minimal r = %.8f\nClassification error = %.4f%%'],...
      nt,nalpha(1),nalpha(2),-ntheta,minr,maxerr*100);
case 2
   % Original Anderson`s solution
   distr.Mean = sets.MI;
   distr.Cov = reshape(sets.SIGMA,2,2,2);

   options.eps = precision;
   options.tmax = t+iter;
   init_model.gamma=lambda;
   init_model.t = t;

   model = androrig( distr, options, init_model);

   nalpha = model.W;
   ntheta = -model.b;
   solution = model.exitflag;
   nt = model.t;
   lambda = model.gamma;
   maxerr = model.err;
   minr = min([model.r1 model.r2]);
  
%   [nalpha,ntheta,solution,nt,lambda,ni,maxerr]=...
%      oanders(sets.MI,sets.SIGMA,sets.I,iter,precision,t,lambda);

  text=sprintf(['Interation(s) t=%d\nLinear rule q(x)=sgn([%f , %f]*x %+f)\n'...
                 'Minimal r = %.8f\nClassification error = %.4f%%'],...
      nt,nalpha(1),nalpha(2),-ntheta,minr,maxerr*100);

%  text=sprintf(...
%      'Step t=%d\nLine [%f , %f]*x=%f\nNi = %f, (1-Lambda)/Lambda = %f, Max error = %f%%',...
%      nt,nalpha(1),nalpha(2),ntheta,ni,(1-lambda)/lambda,maxerr*100);
case 3
   distr = sets;
   options.eps = precision;
   options.tmax = t+iter;
   init_model.W=alpha;
   init_model.b=-theta;
   init_model.t = t;

   % e-Optimal solution
   model = ggradandr( distr, options, init_model );

   nalpha = model.W;
   ntheta = -model.b;
   solution = model.exitflag;
   minr = model.r;
   nt = model.t;
   maxerr = model.err;
  
%   [nalpha,ntheta,solution,nt,alpha1,alpha2]=...
%       eanders(sets.MI,sets.SIGMA,sets.I,iter,precision/100,t,alpha1,alpha2);
%    if sum(nalpha)==0,
%      solution=-1;
%      nalpha=alpha;
%      ntheta=theta;
%   end
% text=sprintf('Step t=%d\nLine [%f , %f]*x=%f',nt,nalpha(1),nalpha(2),ntheta);
  text=sprintf(['Interation(s) t=%d\nLinear rule q(x)=sgn([%f , %f]*x %+f)\n'...
                 'Minimal r = %.8f\nClassification error = %.4f%%'],...
      nt,nalpha(1),nalpha(2),-ntheta,minr,maxerr*100);
 
end

if solution==-1,
    text=sprintf('Solution does not exist.\n');
    play=-1;
    return;
elseif solution==1,
   text=strvcat(text,sprintf('Solution was found in %d step(s)',nt));
   play=0;
else
   play=1;
end

% store new values
h.line.t = nt;
h.line.alpha = nalpha;
h.line.alpha1 = alpha1;
h.line.alpha2 = alpha2;
h.line.lambda = lambda;
h.line.theta = ntheta;

return

%==============================================================
function [handler]=pandr2d(MI,SIGMA,I,alpha1,theta1,handler,anim,alpha2,theta2)
% PANDR2D displays solution of Generalized Anderson's task in 2D.
% [handler]=pandr2d(MI,SIGMA,I,alpha1,theta1,handler,anim,alpha2,theta2)
%
% PANDR2D plots given solution of the Generalized Anderson`s task (GAT) in
%  2-dimensional feature space. This function plots separation line 
%  (2D case of the GAT) and input classes which are symbolized by sets 
%  of ellipsoids.
%
%  Input arguments MI, SIGMA and I describe input mixture of normal 
%  distributions. The pair of input arguments {alpha,theta} describes 
%  separation line which is particular solution of the GAT. 
%
%  For information on the meaning of the arguments refer to help of 
%  functions solving the GAT (ganders,ganders2,eanders etc.).
%
%  When the quadruple of input parameters handler, anim, alpha2 and theta2
%  enter function then a change of the solution is depicted too. The pair
%  alpha1, theta1 is old solution and alpha2, theta2 is new solution. The
%  argument handler contains information about graphics abjects used 
%  in the last call of the function and returned in output variable handle. 
%  When the argument anim=1 then the change is animated.
%
% See also PANDR2D, DEMO_ANDERSON, GANDERS, GANDERS2, EANDERS, OANDERS.
%

% Statistical Pattern Recognition Toolbox, Vojtech Franc, Vaclav Hlavac
% (c) Czech Technical University Prague, http://cmp.felk.cvut.cz
% Written Vojtech Franc (diploma thesis) 23.11.1999, 11.5.2000
% Modifications
% 24. 6.00 V. Hlavac, comments polished.
% 2.8.00 V.Franc, comments changed

% constants
BORDER=0.95;
POINT_SIZE=8;
POINT_COLOR='k';
POINT_WIDTH=2;
LINE_WIDTH=1;
LINE_COLOR='k';
ELLIPSE_WIDTH=1;
INTERPOL=50;
ANIM_DIF=20;


% handles input arguments
if nargin < 7,
   anim=0;
end
if nargin < 6,
   handler=-1;
end

%%handler


alpha1=alpha1(:);     % alpha1 will column

% number of ellipses
N=size(MI,2);
% dimension, must be equal to 2
DIM=size(MI,1);

if handler==-1,
   % proportions of a separation line
   window=axis*BORDER;
else
   % proportions of a separation line
   window=getaxis(handler(2*N+3))*BORDER;
end;

if anim==0,

   % computes minimal distance among the line alpha*x=theta and
   % the elipses (x-MI)'*inv(SIGMA)*(x-MI)
    [R,inx]=min( abs(alpha1'*MI-theta1)...
        ./sqrt( reshape(alpha1'*SIGMA,DIM,N)'*alpha1 )' );


   if handler==-1,
      % first painting

      % ellipses
      for i=1:N,
         mi=MI(:,i);
         sigma=SIGMA(:,(i-1)*DIM+1:i*DIM);
%%%         isg=inv(sigma);

%%%         [x,y]=ellipse(isg,INTERPOL,R,mi);
         [x,y]=ellips(mi,inv(sigma),R,INTERPOL);

         handler(i)=plot(x,y,'Color',marker_color(I(i)),...
            'EraseMode','xor',...
            'LineWidth',ELLIPSE_WIDTH,...
            'UserData',sigma);
      end
      % line
      [x1,y1,x2,y2]=cliplin1(alpha1,theta1,window);
      handler(N*2+1)=line([x1 x2],[y1 y2],...
         'LineWidth',LINE_WIDTH,...
         'Color',LINE_COLOR,...
         'EraseMode','xor');

%%      get(handler(N*2+1))

      % pull point
      mi=MI(:,inx);
      sigma=SIGMA(:,(inx-1)*DIM+1:inx*DIM);
      x0=mi-(alpha1'*mi-theta1)*sigma*alpha1/(alpha1'*sigma*alpha1);
      handler(N*2+2)=line(x0(1),x0(2),...
         'LineStyle','none',...
         'Color',POINT_COLOR,...
         'MarkerSize',POINT_SIZE,...
         'LineWidth',POINT_WIDTH,...
         'Marker','x',...
         'EraseMode','xor');

       handler(N*2+3)=gca;

   else % if handler==-1,
      % at least second painting

      % hide all objects
%      set(handler,'Visible','off');

      % ellipses
      for i=1:N,
         mi=MI(:,i);
         sigma=get(handler(i),'UserData');
         [x,y]=ellips(mi,inv(sigma),R,INTERPOL);
         set(handler(i),'XData',x,'YData',y,'Visible','on');
      end

      % line
      [x1,y1,x2,y2]=cliplin1(alpha1,theta1,window);
%%      N*2+1
%      get(handler(N*2+1))
%      handler
%      N
 %     x1
%      x2
      set(handler(N*2+1),'XData',[x1 x2],'YData',[y1 y2],'Visible','on');

      % pull point
      mi=MI(:,inx);
      sigma=SIGMA(:,(inx-1)*2+1:inx*2);
      x0=mi-(alpha1'*mi-theta1)*sigma*alpha1/(alpha1'*sigma*alpha1);
      set(handler(N*2+2),'XData',x0(1),'YData',x0(2),'Visible','on');

   end % if handler==-1,

else % if anim==0,

   alpha2=alpha2(:);   % alpha2 will column

   % computes number of the animation steps
   [x11,y11,x12,y12]=cliplin1(alpha1,theta1,window);
   [x21,y21,x22,y22]=cliplin1(alpha2,theta2,window);

   difa=max([sqrt((x11-x21)^2+(y11-y21)^2),sqrt((x11-x21)^2+(y11-y21)^2)]);
   difw=min([window(4)-window(3), window(2)-window(1)]);
   nsteps=max([1,ceil(ANIM_DIF*(difa/difw))]);

   % play
   for i=1:nsteps,
      k=i/nsteps;
      alpha=(1-k)*alpha2+k*alpha1;  % smooth transition of alpha
      theta=(1-k)*theta2+k*theta1;  % --//-- theta

      % recursion
      handler=pandr2d(MI,SIGMA,I,alpha,theta,handler,0);
   end

end % if anim==0,


%=========================================================
function [x1,y1,x2,y2,in]=cliplin1(alpha,theta,window)
% [x1,y1,x2,y2,in]=cliplin1(alpha,theta,window)
%
% CLIPLIN1 clips the line given by the equation alpha*x=theta along
%   the window. It returns two points on the border of the window.
%   If the line is in the window then the argument is equal to 1
%   else it returns 0.
%

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

delete(get(handle,'Children'));

return;

function [win]=cmpwin(mins,maxs,xborder,yborder)

dx=max( (maxs(1)-mins(1)), 1 )*xborder;
dy=max( (maxs(2)-mins(2)), 1 )*yborder;
x1=(mins(1)-dx);
x2=(maxs(1)+dx);
y1=(mins(2)-dy);
y2=(maxs(2)+dx);
win=[x1 x2 y1 y2];

%========================================
function [rect]=getaxis(handle)

rect=[get(handle,'XLim'),get(handle,'YLim'),get(handle,'ZLim')];

return;

function []=setaxis(handle,rect)

set(handle,'XLim',rect(1:2));
set(handle,'YLim',rect(3:4));

if size(rect,2)>=6,
   set(handle,'ZLim',rect(5:6));
end

return;
