function demo_emgmm(action,hfigure,varargin)
% DEMO_EMGMM Demo on Expectation-Maximization (EM) algorithm.
%
% Synopsis:
%  demo_emgmm
%
% Description:
%  This demo shows the Expectation-Maximization (EM) algorithm
%  [Schles68][DLR77] for Gaussians mixture model (GMM). The EM 
%  fits the GMM to i.i.d. sample data (in this case only 2D) 
%  such that the likelihood is maximized. 
%
%  The found model is described by ellipsoids (shape of 
%  covariances) and a crosses (mean value vectors). The value
%  of the optimized log-likelihood function for the current estimate 
%  is displayed in the bottom part.
%
% Control:
%  Covariance  - Determines type of the covariance matrix:
%                Diagonal (independent features),
%                Full (correlated features).
%  Components  - Number of components (Gaussians) in the mixture.
%               
%  Iterations  - Number of iterations in one step.
%  Random init - the initial model is randomly generated and/or 
%                first n training samples are taken as the
%                mean vectors.
%
%  FIG2EPS     - Export screen to the PostScript file.
%  Save model  - Save current model to file.
%  Load data   - Load input point sets from file.
%  Create data - Invoke program for creating point sets.
%  Reset       - Set the tested algorithm to the initial state.
%  Play        - Run the tested algorithm.
%  Stop        - Stop the running algorithm.
%  Step        - Perform only one step.
%  Info        - Info box.
%  Close       - Close the program.
%
% See also EMGMM.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 19-sep-2003, VF
% 11-june-2001, V.Franc, comments added.
% 27.02.00 V. Franc
%  5. 4.00 V. Franc
% 23.06.00 V. Hlavac Comments polished. Message when no data loaded.
%                    Export of the solution to global variables.
% 27-mar-2001, V.Franc, Graph og log-likelihood function added


% Used functions: PPOINTS, PNMIX

% == Global variables, used to export results from demo_emgmm ========

global UNSU_MI       % K vectors with mean values
global UNSU_SIGMA    % K covariance matrices
global UNSU_PK       % K apriori probabilities for each distributions.
%global UNSU_eI      % Used only by the next iteration, not globalised
global UNSU_solution % 1 if alg. finished in stationary point, 0 otherwise
global UNSU_t        % Number of iterations the algorithm performed

% == Constants =====================================================

AXIST_ADD=10;
AXISY_ADD=5;

BORDER=0.25;           % space betwean window outer and the points
CENTERSIZE=10;         % size of center point
LINE_WIDTH=1;
AXIST_ADD=10;
DATA_IDENT='Finite sets, enumeration';   % file identifier
randinit=1;

if nargin < 1,
   action = 'initialize';
end

% What action is required ?
switch lower(action)

case 'initialize'
   % == Initialize user interface control and figure window ================

   % == Figure
   % =============================================================
   left=0.2;
   width=0.6;
   bottom=0.1;
   height=0.8;
   hfigure=figure('Name','EM algorithm', ...
      'Visible','off',...
      'NumberTitle','off', ...
      'Units','normalized', ...
      'Position',[left bottom width height],...
      'tag','demo_emgmm',...
      'doublebuffer','on',...
      'backingstore','off');

   
   % == Axes ===============================================================
   left=0.1;
   width=0.65;
   bottom=0.45;
   height=0.5;
   haxes1=axes(...
      'Units','normalized', ...
      'NextPlot','add',...
      'UserData',[],...
      'Position',[left bottom width height]);
   xlabel('feature x\_1');
   ylabel('feature x\_2');

   htitle1=title('No data loaded',...
      'VerticalAlignment','bottom',...
      'Parent',haxes1,...
      'HorizontalAlignment','left',...
      'Units','normalized',...
      'Position',[0 1 0]);

   % axes log-Likelihood graph
   left=0.1;
   width=0.65;
   bottom=0.1;
   height=0.25;
   haxes2=axes(...
      'Units','normalized', ...
      'NextPlot','add',...
      'Position',[left bottom width height]);
   ylabel('logL(t)');
   
   htitle2=title('Log-likelihood function',...
       'Parent',haxes2,...
       'VerticalAlignment','bottom',...
       'Units','normalized',...
       'HorizontalAlignment','left',...
       'Position',[0 1 0]);
    htxsteps=xlabel('step number t=0');
   
   
   
   % == Comment Window frame ==============================================
%%   bottom=0.05;
%%   height=0.16;
%%   uicontrol( ...
%%        'Style','frame', ...
%%        'Units','normalized', ...
%%        'Position',[left bottom width height], ...
%%        'BackgroundColor',[0.5 0.5 0.5]);

   % Text label
%   uicontrol( ...
%        'Style','text', ...
%        'Units','normalized', ...
%        'Position',[left height-0.01 width 0.05], ...
%        'BackgroundColor',[0.5 0.5 0.5], ...
%        'ForegroundColor',[1 1 1], ...
%        'String','Comment Window');

   % Edit window
%%   border=0.01;
%%   hconsole=uicontrol( ...
%%        'Style','edit', ...
%%        'HorizontalAlignment','left', ...
%%        'Units','normalized', ...
%%        'Max',10, ...
%%        'BackgroundColor',[1 1 1], ...
%%        'Position',[left+border bottom width-2*border height-0.05], ...
%%      'Enable','inactive',...
%%        'String','');


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
   height=0.044;
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
      'Callback','demo_emgmm(''info'',gcf)',...
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
      'Callback','demo_emgmm(''step'',gcf)');

   % Stop button: stop process of adaptation
   bottom=bottom+height;
   hbtstop = uicontrol(...
    'Units','Normalized', ...
        'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Stop', ...
      'Callback','set(gcbo,''UserData'',1)',...
      'Enable','off');

   % Play button: begin adaptation
   bottom=bottom+height;
   hbtplay = uicontrol(...
    'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Play', ...
      'Callback','demo_emgmm(''play'',gcf)');

   % Reset button: set up t = 0
   bottom=bottom+height;
    hbtreset = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Reset', ...
      'Callback','demo_emgmm(''reset'',gcf)');

   % Create data
   bottom=bottom+1.5*height;
    hbtcreat = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Create data', ...
      'Callback','demo_emgmm(''creatdata'',gcf)');

   % Load data
   bottom=bottom+1*height;
    hbtload = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Load data', ...
      'Callback','demo_emgmm(''getfile'',gcf)');

   % Save model
   bottom=bottom+1.5*height;
    hbtSaveModel = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Save model', ...
      'Callback','demo_emgmm(''savemodel'',gcf)');

   % Load model
%   bottom=bottom+1*height;
%    hbtLoadModel = uicontrol(...
%      'Units','Normalized', ...
%      'ListboxTop',0, ...
%        'Position',[left bottom width height], ...
%      'String','Load model', ...
%      'Callback','demo_emgmm(''loadmodel'',gcf)');

   % == PopUp Menu =====================================================
   bottom=0.95-height;
   htxfeatures=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width height], ...
      'String','Covariance');
   % popup menu
   bottom=bottom-height;
   hpufeatures=uicontrol( ...
      'Style','popup', ...
      'Units','normalized', ...
      'Position',[left bottom width height], ...
      'String',['Diagonal '; 'Full     ']);

   % == Edit line ==========================================================
   % prior info about number of the classes
   bottom=bottom-1.3*height;
   htxclasses=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'String','Components');
   bottom=bottom-height;
   hedclasses = uicontrol(...
    'Units','normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'Style','edit',...
      'String','2');

   % Iterations
   bottom=bottom-1.3*height;
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

   % == CheckBox =========================================================

   % Should the first algorithm step be random or not ?
   bottom=bottom-height*1.3;
    hxbrandom = uicontrol(...
    'Style','checkbox', ...
       'Units','normalized', ...
    'ListboxTop',0, ...
       'Position',[left bottom width height], ...
    'String','Random init');

   %=====================================================================
   % Store handlers
   handlers=struct(...
      'ellipse',struct('handler',-1,'mi',[],'sigma',[],'t',0,'Pk',[],...
        'solution',0),...
      'center',[],...
      'graph1',struct('handler',-1,'loglik',[],'axist',0,'time',[]),...
      'title1',htitle1,...
      'title2',htitle2,...
      'btSaveModel',hbtSaveModel,...%      'btLoadModel',hbtLoadModel,...
      'btstep',hbtstep,...
      'btstop',hbtstop,... 
      'btclose',hbtclose,...
      'btplay',hbtplay,...
      'btreset',hbtreset,...
      'btinfo',hbtinfo,...
      'btload',hbtload,...
      'btcreat',hbtcreat,... %%%      'console',hconsole,...
      'txsteps',htxsteps,...
      'txclasses',htxclasses,...
      'txiter',htxiter,...
      'txfeatures',htxfeatures,...
      'pufeatures',hpufeatures,...
      'editer',hediter,...
      'xbrandom',hxbrandom,...
      'axes1',haxes1,...
      'axes2',haxes2,...
      'edclasses',hedclasses);
   set(hfigure,'UserData',handlers)

   % Reset
   demo_emgmm('reset',hfigure);

   % Put figure on desktop
   set(hfigure,'Visible','on');
   drawnow;

 case 'savemodel'
   % == Save model ============================================

   h=get(hfigure,'UserData');
   
   if h.ellipse.t == 0,
     errordlg('No model has found yet.','No model to save','modal');
     return;
   end
   
   [name,path]=uiputfile('*.mat','Save model');
   if name ~= 0,
         
     Mean = h.ellipse.mi;
     Cov = reshape(h.ellipse.sigma,2,2,size(Mean,2));
     Prior = h.ellipse.Pk;
     fun = 'pdfgmm';
     save(strcat(path,name),'Mean','Cov','Prior','fun');
     
   end   
   
 case 'loadmodel'
   % == Load model ============================================

   h=get(hfigure,'UserData');
   sets=get(h.axes1,'UserData');

   % Are data sets loaded ?
   if isempty(sets)==1,
      % no warning is needy because of the huge topic on the screen:))).
      return;
   end
   
   [name,path]=uigetfile('*.mat','Load model');
   if name ~= 0,

     fname=strcat(path,name);
%     if checkdat(fname,DATA_IDENT,2,0)==1,

       demo_emgmm('reset',hfigure);
       
       model=load(fname);
       
       if exist('model.Pk')==0,
         % suppose uniformly distributed Pk
         model.Pk=ones(1,sum(model.K))/sum(model.K);
       end
       
       h.ellipse.mi = model.MI;
       h.ellipse.sigma=model.SIGMA;
       h.ellipse.Pk = model.Pk;
       
       h.ellipse.solution=0;
       h.ellipse.t=0;

       set(h.edclasses,'String',num2str(sum(model.K)));
       h.ellipse.classes=sum(model.K);
       h.ellipse.features=get(h.pufeatures,'Value');
                     
       set(hfigure,'UserData',h);                
       
       demo_emgmm('step',hfigure);
 %    else
%      errordlg('This file does not contain required data.','Bad file','modal');
%     end
   end
 
   
 case 'play'
   % == Play ============================================
   h=get(hfigure,'UserData');

   % get data set
   sets=get(h.axes1,'UserData');

   % Are data sets loaded ?
   if isempty(sets)==1 | h.ellipse.solution==1,
%%      text=sprintf(...
%%        'No action performed. No data to work on. Load or create it!');
%%      set(h.console,'String',text);
      return;
   end

   % disable button
   set([h.editer,h.btstep,h.btclose,h.btplay,...
        h.btreset,h.btinfo,h.btload,h.btcreat,h.txiter],...
      'Enable','off');

   % enable stop button
   set(h.btstop,'Enable','on');

   % get # of iterations
   iter=str2num(get(h.editer,'String'));

   % # of classes
   if h.ellipse.t==0,
      h.ellipse.classes=str2num(get(h.edclasses,'String'));
      h.ellipse.features=get(h.pufeatures,'Value');
      set([h.xbrandom,h.edclasses,h.txclasses,h.pufeatures,h.txfeatures],...
        'Enable','off');
   end

   % Shall the init be random, yes or no ?
   randinit=get(h.xbrandom,'Value');

   % set stop button
   set(h.btstop,'UserData',0);


   % Play - adaptation process
   play=1; % flag 1 - not finished, 1 - finished;
   while play==1 & get(h.btstop,'UserData')==0,
      
     options.rand = randinit;
     switch 3-h.ellipse.features,
       case 1,
         options.cov_type = 'full';
       case 2,
         options.cov_type = 'diag';
     end
%     options.cov_type = 3-h.ellipse.features;
     options.ncomp=h.ellipse.classes;
     options.tmax = iter+h.ellipse.t;
     if h.ellipse.t >= 1,
      init_model.logL = h.ellipse.logL;
      init_model.Alpha = h.ellipse.Alpha;
      init_model.t = h.ellipse.t;
      init_model.Mean = h.ellipse.mi;
      init_model.Cov = reshape(h.ellipse.sigma,2,2,size(h.ellipse.mi,2));
      init_model.Prior = h.ellipse.Pk;
      model=emgmm(sets.X,options,init_model);
     else
       model=emgmm(sets.X,options);
     end
     
     h.ellipse.mi = model.Mean;
     h.ellipse.sigma = reshape( model.Cov,2,2*size(model.Mean,2));
     h.ellipse.Pk = model.Prior;
     [tmp,eI]=max(model.Alpha);
     h.ellipse.solution = model.exitflag;
     h.ellipse.t = model.t;
     h.ellipse.Alpha = model.Alpha;
     h.ellipse.logL = model.logL(end);
      
      % perform one learning step
%      if h.ellipse.features==2,   % correlated features
%        [h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk,eI,h.ellipse.solution,...
%         h.ellipse.t]=unsund(sets.X,h.ellipse.classes,iter,randinit,...
%         h.ellipse.t,h.ellipse.mi, h.ellipse.sigma,h.ellipse.Pk);
%      else % independent
%         [h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk,eI,h.ellipse.solution,...
%         h.ellipse.t]=unsuni(sets.X,h.ellipse.classes,iter,randinit,...
%         h.ellipse.t,h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk);
%      end

      text=sprintf('step number t=%d ',h.ellipse.t);
      if h.ellipse.solution==1,
%%         text=strvcat(text,'EM has converged.');
         text=[text ',EM has converged.'];
         play=0;
         set(h.txsteps,'String',text);
      else
         set(h.txsteps,'String',text);

         val= mln(sets.X,h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk);

         text=sprintf('Log-likelihood, logL(t) = %f',val);
         set(h.title2,'String',text);
      
         h.graph1.time=[h.graph1.time,h.ellipse.t];
         h.graph1.loglik=[h.graph1.loglik,val];

         ylimit=get(h.axes2,'YLim');
         if ylimit(2) < val,
            set(h.axes2,'YLim',[ylimit(1) val+AXISY_ADD]);
         end        
      
         % is axis to be changed ?
         if h.ellipse.t > h.graph1.axist,
            h.graph1.axist=h.ellipse.t+iter*AXIST_ADD;
            set(h.axes2,'XLim',[1 h.graph1.axist]);
         end
         set(h.graph1.handler,'XData',h.graph1.time,'YData',h.graph1.loglik,...
           'Visible','on');    
         
      
         if h.ellipse.handler==-1,
            axes(h.axes1);
         end
         [h.ellipse.handler,h.center]=...
          pnmix(sets.X,h.ellipse.mi,h.ellipse.sigma,eI,h.ellipse.handler,h.center);
      end

      % comment
%%      set(h.console,'String',text);

      % store data
      set(hfigure,'UserData',h);

      % flush it on desktop
      drawnow;
   end % of while

   % dissable button
   set([h.editer,h.btstep,h.btclose,h.btplay,...
        h.btreset,h.btinfo,h.btload,h.btcreat,h.txiter],...
      'Enable','on');
   % enable stop button
   set(h.btstop,'Enable','off');

   % copy the solution to global variables to be visible outside demo_emgmm
   UNSU_MI = h.ellipse.mi;
   UNSU_SIGMA = h.ellipse.sigma;
   UNSU_PK = h.ellipse.Pk;
   UNSU_solution = h.ellipse.solution;
   UNSU_t = h.ellipse.t;

case 'step'
   % == One step of learning ===========================================
   h=get(hfigure,'UserData');

   % get data set
   sets=get(h.axes1,'UserData');

   % are data sets loaded ?
   if isempty(sets)==1 | h.ellipse.solution==1,
%%text=sprintf('No action performed. No data to work on. Load or create it!');
%%      set(h.console,'String',text);
%%      set(h.txsteps,'String',text);
      return;
   end

   % get # of iter
   iter=str2num(get(h.editer,'String'));

   % # of classes
   if h.ellipse.t==0,
      h.ellipse.classes=str2num(get(h.edclasses,'String'));
      h.ellipse.features=get(h.pufeatures,'Value');
      set([h.xbrandom,h.edclasses,h.txclasses,h.pufeatures,h.txfeatures],'Enable','off');
   end

   % random init yes or no
   randinit=get(h.xbrandom,'Value');

   % perform one learning step
%   if h.ellipse.features==2,   % correlated
%       [h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk,eI,h.ellipse.solution,...
%        h.ellipse.t]=unsund(sets.X,h.ellipse.classes,iter,randinit,...
%        h.ellipse.t,h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk);
%  else % independent
     options.rand = randinit;
     switch 3-h.ellipse.features,
       case 1,
         options.cov_type = 'full';
       case 2,
         options.cov_type = 'diag';
     end
     
%     options.cov_type = 3-h.ellipse.features;
     options.ncomp=h.ellipse.classes;
     options.tmax = iter+h.ellipse.t;
     if h.ellipse.t >= 1,
      init_model.logL = h.ellipse.logL;
      init_model.Alpha = h.ellipse.Alpha;
      init_model.t = h.ellipse.t;
      init_model.Mean = h.ellipse.mi;
      init_model.Cov = reshape(h.ellipse.sigma,2,2,size(h.ellipse.mi,2));
      init_model.Prior = h.ellipse.Pk;
      model=emgmm(sets.X,options,init_model);
     else
       model=emgmm(sets.X,options);
     end
     
     h.ellipse.mi = model.Mean;
     h.ellipse.sigma = reshape( model.Cov,2,2*size(model.Mean,2));
     h.ellipse.Pk = model.Prior;
     [tmp,eI]=max(model.Alpha);
     h.ellipse.solution = model.exitflag;
     h.ellipse.t = model.t;
     h.ellipse.Alpha = model.Alpha;
     h.ellipse.logL = model.logL(end);
%      [h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk,eI,h.ellipse.solution,...
%         h.ellipse.t]=unsuni(sets.X,h.ellipse.classes,iter,randinit,...
%      h.ellipse.t,h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk);
%   h.ellipse
%   end

   text=sprintf('step number t=%d ',h.ellipse.t);
   if h.ellipse.solution==1,
%%      text=strvcat(text,'Solution is found.');
        text=[text ',EM has converged.'];
      set(h.txsteps,'String',text);
    else
      set(h.txsteps,'String',text);

      val= mln(sets.X,h.ellipse.mi,h.ellipse.sigma,h.ellipse.Pk);

      text=sprintf('Log-likelihood, logL(t) = %f',val);
      set(h.title2,'String',text);
      
      h.graph1.time=[h.graph1.time,h.ellipse.t];
      h.graph1.loglik=[h.graph1.loglik,val];

      ylimit=get(h.axes2,'YLim');
      if ylimit(2) < val,
         set(h.axes2,'YLim',[ylimit(1) val+AXISY_ADD]);
      end        
      
      % is axis to be changed ?
      if h.ellipse.t > h.graph1.axist,
         h.graph1.axist=h.ellipse.t+iter*AXIST_ADD;
         set(h.axes2,'XLim',[1 h.graph1.axist]);
      end
      set(h.graph1.handler,'XData',h.graph1.time,'YData',h.graph1.loglik,...
        'Visible','on');    
      
      if h.ellipse.handler==-1,
         axes(h.axes1);
      end            

      [h.ellipse.handler,h.center]=...
      pnmix(sets.X,h.ellipse.mi,h.ellipse.sigma,eI,h.ellipse.handler,h.center);
   end

   % comment
%%   set(h.console,'String',text);

   % flush it on desktop
   drawnow;

   set(hfigure,'UserData',h);

   % copy the solution to global variables to be visible outside demo_emgmm
   UNSU_MI = h.ellipse.mi;
   UNSU_SIGMA = h.ellipse.sigma;
   UNSU_PK = h.ellipse.Pk;
   UNSU_solution = h.ellipse.solution;
   UNSU_t = h.ellipse.t;


case 'getfile'
   % == Invoke standard open file dialog ====================================
   % Opens file and checks if contains apropriate data, if yes than loads data.

   h=get(hfigure,'UserData');

   % change path to directory
%%   wres=what('unsuper');
%%   cd(wres.path);

   [name,path]=uigetfile('*.mat','Open file');
   if name~=0,
      file.pathname=strcat(path,name);
      file.path=path;
      file.name=name;
%      if checkdat(file.pathname,DATA_IDENT,2,0)==1,
      if check2ddata(file.pathname)==1,
         set(h.btload,'UserData',file);
         demo_emgmm('loadsets',hfigure);
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
   sets.K = size(sets.X,2);
   sets.N = size(sets.X,1);

   % store loaded sets
   set(h.axes1,'UserData',sets);

   % call reset
   demo_emgmm('reset',hfigure);

   drawnow;


case 'reset'
   % == Reset adaptation process, set up t=0 ================

   h=get(hfigure,'UserData');                     % get handlers

   % get file
   file=get(h.btload,'UserData');

   % get data set
   sets=get(h.axes1,'UserData');

   % zeroes parameters of the separation line
   h.ellipse.mi=[];
   h.ellipse.sigma=[];
   h.ellipse.t=0;
   h.ellipse.Pk=[];
   h.center=-1;
   h.ellipse.handler=-1;
   h.ellipse.solution=0;

   h.graph1.time=[];
   h.graph1.axist=0;
   h.graph1.loglik=[];

   % clear axes prob.
   clrchild(h.axes2);
   axes(h.axes2);
%%   setaxis(h.axes2,[0 1 0 1]);
   axis auto;
   h.graph1.handler=plot([0],[0],'b','Parent',h.axes2,...
      'EraseMode','background','Visible','off');
         
   % store h
   set(hfigure,'UserData',h);

   % enable Edit Line 'classes'
   set([h.xbrandom,h.edclasses,h.txclasses,h.pufeatures,h.txfeatures],'Enable','on');

   % comment
   text=sprintf('Log-likelihood, logL(t)');
   set(h.title2,'String',text);

   text=sprintf('step number t=0. ');
   set(h.txsteps,'String',text);

   % clears axes
   set(get(h.axes1,'Children'),'EraseMode','normal');
   %%%   cla;
   clrchild(h.axes1);
   drawnow;

   % set axes and plot mixture
   axes(h.axes1);
   if isempty(sets)==0,
      win=cmpwin(min(sets.X'),max(sets.X'),BORDER,BORDER);
      %%%      axis(win);
      setaxis(h.axes1,win);
%%      axes(h.axes1);

%%      ppoints(sets.X,sets.I);
      ppatterns(sets.X);
   end
   %%%%%%%%%%%%%%%%%%%%%%%%%%%%


   % create comment
   if isempty(sets)==0,
      set(h.title1,'String',sprintf('File: %s, # of points K = %d',file.name,sum(sets.K)));
   else
      set(h.title1,'String','No data loaded');

      pos=get(h.axes1,'Position');
      fsize=min(pos(3),pos(4))/7;
      setaxis(h.axes1,[-1 1 -1 1]);
      builtin('text',0,0,'Press ''Load data'' button.',...
         'Parent',h.axes1,...
         'HorizontalAlignment','center',...
         'FontUnits','normalized',...
         'Clipping','on',...
         'FontSize',fsize);
      builtin('text',0,-fsize*2,...
         'Load sample data from ../toolboxroot/data/gmm\_samples/ ',...
         'Parent',h.axes1,...
         'HorizontalAlignment','center',...
         'FontUnits','normalized',...
         'Clipping','on',...
         'FontSize',fsize*0.65);
   end

   drawnow;

case 'creatdata'
   % == Invoke data set creator ============================================
   createdata('finite',10,'demo_emgmm','created',hfigure);

case 'created'
   % == Load new created data set ===========================================

   % get handler and make this figure active
   figure(hfigure);
   h=get(hfigure,'UserData');

   % get file name
   path=varargin{1};
   name=varargin{2};
   pathname=strcat(path,name);

%   if checkdat(pathname,DATA_IDENT,2,0)==1,
   if check2ddata(pathname)==1,
      file.pathname=pathname;
      file.path=path;
      file.name=name;
      set(h.btload,'UserData',file);
      demo_emgmm('loadsets',hfigure);
   else
      errordlg('This file does not contain required data.','Bad file','modal');
   end

case 'info'
   % == Call standard Matlab`s info box =========================================
   helpwin(mfilename);

end % of switch


function []=clrchild(handle)
delete(get(handle,'Children'));
return;

function []=setaxis(handle,rect)
set(handle,'XLim',rect(1:2));
set(handle,'YLim',rect(3:4));

if size(rect,2)>=6,
   set(handle,'ZLim',rect(5:6));
end

return;

function [win]=cmpwin(mins,maxs,xborder,yborder)
dx=max( (maxs(1)-mins(1)), 1 )*xborder;
dy=max( (maxs(2)-mins(2)), 1 )*yborder;

x1=(mins(1)-dx);
x2=(maxs(1)+dx);
y1=(mins(2)-dy);
y2=(maxs(2)+dx);

win=[x1 x2 y1 y2];

function [logL] = mln(X,MI,SIGMA,Pk)
 

D=size(MI,1);   % dimension
K=size(MI,2);   % % of classes
N=size(X,2);    % # of points

A=zeros(N,K);

for k=1:K,
    pxk=normald(X,MI(:,k),SIGMA(:,1+(k-1)*D:k*D));

    A(:,k)=pxk(:)*Pk(k);
end

logL=sum(log(sum(A,2)));


function [p]=normald(X,mi,sigma)

DIM=size(X,1);
p=exp(-1/2*mahalan(X,mi,sigma))/((2*pi)^(DIM/2) * sqrt(det(sigma)));

function [hellipse,hcenter]=pnmix(X,MI,SIGMA,I,hellipse,hcenter)
%  [hellipse,hcenter]=pnmix(X,MI,SIGMA,I,hellipse,hcenter)
%
% PNMIX vizualizes mixture of normal distributions in 2D space. 
%  Each normal distribution is determined by a pair of mean values and 
%  covariance matrix. The mean value is vizualized as a point and the 
%  covariance matrix as en ellipse. A size of the ellipse is determined 
%  in order to contain all the points belonging to given class which 
%  the ellipse describes.
%
%  Input and output arguments hellipse and hcenter contain
%  handles of graphics objects (ellipses and their centers) and
%  if they enter the function, the old graphics objects vanish and 
%  then new objects are plotted.
%
%  The function is useful for vizualization of results of the unsupervised
%  learning algorithms (see help of UNSUNI, UNSUND).
%
% Input:
%   X [NxK] contains K points which are N-dimensional, X=[x_1,x_2,...,x_K].
%   I [1xK] contains class labels for all the points.
%   MI [NxM] mean values for each class, MI=[mi_1,mi_2,...,mi_M]
%   SIGMA [Nx(MxN)] covariance matrices for each class,
%      SIGMA=[sigma_1,sigma_2,...,sigma_M].
%   hellipse [vector], hcenter [vector] handlers of graphics objects.
%
% Output:
%   hellipse [vector], hcenter [vector] handlers of graphics objects.
%
% See also UNSUNI, UNSUND, UNSUDEMO.
%


% Statistical Pattern Recognition Toolbox, Vojtech Franc, Vaclav Hlavac
% (c) Czech Technical University Prague, http://cmp.felk.cvut.cz
% Written Vojtech Franc (diploma thesis) 10.11.1999, 23.12.1999
% Modifications:
% 22. 6.00 V. Hlavac, comments polished.


if nargin < 5,
  hellipse=-1;
end

DIM=size(X,1);
N=size(X,2);
K=size(MI,2);
maxr=zeros(1,K);
for i=1:N,
   r=sqrt(mahalan(X(:,i),MI(:,I(i)),SIGMA(:,(I(i)-1)*DIM+1:DIM*I(i))));
   if maxr(I(i)) < r,
      maxr(I(i)) = r;
   end
end

if hellipse==-1,
   for i=1:K,
%%%      [x,y]=ellipse(inv(SIGMA(:,(i-1)*DIM+1:DIM*i)),30,maxr(i),MI(:,i));
      [x,y]=ellips(MI(:,i),inv(SIGMA(:,(i-1)*DIM+1:DIM*i)),maxr(i),30);
      hellipse(i)=plot(x,y,'k','EraseMode','xor');
      hcenter(i)=plot(MI(1,i),MI(2,i),'+k','EraseMode','xor');
      drawnow;
   end
else
   for i=1:K,
%%%      [x,y]=ellipse(inv(SIGMA(:,(i-1)*DIM+1:DIM*i)),30,maxr(i),MI(:,i));
      [x,y]=ellips(MI(:,i),inv(SIGMA(:,(i-1)*DIM+1:DIM*i)),maxr(i),30);
      set(hellipse(i),'XData',x,'YData',y,'Visible','on');
      set(hcenter(i),'XData',MI(1,i),'YData',MI(2,i),'Visible','on');
      drawnow;
   end
end
