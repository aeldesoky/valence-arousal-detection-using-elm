function result = demo_svm(action,hfigure,varargin)
% DEMO_SVM Demo on Support Vector Machines.
%
% Synopsis:
%  demo_svm
%
% Description:
%  DEMO_SVM demonstrates algorithms training the binary 
%   SVM classifier L1-soft and L2-soft margin [Vapnik95]
%   [Cris00]. The input training vectors must be 2-dimensional 
%   and can be interactively created by the user.
%
%  Following algorithms can be tested:
%
%  - Sequential Minimal Optimizer (SMO) for L1-norm soft margin.
%  - QP solver (quadprog) used to train SVM with L2-norm soft margin.
%  - Kernel Perceptron for separable hyperplane.
%
% Control:
%  Algorithm       - algorithm for testing.
%  Kernel          - non-linear kernel.
%  Kernel argument - argument of the non-linear kernel.
%  C-constant      - trade-off (regularization) constant.
%  parameters      - parameters of the selected algorithm.
%  background      - if selected then the background color
%     denotes the sign and the intenzity denotes the value 
%     of the found decision function.
%
%  FIG2EPS     - exports screen to the PostScript file.
%  Load data   - loads input training sets from file.
%  Create data - calls program for creating point sets.
%  Reset       - clears the screen.
%  Train SVM   - trains and displays the SVM classifer.
%  Info        - calls the info box.
%  Close       - close the program.
%
% See also 
%  SMO, SVMQUADPROG, KPERCEPTR.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 2-june-2004, VF
% 18-July-2003, VF
% 6-march-2002, V.Franc
% 23-oct-2001, V.Franc

BORDER=0.2;          % minimal space between axis and points
DATA_IDENT='Finite sets, Enumeration';  % file identifier
ALGOS=['SMO (L1)          ';...
       'QUADPROG (L2)     ';...
       'Kernel-Perceptron '];
KERNELS=['Linear    ';...
         'Polynomial';...
         'RBF       '];

SMO_PARAM = 'epsilon,tolerance';
DEF_SMO_PARAM = '1e-3,1e-3';
KERNELSK_PARAM = 'epsilon,iter_limit';
DEF_KERNELSK_PARAM = '1e-3,inf';
KPERCEPTR_PARAM = 'tmax';
DEF_KPERCEPTR_PARAM = 'inf';

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
   hfigure=figure('Name','Support Vector Machines', ...
      'Visible','off',...
      'NumberTitle','off', ...
      'Units','normalized', ...
      'Position',[left bottom width height],...
      'Units','normalized', ...
      'tag','demo_svm');

   % == Axes =========================================
   left=0.1;
   width=0.65;
   bottom=0.35;
   height=0.60;
   haxes=axes(...
      'Units','normalized', ...
      'UserData',[],...
      'Position',[left bottom width height]);
   xlabel('feature x');
   ylabel('feature y');

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
   
   % Export to EPS 
   width=0.1;
   left=0.75-width;
   bottom=0.95;
   height=0.04;
   hbt_close = uicontrol(...
    'Units','Normalized', ...
      'Callback','fig2eps(gcf)',...
        'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','FIG2EPS');

   % Close button
   left=0.8;
   bottom=0.05;
   height=0.045;
   width=0.15;
   hbt_close = uicontrol(...
    'Units','Normalized', ...
      'Callback','close(gcf)',...
        'ListboxTop',0, ...
        'Position',[left bottom width height], ...
        'String','Close');

   % Info button: call stanard info box
   bottom=bottom+1.5*height;
   hbt_info = uicontrol(...
    'Units','Normalized', ...
      'Callback','demo_svm(''info'',gcf)',...
        'ListboxTop',0, ...
        'Position',[left bottom width height], ...
        'String','Info');

   % Train SVM button
   bottom=bottom+1.5*height;
   hbt_train = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Train SVM', ...
      'Callback','demo_svm(''train'',gcf)');

   % Reset button
   bottom=bottom+height;
    hbt_reset = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Reset', ...
      'Callback','demo_svm(''reset'',gcf)');

   % Creat data
   bottom=bottom+1.5*height;
    hbt_creat = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Create data', ...
      'Callback','demo_svm(''creatdata'',gcf)');

   % Load data
   bottom=bottom+1*height;
    hbt_load = uicontrol(...
      'Units','Normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'String','Load data', ...
      'Callback','demo_svm(''getfile'',gcf)');


   % == Popup menus ======================================
   bottom=0.95-height;
   htx_algo=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width height], ...
      'String','Algorithm');
   % popup menu
   bottom=bottom-height*0.8;
   hpu_algo=uicontrol( ...
      'Style','popup', ...
      'Units','normalized', ...
      'CallBack','demo_svm(''algo_handler'',gcf)',...
      'Position',[left bottom width height], ...
      'String',ALGOS);

   % pop menu - kernel
   bottom=bottom-height*1.2;
   htx_kernel=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width height], ...
      'String','Kernel');
   % popup menu
   bottom=bottom-height*.8;
   hpu_kernel=uicontrol( ...
      'Style','popup', ...
      'Units','normalized', ...
      'CallBack','demo_svm(''kernel_handler'',gcf)',...
      'Position',[left bottom width height], ...
      'String',KERNELS);


   % == Edit line ========================================
   % kernel argument
   bottom=bottom-1.2*height;
   htx_arg=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'Enable','off',...
      'String','Kernel argument');
   bottom=bottom-height*.8;
   hed_arg = uicontrol(...
    'Units','normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'Style','edit',...
      'Enable','off',...
      'CallBack','demo_svm(''arg_handler'',gcf)',...
      'String','1');
   
   % C const
   bottom=bottom-1.2*height;
   htx_cconst=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'Enable','on',...
      'String','C-constant');
   bottom=bottom-height*.8;
   hed_cconst = uicontrol(...
    'Units','normalized', ...
      'ListboxTop',0, ...
        'Position',[left bottom width height], ...
      'Style','edit',...
      'Enable','on',...
      'CallBack','demo_svm(''cconst_handler'',gcf)',...
      'String','100');

   
   % parameters of the algortihm
   bottom=bottom-1.2*height;
   htx_param=uicontrol( ...
      'Style','text', ...
      'Units','normalized', ...
      'Position',[left bottom width 0.9*height], ...
      'Enable','on',...
      'String',SMO_PARAM);
   bottom=bottom-height*.8;
   hed_param = uicontrol(...
      'Units','normalized', ...
      'ListboxTop',0, ...
      'Position',[left bottom width height], ...
      'Style','edit',...
      'Enable','on',...
      'CallBack','demo_svm(''param_handler'',gcf)',...
      'String',DEF_SMO_PARAM);
   

   % == Check boxes ==============================================
   bottom=bottom-height*1.2;
    hxb_background = uicontrol(...
    'Style','checkbox', ...
    'Units','normalized', ...
    'ListboxTop',0, ...
    'Position',[left bottom width height], ...
    'String','Background');

   % ============================================================
   data=struct(...
      'bt_close',hbt_close,...
      'bt_train',hbt_train,...
      'bt_reset',hbt_reset,...
      'bt_info',hbt_info,...
      'bt_load',hbt_load,...
      'bt_creat',hbt_creat,...
      'pu_algo',hpu_algo,...
      'pu_kernel', hpu_kernel,...
      'ed_arg', hed_arg,...
      'tx_arg', htx_arg,...
      'tx_cconst', htx_cconst,...
      'ed_cconst', hed_cconst,...
      'tx_param', htx_param,...
      'ed_param', hed_param,...
      'console',hconsole,...
      'axes',haxes,...
      'xb_background',hxb_background);
   set(hfigure,'UserData',data );

   % Reset
   demo_svm('reset',hfigure);

   % Put figure on desktop
   set(hfigure,'Visible','on');
   drawnow;
   
%== Trains SVM and displays result ================================
case 'train'
   data = get( hfigure, 'UserData');
   
   trn = get( data.axes, 'UserData' );
   
   if isempty( trn ),
      return;
   end
   
   C = str2num(get( data.ed_cconst, 'String' ));
   ker_inx = get( data.pu_kernel, 'Value' );
   if ker_inx == 1,
      ker = 'linear';
   elseif ker_inx == 2;
      ker = 'poly';
   else 
      ker = 'rbf';
   end
   arg = str2num(get( data.ed_arg, 'String' ));
   
   [Alpha,bias] = svm_train( data, trn, ker, arg, C );
   
   % focus on the axes
   axes( data.axes );
   
   % Clear axes
   clrchild( data.axes);
   
   % get options
   options.background = get( data.xb_background, 'Value'); 
   
   % plot decision function
   model.sv.X = trn.X;
   model.sv.y = trn.I;
   ppatterns(model.sv); 
   inx=find( Alpha ~= 0 );
   model.sv.X = trn.X(:,inx);
   model.sv.y = trn.I(:,inx);
   model.Alpha = Alpha(:,inx)';
   model.b = bias;
   model.options.ker = ker;
   model.options.arg = arg;
   psvm( model, options );
   
%== Handler for Algorithm popup menu ==========================
case 'algo_handler'
   data=get(hfigure,'UserData');
   
   % which algorithm ?
   switch get(data.pu_algo, 'Value' )
     
    case 1   % SMO
        set(data.tx_cconst,'Enable','on');
        set(data.ed_cconst,'Enable','on');
        set(data.tx_param,'Enable','on');
        set(data.tx_param,'String',SMO_PARAM);
        set(data.ed_param,'Enable','on');
        set(data.ed_param,'String',DEF_SMO_PARAM);
        
%    case 2    % Matlab toolbox
%        set(data.tx_cconst,'Enable','on');
%        set(data.ed_cconst,'Enable','on');
%        set(data.tx_param,'Enable','off');
%        set(data.ed_param,'Enable','off');
     
%    case 3   % Kernel-SK
    case 2   % Kernel-SK
        set(data.tx_cconst,'Enable','on');
        set(data.ed_cconst,'Enable','on');
        set(data.tx_param,'Enable','on');
        set(data.ed_param,'Enable','on');
        set(data.tx_param,'String',KERNELSK_PARAM);
        set(data.ed_param,'String',DEF_KERNELSK_PARAM);
%    case 4    % Matlab toolbox
%        set(data.tx_cconst,'Enable','on');
%        set(data.ed_cconst,'Enable','on');
%        set(data.tx_param,'Enable','off');
%        set(data.ed_param,'Enable','off');
%    case 5    % Kernel Perceptron
    case 3    % Kernel Perceptron
        set(data.tx_cconst,'Enable','off');
        set(data.ed_cconst,'Enable','off');
        set(data.tx_param,'Enable','on');
        set(data.ed_param,'Enable','on');
        set(data.tx_param,'String',KPERCEPTR_PARAM);
        set(data.ed_param,'String',DEF_KPERCEPTR_PARAM);
   end
   
   
%== Handler for C-const edit line ===========================
case 'cconst_handler'
   data=get(hfigure,'UserData');

   C = str2num(get(data.ed_cconst,'String'));
   
   if C <= 0,
      C = 100;
   end
   set( data.ed_cconst,'String',num2str(C));
   
%== Handle for arg edit line ===========================
case 'arg_handler'
   data=get(hfigure,'UserData');

   arg = str2num(get(data.ed_arg,'String'));
   
   if arg < 0,
      arg = 1;
   end
   
   set( data.ed_arg, 'String',num2str( arg));

   
%== Handle for kernel pop up menu ===========================
case 'kernel_handler'
   data=get(hfigure,'UserData');

   ker_inx = get( data.pu_kernel,'Value');
   if ker_inx >= 2,
      set( data.ed_arg,'Enable','on');
      set( data.tx_arg,'Enable','on');
   else
      set( data.ed_arg,'Enable','off');
      set( data.tx_arg,'Enable','off');
   end

%== Calls data creator ==========================================
case 'creatdata'
   createdata('finite',2,'demo_svm','created',hfigure);
   
% == Loads recently created data ================================
case 'created'

   % get handler and make this figure active
   figure(hfigure);
   data = get(hfigure,'UserData');

   % get file name
   path=varargin{1};
   name=varargin{2};
   pathname=strcat(path,name);

   if check2ddata(pathname),                                                    
      file.pathname=pathname;                                                   
      file.path=path;                                                           
      file.name=name;                                                           
      set(data.bt_load,'UserData',file);
      demo_svm('loadsets',hfigure); 
      demo_svm('reset',hfigure); 
   else 
      errordlg('This file does not contain required data.','Bad file','modal'); 
   end 
    

% == Calls standard open file dialog ==========================
case 'getfile'
   
   data=get(hfigure,'UserData');

   [name,path]=uigetfile('*.mat','Open file');
   if name~=0,
      file.pathname=strcat(path,name);
      file.path=path;
      file.name=name;
%      if checkdat(file.pathname,DATA_IDENT,2,[0 0])==1,
      if check2ddata( file.pathname ),
         set(data.bt_load,'UserData',file);
         demo_svm('loadsets',hfigure);
         demo_svm('reset',hfigure);
      else
         errordlg('This file does not contain required data.','Bad file','modal');
      end
   end

case 'loadsets'
   % == Load data from file ==========================================

   data = get( hfigure,'UserData' );

   % Clear axes
   clrchild( data.axes);
   
   % set x and y axes labels
   xlabel('feature x');
   ylabel('feature y');

   % Get file name with sets
   file=get( data.bt_load,'UserData');

   % Load sets
   trn = load(file.pathname );
   trn.I=trn.y;
   trn.N= 2;
   trn.K = [length( find(trn.y==1)),length(find(trn.y==2))];

   % store loaded sets
   set( data.axes,'UserData', trn);

   % focus on axes
   axes( data.axes );

   % plots points
   ppatterns( trn );

   drawnow;
  

% == Reset ==========================================================
case 'reset'

   data = get(hfigure,'UserData');                 
   
   % Clear axes
   clrchild( data.axes);

   % get data set
   trn = get( data.axes, 'UserData');

   % get file
   file = get( data.bt_load,'UserData');
   
   % create comment
   if isempty( trn ) == 0,
     consoletext=sprintf('Data loaded.\nSelect algorithm and press Train SVM button.\n');
     titletext=sprintf('File: %s, # of points K = %d', file.name , size(trn.X,2));
     
     set( data.axes, 'XLimMode','auto', 'YLimMode','auto');
     ppatterns( trn);
     
   else
     consoletext=sprintf(['No data loaded.\n' ...
           'Press Create data button to create your own data.\n'...
           'Press Load data button to load data.\n' ...
           'Load sample data from ../data/binary/']);
     titletext='';

     pos=get( data.axes,'Position');
     fsize=min(pos(3),pos(4))/10;
     setaxis( data.axes,[-1 1 -1 1]);
%     axis([-1 1 -1 1]);
      
     builtin('text',0,0,'Press ''Load data'' button.',...
        'HorizontalAlignment','center',...
        'FontUnits','normalized',...
        'Clipping','on',...
        'FontSize',fsize);
   end

   % show comment
   set( data.console,'String',consoletext );

   % print title
   pos=get( data.axes,'Position');
   fsize=(1-pos(2)-pos(4))*1;
   title(titletext,...
      'Parent', data.axes,...
      'VerticalAlignment','bottom',...
      'HorizontalAlignment','left',...
      'FontUnits','normalized',...
      'Units','normalized',...
      'Position',[0 1 0],...
      'FontSize',fsize);


% == Calls Matlab`s info box ==================================
case 'info'
   helpwin(mfilename);
end

return;

%===============================================================
function [Alpha,bias] = svm_train( data, trn, ker, arg, C )

if strcmpi( ker, 'linear'),
  strarg = '-';
else
  strarg = num2str( arg);
end

param = str2num( get( data.ed_param,'String'));

switch get( data.pu_algo, 'Value' ),
  case 1   % Sequential Minimal Optimizer

%   if length(param) ~= length( str2num(DEF_SMO_PARAM)),
%     param = str2num(DEF_SMO_PARAM);
%   end
   
%   [Alpha,bias,nsv,kercnt,trn_err,margin]=...
%         smo(trn.X,trn.I,ker,arg,C, param(1),param(2) );
     options.ker = ker;
     options.arg = arg;
     options.C = C;
     options.eps = param(1);
     options.tol = param(2);
     
     model = smo( trn, options);
     
     Alpha = zeros(1,size(trn.X,2));
     Alpha(model.sv.inx) = model.Alpha(:)';
     bias = model.b;
     nsv = model.nsv;
     kercnt=model.kercnt;
     trn_err = model.trnerr;
     margin = model.margin;
          
     text = sprintf(...
        ['SVM (L1) by Sequential Minimal Optimizer\n',...
        'Kernel: %s (%s), C: %.4f\n',...
        'Kernel evaluations: %d\n',...
        'Number of Support Vectors: %d\n',...
        'Margin: %.4f\n',...
        'Training error: %.2f%%'],...
        ker, strarg, C, kercnt, nsv, margin, 100*trn_err );

%  case 2   % Matlab Optimization toolbox (L1)
%     model=svmmot(trn, {'ker',ker,'arg',arg,'C',C,'norm',1});
%
%     Alpha = zeros(1,size(trn.X,2));
%     Alpha(model.sv.inx) = model.Alpha;
%     bias = model.b;
%     nsv = model.nsv;
%     kercnt=model.kercnt;
%     trn_err = model.trnerr;
%     margin = model.margin;
%     exitflag=model.exitflag;
%     
%     text = sprintf(...
%        ['SVM (L1) using QUADPROG of the Optimization Toolbox\n',...
%        'Kernel: %s (%s), C: %.4f\n',...
%        'Exitflag: %d\n',...
%        'Kernel evaluations: %d\n',...
%        'Number of Support Vectors: %d\n',...
%        'Margin: %.4f\n',...
%        'Training error: %.2f%%'],...
%        ker, strarg, C, exitflag, kercnt, nsv, margin,100*trn_err);
%
% case 3 % Kernel Schlesinger-Kozinec's algorithm (L2)
  case 2   % Kernel Schlesinger-Kozinec's algorithm (L2)

      model=svmquadprog(trn, {'ker',ker,'arg',arg,'C',C,...
          'tolabs',param(1),'tolrel',0,'tmax',param(2) });

     Alpha = zeros(1,size(trn.X,2));
     Alpha(model.sv.inx) = model.Alpha(:)';
     bias = model.b;
     nsv = model.nsv;
     kercnt=model.kercnt;
     trn_err = model.trnerr;
     margin = model.margin;
     exitflag=model.exitflag;
    
     if( exitflag == 1), exitflag = 'found'; else exitflag = 'not found'; end
     
     text = sprintf(...
        ['SVM (L2) by Kernel Schlesinger-Kozinec`s algorithm\n',...
        'Kernel: %s (%s), C: %.4f\n',...
        'Solution: %s\n', ...
        'Number of Support Vectors: %d\n',...
        'Kernel evaluations: %d\n',...
        'Margin: %.4f\n',...
        'Training error: %.2f%%'],...
        ker, strarg, C, exitflag, nsv, kercnt, margin, trn_err*100 );

%  case 4   % Matlab Optimization toolbox (L2)
%     [Alpha,bias,nsv,exitflag,flps,margin,trn_err]=...
%         svm2mot(trn.X,trn.I,ker,arg,C);
%
%     model=svmmot(trn, {'ker',ker,'arg',arg,'C',C,'norm',2});
%
%     Alpha = zeros(1,size(trn.X,2));
%     Alpha(model.sv.inx) = model.Alpha;
%     bias = model.b;
%     nsv = model.nsv;
%     kercnt=model.kercnt;
%     trn_err = model.trnerr;
%     margin = model.margin;
%     exitflag=model.exitflag;
%
%     text = sprintf(...
%        ['SVM (L2) using QUADPROG of the Optimization Toolbox\n',...
%        'Kernel: %s (%s), C: %.4f\n',...
%        'Exitflag: %d\n',...
%        'Kernel evaluations: %d\n',...
%        'Number of Support Vectors: %d\n',...
%        'Margin: %.4f\n',...
%        'Training error: %.2f%%'],...
%        ker, strarg, C, exitflag, kercnt, nsv, margin,100*trn_err);
  
%  case 5   % Kernel Perceptron
  case 3   % Kernel Perceptron
     model=kperceptr(trn, {'ker',ker,'arg',arg,'tmax',param(1)});

     Alpha = zeros(1,size(trn.X,2));
     Alpha(model.sv.inx) = model.Alpha;
     bias = model.b;
     nsv = model.nsv;
     kercnt=model.kercnt;
     trn_err = model.trnerr;
     exitflag=model.exitflag;

     text = sprintf(...
        ['Kernel Perceptron\n',...
        'Kernel: %s (%s)\n',...
        'Exitflag: %d\n',...
        'Kernel evaluations: %d\n',...
        'Number of Support Vectors: %d\n',...
        'Training error: %.2f%%'],...
        ker, strarg, exitflag, kercnt, nsv, 100*trn_err);
end

% show comment
set( data.console,'String', text );

return;
     
function []=clrchild(handle)
% function []=clraxis(handle)
%
% CLRCHILD clears children of an object with the given handle.

delete(get(handle,'Children'));

return;

function []=setaxis(handle,rect)
% function []=setaxis(handle,rect)
%

set(handle,'XLim',rect(1:2));
set(handle,'YLim',rect(3:4));

if size(rect,2)>=6,
   set(handle,'ZLim',rect(5:6));
end

return;
