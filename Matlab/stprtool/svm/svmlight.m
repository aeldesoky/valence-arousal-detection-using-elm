function model=svmlight(data,options)
% SVMLIGHT Interface to SVM^{light} software.
%
% Synopsis:
%  model = svmlight(data)
%  model = svmlight(data,options)
%
% Description:
%  This function serves as an interface between Matlab 
%  and SVM^{light} (Version: 5.00) optimizer which trains 
%  the Support Vector Machines classifier.
%
%  The executable file 'svm_learn' must be in the path. 
%  The SVM^{light} software can be downloaded from:
%      http://svmlight.joachims.org/
%
%  This function creates temporary files 'tmp_alphaXX.txt', 
%  'tmp_examplesXX.txt', 'tmp_modelXX.txt' and 'tmp_verbXX.txt' for 
%  comunication with the SVM^{light} software. The XX=datestr(now)
%  is string consisting of current date and time.
%           
% Input:
%  data [struct] Labeled binary data:
%   .X [dim x num_data] Training vectors.
%   .y [1 x num_data] Labels of training data (1 or 2).
%  
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier: 
%      'linear' (default),'rbf' and 'poly'. 
%   .arg [1x1] Kernel argument (default []).
%   .C [1x1] SVM regularization constant (default C=inf).
%   .mC [1x1] if mC is given then C is set to mC/length(data.y). 
%   .j [1x1] Cost-factor, by which training errors on 
%     positive examples outweight errors on negative examples (default 1).
%   .eps [1x1] Tolerance of KKT-conditions (default eps=0.001).
%   .b [1x1] if 1 (default) then finds w'*x +b else b = 0;
%   .keep_files [1x1] If ==1 then keeps temporary files otherwise
%     erase them.
%   .svm_command [string] Path to SVM^{light} solver (default "svm_learn")
%
% Output:
%  model [struct] Binary SVM classifier:
%   .Alpha [nsv x 1] Weights of support vectors.
%   .b [1x1] Bias of decision function.
%   .sv.X [dim x nsv] Support vectors.
%   .sv.inx [1 x nsv] Indices of SVs (model.sv.X = data.X(:,inx)).
%   .nsv [int] Number of Support Vectors.
%   .kercnt [int] Number of kernel evaluations used by the SVM^{light}.
%   .trnerr [real] Classification error on training data.
%   .margin [real] Margin of found classifier.
%   .options [struct] Copy of used options.
%   .cputime [real] Used CPU time in seconds.
%
% Example:
%  data=load('riply_trn');  
%  model=svmlight(data,struct('ker','rbf','C',10,'arg',1))
%  figure; ppatterns(data); psvm(model);
%
% See also 
%  SVMCLASS, XY2SVMLIGHT.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>


% Modifications:
% 09-sep-2007, VF, -b option added
% 21-may-2007, VF, -q 42 (size of QP subproblem) added based on Soeren's suggestion
% 20-nov-2006, VF, added optional parameter mC
% 10-oct-2006, VF, "svm_command" option added
% 09-feb-2006, VF, added date_str(findstr(date_str,':')) = '.'; based on
%   M.Urban comment.
% 16-may-2004, VF
% 15-jan-2004, VF, handling argument of poly kernel repared
% 10-oct-2003, VF, computation of lin model added
% 29-aug-2003, VF, seconds are added to the name of temporary files
% 12-may-2003, VF, 1st 3 lines of verb_file are skiped
% 31-jan-2003, VF, added option 'j' 
% 28-Jan-2003, VF
% 20-jan-2003, VF, temporary files are unique and are deleted at the end
% 14-Jan-2003, VF
% 26-sep-2002, VF
% 3-Jun-2002, V.Franc

tic;

data=c2s(data);

% gets current date and time
date_str=datestr(now);
date_str(findstr(date_str,' ')) = '-';
date_str(findstr(date_str,':')) = '.';
sec=clock;
date_str = [date_str '-' num2str(sec(end))];


% names of temporary files 
examples_file = ['tmp_examples' date_str '.txt'];
model_file = ['tmp_model' date_str '.txt'];
verb_file = ['tmp_verb' date_str '.txt'];
alpha_file = ['tmp_alpha' date_str '.txt'];

% make model
model.name = 'SVM classifier';

% -- Process input arguments --------------------------
if nargin < 2, options = []; else options=c2s(options); end

if isfield(options,'mC'), options.C = options.mC/length(data.y); end
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = '[]'; end
if ~isfield(options,'C'), options.C = inf; end
if ~isfield(options,'eps'), options.eps = 0.001; end
if ~isfield(options,'keep_files'), options.keep_files = 0; end
if ~isfield(options,'j'), options.j = 1; end
if ~isfield(options,'b'), options.b = 1; end
if ~isfield(options,'svm_command'), options.svm_command = 'svm_learn'; end

% gets data dimensions
[dim,num_data ] = size(data);

%--------------------------------
switch options.ker
  case 'linear'
    ker='-t 0';
  case 'rbf'
    ker=['-t 2 -g ' num2str(1/(2*options.arg^2))]; 
  case 'poly' 
    if length(options.arg) == 1,
      ker=['-t 1 -r 1 -s 1 -d ' num2str(options.arg)];  
    else
      ker=['-t 1 -s 1 -r ' num2str(options.arg(2)) ' -d ' ...
           num2str(options.arg(1))];  
    end
end
command=[options.svm_command ' ' ...
         '-c ' num2str(options.C) ' '...
         ker ' ' ...
         '-v 1 ' ...
         '-m 40 ' ...
         '-q 42 ' ...
         '-j ' num2str(options.j) ' '...
         '-e ' num2str(options.eps) ' '...
         '-b ' num2str(options.b) ' ' ...
         '-a ' alpha_file ' ' examples_file ' ' model_file ' > ' verb_file];
   
% converts data to SVM^light format
xy2svmlight(data,examples_file);
    
% call svm_learn
[a,b]=system(command);   

% parses model file
checkfile(model_file); [lines]=textread(model_file,'%s');
for i=1:size(lines,1)
   if strcmpi( lines(i), 'threshold' )==1,
     model.b=-str2num( lines{i-2});
     break;
   end
end
    
checkfile(alpha_file); model.Alpha=textread(alpha_file,'%f');
model.Alpha=model.Alpha(:);
model.Alpha(find(data.y==2)) = -model.Alpha(find(data.y==2));

checkfile(verb_file);
[lines]=textread(verb_file,'%s',-1,'bufsize',5000000,'headerlines',3);
for i=1:size(lines,1)
  if strcmpi( lines{i}, 'misclassified,' ),
    model.trnerr=str2num( lines{i-1}(2:end));
    model.trnerr=model.trnerr/length(model.Alpha);
  end
  if strcmpi( lines(i), 'vector:' ) & strcmpi( lines(i-1), 'weight' )==1,
    tmp=str2num( lines{i+1}(5:end));
    if tmp~=0, model.margin=1/tmp; else model.margin=[]; end
  end
  if strcmpi( lines(i), 'SV:' )==1,
    model.nsv=str2num( lines{i+1});
  end
  if strcmpi( lines(i), 'evaluations:' )==1,
    model.kercnt=str2num( lines{i+1});
  end
end

model.nsv = length(find(model.Alpha~=0));

inx=find(model.Alpha);
model.sv.X = data.X(:,inx);
model.sv.y = data.y(inx);
model.sv.inx = inx;
model.Alpha = model.Alpha(inx);
model.Alpha(find(model.sv.y==2)) = -model.Alpha(find(model.sv.y==2));

if strcmp( options.ker, 'linear'),
  model.W = model.sv.X * model.Alpha;
end

model.options = options;
model.fun = 'svmclass';

% deletes temporary files
if options.keep_files == 0,
  delete(examples_file);
  delete(model_file);
  delete(verb_file);
  delete(alpha_file);
end

model.cputime=toc;

return;


function checkfile(fname)
% Check if file of given name exists. If not then prints
% error.
% 

attempts=0;
found = exist(fname);
while attempts < 5 && ~found
  found = exist(fname);
  attempts = attempts+1;  
end

if found == 0,
  error('File %s not found.\n', fname);
end

return;

%EOF



