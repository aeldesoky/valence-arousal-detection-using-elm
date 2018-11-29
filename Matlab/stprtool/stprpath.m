function stprpath(toolboxroot)
% STPRPATH sets path to Statistical Pattern Recognition Toolbox. 
%
% Synopsis:
%  stprpath
%  stprpath(toolboxroot)
%
% Description:
%  stprpath(toolboxroot) sets path to the Statistical Pattern 
%   Recognition Toolbox stored in given root directory toolboxroot.
%
%  stprpath uses toolboxroot = pwd .
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2005, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 23-nov-2005, VF, MAC computers problem removed - updated version received 
%                  from Vivek Varshney
% 24-jan-2005, VF, added a new directory 'optimization'
% 28-apr-2004, VF, renamed to stprpath
% 22-oct-2003, FV, command addpath used.
% 18-July-2003, VF
% 9-Feb-2003, VF
% 23-Jan-2003, VF
% 7-jan-2003, VF, A new coat.
% 23-mar-2001, V.Franc, added new directories

if nargin < 1
   toolboxroot=pwd;              % get current directory
end

disp('Adding path for the Statistical Pattern Recognition Toolbox...');

% path for UNIX
p = ['$:',...
     '$bayes:',...
     '$data:',...
     '$data/anderson_task:',...
     '$data/binary_separable:',...
     '$data/gmm_samples:',...
     '$data/mm_samples:',...
     '$data/multi_separable:',...
     '$data/svm_samples:',...
     '$data/iris_data:',...
     '$data/riply_data:',...
     '$demos:',...
     '$demos/ocr:',...
     '$demos/image_denoising:',...
     '$kernels:',...
     '$kernels/extraction:',...
     '$kernels/preimage:',...
     '$linear:',...
     '$linear/anderson:',...
     '$linear/finite:',...
     '$linear/fisher:',...
     '$linear/extraction:',...
     '$misc:',...
     '$optimization:',...
     '$probab:',...
     '$probab/estimation:',...
     '$quadrat:',...
     '$svm:',...
     '$visual:',...
     '$xtal_regression:',...
    ];

p=translate(p,toolboxroot);

% adds path at the start
addpath(p);



%--translate ---------------------------------------------------------
function p = translate(p,toolboxroot);
%TRANSLATE Translate unix path to platform specific path
%   TRANSLATE fixes up the path so that it's valid on non-UNIX platforms
%
% This function was derived from MathWork M-file "pathdef.m"

cname = computer;
% Look for VMS, this covers VAX_VMSxx as well as AXP_VMSxx.
%if (length (cname) >= 7) & strcmp(cname(4:7),'_VMS')
%  p = strrep(p,'/','.');
%  p = strrep(p,':','],');
%  p = strrep(p,'$toolbox.','toolbox:[');
%  p = strrep(p,'$','matlab:[');
%  p = [p ']']; % Append a final ']'

% Look for PC
if strncmp(cname,'PC',2)
  p = strrep(p,'/','\');
  p = strrep(p,':',';');
  p = strrep(p,'$',[toolboxroot '\']);

% Look for MAC
%elseif strncmp(cname,'MAC',3)
%  p = strrep(p,':',':;');
%  p = strrep(p,'/',':');
%  m = toolboxroot;
%  if m(end) ~= ':'
%    p = strrep(p,'$',[toolboxroot ':']);
%  else
%    p = strrep(p,'$',toolboxroot);
%  end
else
  p = strrep(p,'$',[toolboxroot '/']);
end
