function compilemex(root)
% COMPILEMEX Compiles all MEX files of the STPRtoolbox.
%
% Synopsis:
%  compilemex
%  compilemex( toolboxroot )
%
% Description:
%  It calls MEX complier on all C-codes of the STPRtoolbox.
%  Run this function from the STPRtoolbox root directory or
%  specify the root directory as an input argument.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 23-sep-2008, VF: fix of path traslation of MacOS; suggested by dkim@mrn.org
% 30-mar-2007, VF
% 26-mar-2007, VF
% 20-nov-2006, VF, added compilation of "qpbsvm_mex.c" and "qpssvm_mex.c"
% 09-sep-2005, VF, added compilation of "gnnls_mex.c" and "gnpp_mex.c"
% 25-aug-2005, VF
% 24-jan-2005, VF
% 29-dec-2004, VF, inconsistent variables ("root" and "Root") unified
% 29-nov-2004, VF
% 19-sep-2004, VF
% 16-may-2004, VF
% 5-July-2003, VF
% 20-June-2003, VF
% 23-Jan-2003, VF

fprintf('Compiling MEX files of STPRtool...\n');

if nargin < 1
   root=pwd;              % get current directory
end

% -- List of functions to be complied ---------------------------
fun(1).source={'$kernels/kernel.c','$kernels/kernel_fun.c'};
fun(1).outdir='$kernels';
fun(1).include='$kernels';

fun(2).source={'$kernels/diagker.c','$kernels/kernel_fun.c'};
fun(2).outdir='$kernels';
fun(2).include='$kernels';

fun(3).source={'$svm/smo1d_mex.c'};
fun(3).outdir='$svm';
fun(3).include='$kernels';

fun(4).source={'$svm/smo_mex.c','$kernels/kernel_fun.c'};
fun(4).outdir='$svm';
fun(4).include='$kernels';

fun(5).source = {'$svm/bsvm2_mex.c','$kernels/kernel_fun.c',...
    '$optimization/gmnplib.c'};
fun(5).outdir = '$svm';
fun(5).include = '$kernels';

fun(6).source = {'$misc/knnclass_mex.c'};
fun(6).outdir = '$misc';
fun(6).include = '$';

fun(7).source = {'$svm/svm2_mex.c','$kernels/kernel_fun.c',...
    '$optimization/gnpplib.c'};
fun(7).outdir = '$svm';
fun(7).include = '$kernels';

fun(8).source={'$kernels/kernelproj_mex.c','$kernels/kernel_fun.c'};
fun(8).outdir='$kernels';
fun(8).include='$kernels';

fun(9).source={'$optimization/gmnp_mex.c','$optimization/gmnplib.c'};
fun(9).outdir='$optimization';
fun(9).include='$';

fun(10).source={'$optimization/gnnls_mex.c','$optimization/gnnlslib.c'};
fun(10).outdir='$optimization';
fun(10).include='$';

fun(11).source={'$optimization/gnpp_mex.c','$optimization/gnpplib.c'};
fun(11).outdir='$optimization';
fun(11).include='$';

fun(12).source={'$optimization/qpssvm_mex.c','$optimization/qpssvmlib.c'};
fun(12).outdir='$optimization';
fun(12).include='$';

fun(13).source={'$optimization/qpbsvm_mex.c','$optimization/qpbsvmlib.c'};
fun(13).outdir='$optimization';
fun(13).include='$';

fun(14).source={'$optimization/gsmo_mex.c','$optimization/gsmolib.c'};
fun(14).outdir='$optimization';
fun(14).include='$';


% -- Compile functions -----------------------------
for i=1:length(fun),
   mexstr = ['mex -O -I' translate(fun(i).include,root) ...
             ' -outdir ' translate(fun(i).outdir, root) ' '];
   
  for j=1:length(fun(i).source),    
    mexstr = [mexstr translate(char(fun(i).source(j)),root) ' '];
  end

  fprintf('%s\n',mexstr);
  
  eval(mexstr);
end

fprintf('MEX-files compiled successfully.\n');

return;

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
else
  p = strrep(p,'$',[toolboxroot '/']);
end
