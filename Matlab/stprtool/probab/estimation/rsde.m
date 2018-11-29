function model = rsde(X,options)
% RSDE Reduced Set Density Estimator.
%
% Synopsis:
%  model = rsde(X,options)
%
% Description:
%  This function implements the Reduced Set Density Estimator 
%  [Girol03] which provides kernel density estimate optimal 
%  in the L2 sense. The density is modeled as the weighted sum 
%  of Gaussians (RBF kernel) centered in selected subset of 
%  training data. 
%
%  The estimation is expressed as a special instance of the
%  Quadratic Programming task (see 'help gmnp').
%  
% Input:
%  X [dim x num_data] Input data sample.
%  options [struct] Control parameters:
%   .arg [1x1] Standard deviation of the Gaussian kernel.
%   .solver [string] QP solver (see 'help gmnp'); 'imdm' default.
%
% Output:
%  model [struct] Output density model:
%   .Alpha [nsv x 1] Weights of the kernel functions.
%   .sv.X [dim x nsv] Selected centers of kernel functions.
%   .nsv [1x1] Number of selected centers.
%   .options.arg = options.arg.
%   .options.ker = 'rbf'
%   .stat [struct] Statistics about optimization:
%     .access [1x1] Number of requested columns of matrix H.
%     .t [1x1] Number of iterations.
%     .UB [1x1] Upper bound on the optimal value of criterion. 
%     .LB [1x1] Lower bound on the optimal value of criterion. 
%     .LB_History [1x(t+1)] LB with respect to iteration.
%     .UB_History [1x(t+1)] UB with respect to iteration.
%     .NA [1x1] Number of non-zero entries in solution.
%   
% Example:
%  gnd = struct('Mean',[-2 3],'Cov',[1 0.5],'Prior',[0.4 0.6]);
%  sample = gmmsamp( gnd, 1000 );
%  figure; hold on; ppatterns(sample.X);
%  plot([-4:0.1:8], pdfgmm([-4:0.1:8],gnd),'r');
%
%  model = rsde(sample.X,struct('arg',0.7));
%  x = linspace(-4,8,100);
%  plot(x,kernelproj(x,model),'g'); 
%  ppatterns(model.sv.X,'ob',13);
%  Reduction = model.nsv/size(sample.X,2)
%
% See also 
%  KERNELPROJ, EMGMM, MLCGMM, GMNP.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2005, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 24-jan-2005, VF, Fast QP solver (GMNP) was used instead of QUADPROG.
% 17-sep-2004, VF, revised


% Input arguments
%-------------------------------------------------------
[dim,num_data] = size(X);
if nargin < 2, options = []; else options = c2s(options); end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'solver'), options.solver = 'imdm'; end

options.ker = 'rbf';

% Solve associted QP task
%-------------------------------------------------------

%G2h = 1/((2*pi)^(dim/2)*(sqrt(2)*options.arg)^dim)*...
%    kernel(X,options.ker,options.arg*sqrt(2));
%Ph = -1/((2*pi)^(dim/2)*options.arg^dim)*...
%    sum(kernel(X,options.ker,options.arg),2)/num_data;
G2h = kernel(X,options.ker,options.arg*sqrt(2))/sqrt(2);
Ph = -sum(kernel(X,options.ker,options.arg),2)/num_data;

[Alpha,fval,stat]= gmnp(G2h,Ph,options);
 
inx = find(Alpha > 0);
model.Alpha = Alpha(inx)/((2*pi)^(dim/2)*options.arg^dim);
model.b = 0;
model.sv.X = X(:,inx);
model.nsv = length(inx);
model.pr2 = model.Alpha'*G2h(inx,inx)*model.Alpha;
model.options = options;
model.stat = stat;
model.fun = 'kernelproj';

return;
% EOF
