function model=kfd(data,options)
% KFD Kernel Fisher Discriminat.
%
% Synopsis:
%  model = kfd( data )
%  model = kfd( data, options )
%
% Description:
%  This function is an implementation of the Kernel Fisher
%  Discriminant (KFD) [Mika99a]. The aim is to find a binary 
%  kernel classifier which is the linear decision function in a 
%  feature space induced by the selected kernel function. 
%  The bias is found decision function is trainined by the 
%  linear SVM on the data projected on the optimal direction.
%
% Input:
%  data [struct] Training binary labeled data:
%   .X [dim x num_data] Vectors.
%   .y [1 x num_data] Labels (1 or 2).
%
%  options [struct] Control parameters:
%   .ker [string] Kernel identifier (default 'linear'). 
%     See 'help kernel' for more info.
%   .arg [1 x nargs] Kernel argument(s).
%   .C [1x1] Regularization constant of the linear 1-D SVM 
%     used to optimize the bias (default C=inf).
%   .mu [1x1] Regularization constant added to the diagonal of 
%     the within scatter matrix (default 1e-4).
% 
% Output:
%  model [struct] Binary SVM classifier:
%   .Alpha [num_data x 1] Weight vector.
%   .b [1x1] Bias of decision function.
%   .sv.X [dim x num_data] Training data (support vectors).
%
%   .trnerr [1x1] Training classification error.
%   .kercnt [1x1] Number of kernel evaluations used during training.
%   .nsv [1x1] Number of support vectors.
%   .options [struct] Copy of options.
%   .cputime [1x1] Used cputime.
%
% Example:
%  trn = load('riply_trn');
%  options = struct('ker','rbf','arg',1,'C',10,'mu',0.001);
%  model = kfd(trn, options)
%  figure; ppatterns(trn); psvm(model);
%
% See also 
%  SVMCLASS, FLD, SVM.
%

% Modifications:
% 17-may-2004, VF
% 14-may-2004, VF
% 7-july-2003, VF

% timer
tic;

% processing inputs
% ======================================
[dim,num_data]=size(data.X);
if nargin < 2, options=[]; else options=c2s(options); end
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'mu'), options.mu = 1e-4; end
if ~isfield(options,'C'), options.C = inf; end


% creates matrices M and N 
%=================================
inx1=find(data.y==1);
inx2=find(data.y==2);
l1=length(inx1);
l2=length(inx2);

K = kernel(data.X,options.ker,options.arg);

M1=sum(K(:,inx1),2)/l1;
M2=sum(K(:,inx2),2)/l2;
M=(M1-M2)*(M1-M2)';

E1=eye(l1,l1);
E2=eye(l2,l2);
J1=ones(l1,l1)/l1;
J2=ones(l2,l2)/l2;

N = K(:,inx1)*(E1-J1)*K(:,inx1)' + K(:,inx2)*(E2-J2)*K(:,inx2)';

% regularization
N = N + options.mu * eye(num_data,num_data);

% Optimization
%==============================
%%[Alpha,V,U] = svds( inv(N)*M,1);
Alpha=inv(N)*(M1-M2); % It yields the same Alpha up to scale

% project data on the found direction
projx=(K*Alpha)';

% training bias of decision rule
lin_model = svm1d(struct('X',projx,'y',data.y),struct('C',options.C));

% fill output structure
%===============================
model.Alpha = lin_model.W*Alpha(:);
model.b = lin_model.b;
model.sv = data;
if strcmp(options.ker,'linear'),
  % in the linar case compute normal vector explicitely
  model.W = model.sv.X*model.Alpha;
end  
model.trnerr = lin_model.trnerr;
model.nsv = num_data;
model.kercnt = num_data*(num_data+1)/2;
model.options = options;
model.fun = 'svmclass';
model.cputime=toc;

return;