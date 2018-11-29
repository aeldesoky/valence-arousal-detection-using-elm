function model = gda(data,options)
% GDA Generalized Discriminant Analysis.
% 
% Synopsis:
%  model = gda(data)
%  model = gda(data,options)
% 
% Description:
%  This function is implimentation of the Generalized Discriminant
%  Analysis (GDA) [Baudat01]. The GDA is kernelized version of
%  the Linear Discriminant Analysis (LDA). It produce the kernel data
%  projection which increases class separability of the projected 
%  training data.
%
% Input:
%  data [struct] Labeled training data:
%   .X [dim x num_data] Training vectors.
%   .y [1 x num_data] Labels (1,2,..,mclass).
%  
%  options [struct] Defines kernel and a output dimension:
%   .ker [string] Kernel identifier (default 'linear'); 
%     see 'help kernel' for more info.
%   .arg [1 x nargs] Kernel arguments (default 1).
%   .new_dim [1x1] Output dimension (default dim).
%
% Output:
%  model [struct] Kernel projection:
%   .Alpha [num_data x new_dim] Multipliers.
%   .b [new_dim x 1] Bias.
%   .sv.X [dim x num_data] Training data.
%   .options [struct] Copy of used options.
%   .rankK [int] Rank of centered kernel matrix.
%   .nsv [int] Number of training data.
%
% Example:
%  in_data = load('iris');
%  model = gda(in_data,struct('ker','rbf','arg',1));
%  out_data = kernelproj( in_data, model );
%  figure; ppatterns( out_data );
%
% See also 
%  KERNELPROJ, KPCA.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 24-may-2004, VF
% 4-may-2004, VF


% process input arguments
%-----------------------------

% allos data to be given as a cell
data=c2s(data);

% get dimensions
[dim,num_data]=size(data.X);
nclass = max(data.y);

if nargin < 2, options=[]; else options=c2s(options); end
if ~isfield(options, 'ker'), options.ker = 'linear'; end
if ~isfield(options, 'arg'), options.arg = 1; end
if ~isfield(options, 'new_dim'), options.new_dim = dim; end

% sort data according to labels
[tmp,inx] = sort(data.y);
data.y=data.y(inx);
data.X=data.X(:,inx);

% kernel matrix
K = kernel( data.X, options.ker, options.arg );

% centering matrix
J=ones(num_data,num_data)/num_data;
JK = J*K;

% centering data in non-linear space
Kc = K - JK' - JK + JK*J;

% Kc decomposition; Kc = P*Gamma*P'
[P, Gamma]=eig( Kc );
Gamma=diag(Gamma);
[tmp,inx]=sort(Gamma); % sort eigenvalues in ascending order
inx=inx([num_data:-1:1]); % swap indices
Gamma=Gamma(inx);
P=P(:,inx);

% removes eigenvectors with small value
minEigv = Gamma(1,1)/1000;
inx = find( Gamma >= minEigv );
P=P(:,inx);
Gamma=Gamma(inx);
rankKc = length(inx);

Kc = P*diag(Gamma)*P';

% make diagonal block matrix W
W=[];
for i=1:nclass,
  num_data_class=length(find(data.y==i));
  W=blkdiag(W,ones(num_data_class)/num_data_class);
end  

% new dimension of data
model.new_dim=min([options.new_dim, rankKc, nclass-1]);

% compute vector alpha and its normalization 
[Beta, Lambda] = eig( P'*W*P );
Lambda=diag(Lambda);
[tmp,inx]=sort(Lambda);  % sort eigenvalues in ascending order
inx=inx([length(Lambda):-1:1]);  % swap indices
Lambda=Lambda(inx);
Beta=Beta(:,inx(1:model.new_dim));

%model.Alpha=P*inv(diag(Gamma))*Beta;
model.Alpha=P*diag(1./Gamma)*Beta;

% normalization of vectors Alpha
for i=1:model.new_dim,
  model.Alpha(:,i) = model.Alpha(:,i)/...
      sqrt(model.Alpha(:,i)'* Kc * model.Alpha(:,i));
end

% centering Alpha and computing Bias
sumK=sum(K);
model.b=(-sumK*model.Alpha/num_data+...
     sum(model.Alpha)*sum(sumK)/num_data^2)'; 

for i=1:size(model.Alpha,2),
  model.Alpha(:,i) = model.Alpha(:,i)-sum(model.Alpha(:,i))/num_data;
end

% fill model
model.options = options;
model.sv = data;
model.rankK = rankKc;
model.nsv = num_data;
model.fun = 'kernelproj';

return;
