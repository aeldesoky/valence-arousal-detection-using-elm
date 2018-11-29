function red_model = redquadh(model)
% REDQUADH reduced SVM classifier with homogeneous quadratic kernel.
%
% Synopsis:
%  red_model = redquadh(model)
%
% Description:
%  It uses reduced set techique (Burges) to compute 
%  simpler SVM binary rule with homogeneous quadratic kernel (x'*y)^2.
%  
% Input:
%  model.Alpha [nsv x 1] Weights of kernel expansion.
%  model.b [scalar] Bias.
%  model.sv.X [dim x nsv] Support vectors.
%  model.options.ker = 'poly'
%  model.options.arg = [2 0]
%
% Output:
%  red_model.Alpha [new_nsv x 1] New weights.
%  red_model.b [scalar] Bias.
%  red_model.sv.X [dim x new_nsv] New "support vectors".
%  ...
%
% Example:
%  trn = load('riply_trn');
%  model = smo(trn,{'ker','poly','arg',[2 0],'C',10});
%  red_model = redquadh( model );
%  figure; ppatterns(trn); psvm(model);
%  figure; ppatterns(trn); psvm(red_model);
%

% Modifications:
% 28-nov-2003, VF

dim=size(model.sv.X,1);
nsv = model.nsv;

S = zeros(dim,dim);

for i=1:dim,
  for j=i:dim,
    S(i,j) = (model.sv.X(i,:).*model.sv.X(j,:) )*model.Alpha(:);
    S(j,i) = S(i,j);
  end
end

[V,D] = eig(S);
D=real(diag(D));
[dummy,inx] = sort(-abs(D));
D=D(inx);
V=V(:,inx);

inx = find(D ~= 0);

red_model.nsv = length(inx);
red_model.Alpha = zeros(red_model.nsv,1);
red_model.b = model.b;
red_model.sv.X = zeros(dim,red_model.nsv);
red_model.options = model.options;
red_model.classifier = 'svmclass';
red_model.eigval = D(inx);

cnt = 0;
for i=inx(:)',
  cnt = cnt+1;
  red_model.sv.X(:,cnt) = V(:,i);
  red_model.Alpha(cnt) = D(i)/(red_model.sv.X(:,cnt)'*red_model.sv.X(:,cnt));
end

return;


