function red_model = rspoly2(model,max_nsv)
% RSPOLY2 Reduced set method for second order homogeneous polynomial kernel.
%
% Synopsis:
%  red_model = rspoly2(model)
%  red_model = rspoly2(model,max_nsv)
%
% Description:
%  It uses reduced set techique to reduce complexity
%  of the kernel expansion with second order homogeneous polynomial 
%  kernel k(x,y) = (x'*y)^2 = kernel(x,y,'poly',2) .
%
%  The method was published in 
%  J.C.Burges: Simplified Support Vector Decision Rules. ICML, 1996.
%  
% Input:
%  model [struct] Kernel expansion:
%   .Alpha [nsv x 1] Weights of kernel expansion.
%   .b [1x1] Bias.
%   .sv.X [dim x nsv] Support vectors.
%   .options.ker = 'poly'
%   .options.arg = [2 0]
%
%  max_nsv [1x1] Maximal number of new support vectors. If not given 
%   then the new expansion approximates the original one exactly with
%   at most dim support vectors. 
%
% Output:
%  red_model [struct] Reduced kernel expansion:
%  red_model.Alpha [new_nsv x 1] New weights.
%  red_model.b [scalar] Bias.
%  red_model.sv.X [dim x new_nsv] New support vectors.
%  ...
%
% Example:
%  trn = load('riply_trn');
%  model = smo(trn,struct('ker','poly','arg',[2 0],'C',10));
%  red_model = rspoly2( model );
%  figure; 
%  subplot(1,2,1); axis square; ppatterns(trn); psvm(model);
%  subplot(1,2,2); axis square; ppatterns(trn); psvm(red_model);
%  
% See also 
%  RSRBF.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2004, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 22-dec-2004, VF, header and comments added
% 28-nov-2003, VF

% check inputs
if nargin < 2, max_nsv = inf; end
if strcmpi(model.options.ker,'poly') ~= 1 | ...
   (model.options.arg ~= [2 0] & model.options.arg ~= 2)
  error('Kernel must be homogeneous second order polynomial.');
end

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
D = real(diag(D));
[dummy,inx] = sort(-abs(D));
D=D(inx);
V=V(:,inx);
inx = find(D ~= 0);

% take at most max_nsv support vectors
inx = inx(1:min(max_nsv,length(inx)));

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
% EOF