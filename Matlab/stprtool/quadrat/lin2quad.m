function quad_model = lin2quad(lin_model)
% LIN2QUAD Merges linear rule and quadratic mapping.
% 
% Synopsis:
%  quad_model = lin2quad(lin_model)
%
% Description:
%  It recomputes parameters of input linear classifier
%  onto parameters of output quadratic classisifier.
%  The linear classifier is assumed to be trained on 
%  the n-dimensional data quadraticaly mapped (see 'help qmap') 
%  into the (n*(n+3)/2)-dimensional feature space. The linear 
%  classifier in the feature space appears as the quadratic 
%  classifier in the original data space.
%
% Description:
%  orig_data = load('vltava');
%  map_data = qmap(orig_data);
%  lin_model = perceptron(map_data);
%  quad_model = lin2quad(lin_model);
%  figure; ppatterns(orig_data); 
%  pboundary(quad_model);
%
% See also 
%  QUACLASS, QMAP
%

% Modifications:
% 09-jun-2004, VF
% 17-may-2004, VF

% allows input model to be cell
lin_model = c2s( lin_model );

% check dimension
[m, nfun] = size( lin_model.W );
n = (-3 + sqrt( 9 + 8*m ))/2;
if round(n) ~= n,
  error('Wrong dimension of input linear classifier.');
end

quad_model.A = zeros(n,n,nfun);
quad_model.B = lin_model.W(1:n,:);
quad_model.C = lin_model.b(:)';

for k=1:nfun

 cnt = n;
 for i=1:n,
   for j=i:n,
     cnt = cnt + 1;
     if i == j,
       quad_model.A(i,j,k) = lin_model.W(cnt,k);
     else
       quad_model.A(i,j,k) = 0.5*lin_model.W(cnt,k);
       quad_model.A(j,i,k) = 0.5*lin_model.W(cnt,k);
     end
   end
 end
end

quad_model.fun = 'quadclass';

return;
%EOF