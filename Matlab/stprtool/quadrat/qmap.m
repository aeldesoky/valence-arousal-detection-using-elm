function Y=qmap(X)
% QMAP Quadratic data mapping.
%
% Synopsis:
%  Y = qmap(X)
%  out_data = qmap(in_data)
%
% Description:
%  Y = qmap(X) mapps input n-dimensional data X into a new 
%   (n*(n+3)/2)-dimensional space using the quadratic mapping. An 
%   input vector x=X(:,i) is mapped to its image y=Y(:,i) such that
%
%    y = [x(1)     ;x(2)     ;...;x(n);
%         x(1)*x(1);x(1)*x(2);...;x(1)*x(n);
%                   x(2)*x(2);...;x(2)*x(n);
%                             ...
%                                 x(n)*x(n)]
%
%  out_data = qmap(in_data) use to map labeled data such that
%    out_data.X = qmap(in_data.X)
%    out_data.y = in_data.y
%
% Example:
%  orig_data = load('vltava');
%  map_data = qmap(orig_data);
%  lin_model = perceptron(map_data);
%  quad_model = lin2quad(lin_model);
%  figure; ppatterns(orig_data); 
%  pboundary(quad_model);
%
% See also 
%  LIN2QUAD, QUADCLASS.
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 17-may-2004, VF
% 22-Oct-2003, VF
% 26-June-2001, V.Franc, comments repared.
% 24. 6.00 V. Hlavac, comments polished.

% if labeled data (struct (X,y)) enters the function 
% then the mapping is performed on the vectors and the 
% labels are preserved
if isstruct(X), tmp = X.y; X = X.X; else tmp = []; end

% dimension
[dim,num_data]=size(X);

new_dim=dim*(dim+3)/2;

Y=zeros(new_dim,num_data);

Y(1:dim,:)=X;

inx=dim+1;
for i=1:dim,
   for j=i:dim,
      Y(inx,:)=X(i,:).*X(j,:);
      inx=inx+1;
   end
end

% deals with labeled data
if ~isempty(tmp) Y.X = Y; Y.y = tmp; end
  
return;
% EOF