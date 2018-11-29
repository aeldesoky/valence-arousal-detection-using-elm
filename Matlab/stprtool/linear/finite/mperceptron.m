function model = mpeceptron( data, options, init_model )
% MPERCEPTRON Perceptron algorithm to train linear machine.
%
% Synopsis:
%  model = mpeceptron(data)
%  model = mpeceptron(data,options)
%  model = mpeceptron(data,options,init_model)
%
% Description:
%  model = mperceptron(data) uses the Perceptron learning rule
%   to train linear machine (multi-class linear classitfier).
%   The multi-class problem is transformed to the single-class
%   one using the Kessler's construction [DHS01][SH10].
%
%  model = mperceptron(data,options) specifies stopping condition of
%   the algorithm in structure options:
%    .tmax [1x1]... maximal number of iterations.
%
%  model = mperceptron(data,options,init_model) specifies initial 
%   model which must contain:
%    .W [dim x nclass] ... Normal vectors.
%    .b [nclass x 1] ... Biases.
%
% Input:
%  data [struct] Labeled training data:
%   .X [dim x num_data] Training vectors.
%   .y [1 x num_data] Labels (1,2,...,nclass).
%
%  options [struct] 
%   .tmax [1x1] Maximal number of iterations (default tmax=inf).
%  
%  init_model [struct] Initial model; must contain items .W, .b.
%
% Output:
%  model [struct] Multi-class linear classifier:
%   .W [dim x nclass] Normal vectors.
%   .b [nclass x 1] Biases.
%
%   .exitflag [1x1] 1 ... perceptron has converged.
%                   0 ... number of iterations exceeded tmax.
%   .t [1x1] Number of iterations.
%
% Example:
%  data = load('pentagon');
%  model = mperceptron( data );
%  figure; ppatterns( data ); pboundary( model );
%
% See also 
%  PERCEPTRON, LINCLASS, EKOZINEC.
%

% Modifications:
% 21-may-2004, VF
% 18-may-2004, VF

% input arguments
%----------------------------------------
[dim,num_data] = size(data.X);
nclass = max(data.y);

if nargin < 2, options = []; else options = c2s(options); end
if ~isfield(options,'tmax'), options.tmax = inf; end 

if nargin == 3,
  model = init_model;
else
  model.W = zeros(dim,nclass);
  model.b = zeros(nclass,1);
end
model.t = 0;

% main loop
% -----------------------------------
model.exitflag = 0;
while options.tmax > model.t & model.exitflag == 0,

  model.t = model.t+1;

  model.exitflag = 1;
  
  % search for misclassified vector
  for i=1:nclass,

    class_i = find( data.y == i);
    dfce_i = model.W(:,i)'*data.X(:,class_i) + model.b(i);

    for j=setdiff([1:nclass], i),
     
      dfce_j = model.W(:,j)'*data.X(:,class_i) + model.b(j);
      
      [min_diff,inx] = min( dfce_i - dfce_j);
      
      if min_diff <= 0,
        % Perceptron rule
        
        % take index of misclassified vector
        inx=class_i(inx(1));

        model.W(:,i) = model.W(:,i) + data.X(:,inx);
        model.b(i) = model.b(i) + 1;

        model.W(:,j) = model.W(:,j) - data.X(:,inx);
        model.b(j) = model.b(j) - 1;
        
        % error was found
        model.exitflag = 0;
        break;
      end
      
    end
    if model.exitflag == 0, break; end
  end
  
end

model.fun = 'linclass';

return;
% EOF
