function model=perceptron(data,options,init_model)
% PERCEPTRON Perceptron algorithm to train binary linear classifier. 
%
% Synopsis:
%  model = perceptron(data)
%  model = perceptron(data,options)
%  model = perceptron(data,options,init_model)
%
% Description:
%  model = perceptron(data) uses the Perceptron learning rule
%   to find separating hyperplane from given binary labeled data.
%
%  model = perceptron(data,options) specifies stopping condition of
%   the algorithm in structure options:
%    .tmax [1x1]... maximal number of iterations.
%
%   If tmax==-1 then it only returns index (model.last_update)
%   of data vector which should be used by the algorithm for updating
%   the linear rule in the next iteration.
%
%  model = perceptron(data,options,init_model) specifies initial model
%   which must contain:
%    .W [dim x 1] ... normal vector.
%    .b [1x1] ... bias of hyperplane.
%    .t [1x1] (optional) ... iteration number.
%
% Input:
%  data [struct] Labeled (binary) training data. 
%   .X [dim x num_data] Input vectors.
%   .y [1 x num_data] Labels (1 or 2).
%
%  options [struct] 
%   .tmax [1x1] Maximal number of iterations (default tmax=inf).
%     If tmax==-1 then it does not perform any iteration but returns only 
%     index of the point which should be used to update linear rule.
%  
%  init_model [struct] Initial model; must contain items
%    .W, .b and .t (see above).
%
% Output:
%  model [struct] Binary linear classifier:
%   .W [dim x 1] Normal vector of hyperplane.
%   .b [1x1] Bias of hyperplane.
%  
%   .exitflag [1x1] 1 ... perceptron has converged.
%                   0 ... number of iterations exceeded tmax.
%   .t [int] Number of iterations.
%   .last_update [d x 1] Index of the last point used for update.
%
% Example:
%  data = genlsdata( 2, 50, 1);
%  model = perceptron(data)
%  figure; ppatterns(data); pline(model); 
%
% See also 
%  EKOZINEC, MPERCEPTRON, LINCLASS.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 17-sep-2003, VF
% 16-Feb-2003, VF
% 20-Jan-2003, VF
% 7-jan-2002, VF. A new coat.
% 24. 6.00 V. Hlavac, comments polished.
% 15-dec-2000, texts, returns bad point

% get data dimensions
[dim,num_data] = size(data.X);

% Process input arguments 
% -----------------------------
if nargin < 2, options = []; else options = c2s(options); end
if ~isfield(options,'tmax'), options.tmax = inf; end

if nargin < 3,
  % create init model
  %----------------------------
  model.W = zeros(dim,1);
  model.b = 0;
else
  % take init model from input
  %----------------------------
  model = init_model;
end
if ~isfield(model,'t'), model.t = 0; end

model.exitflag = 0;
model.last_update = 0;

% Add one constant coordinates to the data and swap
% points from the second class along the origin.
% ----------------------------------------------------
data.X = [data.X; ones(1,num_data)];
dim=dim+1;
inx = find(data.y==2);
data.X(:,inx) = -data.X(:,inx);
W = [model.W;model.b];                                                         

if options.tmax == -1,
  % return index of point which should be used to update linear rule
  %----------------------------------------------------------------------
  fvalue = data.X'*W;
  inx = find( fvalue <= 0 );
  
  if length(inx)==0,
    model.exitflag = 1;
  else
    model.exitflag = 0;
    model.last_update = inx(1);
  end
else

  % main loop 
  % -----------------------------------
  while options.tmax > model.t & model.exitflag == 0,
    model.t = model.t+1;
    
    % Compute discriminant function for all data
    fvalue = data.X'*W;
    inx = find( fvalue <= 0 );
  
    if length(inx)==0,
      model.exitflag = 1;
    else
      inx = inx(1);
      model.exitflag = 0;
      model.last_update = inx;

      % Update model using the Perceptron rule
      W = W + data.X(:,inx);
    end
  end
   
  % separates normal vector and bias
  model.W = W(1:dim-1);
  model.b = W(dim);
end

model.fun = 'linclass';

return;
