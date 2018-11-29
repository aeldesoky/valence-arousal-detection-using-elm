function model = oaosvm(data,options)
% OAOSVM Multi-class SVM using One-Against-One decomposition.
% 
% Synopsis:
%  model = oaosvm( data )
%  model = oaosvm( data, options )
%
% Description:
%  model = oaosvm( data ) uses one-agains-one deconposition
%   to train the multi-class Support Vector Machines (SVM)
%   classifier. The classification into nclass classes 
%   is decomposed into nrule = (nclass-1)*nclass/2 binary 
%   problems.
%
%  model = oaosvm( data, options) allows to specify the
%   binary SVM solver and its paramaters.
%
% Input:
%  data [struct] Training data:
%   .X [dim x num_data] Training vectors.
%   .y [1 x num_data] Labels of training data (1,2,...,nclass). 
%
%  options [struct] Control parameters:
%   .bin_svm [string] Function which implements the binary SVM 
%     solver; (default 'smo').
%   .verb [1x1] If 1 then a progress info is displayed (default 0).
%  The other fields of options specifies the options of the binary
%  solver (e.g., ker, arg, C). See help of the selected solver.
%
% Output:
%  model [struct] Multi-class SVM majority voting classifier:
%   .Alpha [nsv x nrule] Weights (Lagrangeans).
%   .bin_y [2 x nrule] Translation between binary responses of
%     the discriminant functions and class labels.
%   .b [nrule x 1] Biases of discriminant functions.
%   .sv.X [dim x nsv] Support vectors.
%   .nsv [1x1] Number of support vectors.
%   .trnerr [1x1] Training error.
%   .kercnt [1x1] Number of kernel evaluations.
%   .options [struct[ Copy of input argument options.
%
% Example:
%  data = load('pentagon');
%  options = struct('ker','rbf','arg',1,'C',1000,'verb',1);
%  model = oaosvm( data, options );
%  figure; 
%  ppatterns(data); ppatterns(model.sv.X,'ok',13);
%  pboundary( model );
%  
% See also 
%  MVSVMCLASS, OAASVM.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2005, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 25-jan-2005, VF, option solver replaced by bin_svm 
% 26-may-2004, VF
% 4-feb-2004, VF
% 9-Feb-2003, VF

% Process inputs
%-----------------------------
if nargin < 2, options = []; else options=c2s(options); end
if ~isfield(options,'verb'), options.verb = 0; end
if ~isfield(options,'bin_svm'), options.bin_svm = 'smo'; end
if ~isfield(options,'ker'), options.ker = 'linear'; end
if ~isfield(options,'arg'), options.arg = 1; end
if ~isfield(options,'C'), options.C = inf; end

[dim,num_data] = size(data.X);
nclass = max(data.y);
nrule = (nclass-1)*nclass/2;

% display info
%---------------------
if options.verb == 1,
  fprintf('Binary rules: %d\n', nrule);
  fprintf('Training data: %d\n', num_data);
  fprintf('Dimension: %d \n', dim);
  if isfield( options, 'ker'), fprintf('Kernel: %s\n', options.ker); end
  if isfield( options, 'arg'), fprintf('arg: %f\n', options.arg(1)); end
  if isfield( options, 'C'), fprintf('C: %f\n', options.C); end
end

%----------------------------------------
Alpha = zeros(num_data,nrule);
b = zeros(nrule,1);
bin_y = zeros(2,nrule);
kercnt = 0;

% One-Against-One decomposition
%-----------------------------------
rule = 0;
for class1 = 1:nclass-1,
  for class2 = class1+1:nclass,
  
    rule = rule + 1;
    
    if options.verb == 1,
      fprintf('building rule %d-%d (%d of %d)', ...
        class1, class2, rule, nrule );
    end
    
    % set binary subtask
    %---------------------------------------------
    bin_y(1,rule) = class1;
    bin_y(2,rule) = class2;
    data_inx = find(data.y==class1 | data.y==class2);
    bin_data.X = data.X(:, data_inx);
    bin_data.y = data.y(data_inx);
    bin_data.y(find(bin_data.y == class1)) = 1;
    bin_data.y(find(bin_data.y == class2)) = 2;
    
    % solve binary subtask
    %---------------------------------------------
    bin_model = feval( options.bin_svm, bin_data, options );

    Alpha(data_inx(bin_model.sv.inx),rule) = bin_model.Alpha(:);
    b(rule) = bin_model.b;

    kercnt = kercnt + bin_model.kercnt;
    
    % progress info
    %-----------------------------
    if options.verb ==1,
     if isfield(bin_model, 'trnerr'),
       fprintf(': trnerr = %.4f', bin_model.trnerr);
     end
     if isfield(bin_model, 'margin'),
       fprintf(', margin = %f', bin_model.margin );
     end
     fprintf('\n');
    end
  
  end
end

% set output model
%---------------------------------

% indices of all support vectors
inx = find(sum(abs(Alpha),2)~= 0);

model.Alpha = Alpha(inx,:);
model.b = b;
model.bin_y = bin_y;
model.sv.X = data.X(:,inx);
model.sv.y = data.y(inx);
model.sv.inx = inx;
model.nsv = length(inx);
model.kercnt = kercnt;
model.options = options;
model.fun = 'mvsvmclass';
model.trnerr = cerror( mvsvmclass(data.X, model), data.y );

% display info
%--------------------
if options.verb == 1,
  fprintf('Total training error = %.4f\n', model.trnerr);
end

return;
% EOF



