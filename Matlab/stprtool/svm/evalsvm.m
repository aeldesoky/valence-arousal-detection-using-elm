function [best_model,Errors] = evalsvm(arg1,arg2,arg3)
% EVALSVM Trains and evaluates Support Vector Machines classifier.
%
% Synopsis:
%  [model,Errors] = evalsvm(data,options)
%  [model,Errors] = evalsvm(trn_data,val_data,options)
%
% Description:
%  [model,Errors] = evalsvm(data,options) uses cross-validation
%    to assess SVM classifiers with given kernel arguments and 
%    regularization constants.
%   
%    The kernel type is given in options.ker (see 'help kernel').
%    The SVM solver to be used is specified by field options.solver 
%    (default 'smo'). Both binary and multi-class SVM solvers are 
%    allowed. The input data have the following format:
%    be used with regards to the number of labels of training data 
%      data.X [dim x num_data] ... training vectors.
%      data.y [1 x num_data] ... labels.
%    The set of SVM parameters to be evaluated are specified in:
%      options.arg [dimarg x nargs] ... enumeration of  kernel arguments; 
%        dimarg determins number of kernel argumens (e.g., dimarg = 1 
%        for 'rbf' kernel and dimarg = 2 for 'sigmoid').
%      options.C [1 x nc] ... enumeration of regularization constants.
%
%    Some extra parameters for the selected SVM solver can be
%    specified in the field options.solver_options.
%
%    Each configuration of SVM paramaters is evaluated using the
%    cross-validation. The number of folds is given in 
%    optios.num_folds (default 5). The trained SVM model with 
%    the smallest cross-validation error is returned. The computed
%    cross-validation errors with respect to SVM parametes are 
%    returned in Errors [nc x nargs].
%
%    The progress info is displayed if options.verb is set to 1
%    (default 0).
%
%  [model,Errors] = evalsvm(trn_data,val_data,options) each
%    SVM is trained on the trn_data and evaluated on the 
%    validation val_data instead of using cross-validation.
%   
% Example:
%  trn = load('riply_trn');
%  tst = load('riply_tst');
%  options.ker = 'rbf';
%  options.arg = [0.1 0.5 1 5];
%  options.C = [1 10 100];
%  options.solver = 'smo';
%  options.num_folds = 5;
%  options.verb = 1;
%  [model,Errors] = evalsvm(trn,options);
%  figure; mesh(options.arg,options.C,Errors);
%  hold on; xlabel('arg'); ylabel('C');
%  ypred = svmclass(tst.X,model);
%  cerror(ypred,tst.y)
%
% See also: 
%   SMO, SVMLIGHT, SVMCLASS.
%

% Modifications:
% 17-sep-2004, VF, Help improved. Info about training stage added.
% 18-aug-2004, VF, svm_options changed to solver_options
% 4-june-2004, VF
% 3-jun-2004, VF

if nargin == 2,
  % evaluation by cross-validation
  %=================================================================
  options = c2s(arg2);
  if ~isfield(options,'verb'), options.verb = 0; end
  if ~isfield(options,'solver'), options.solver = 'smo'; end
  if ~isfield(options,'ker'), options.ker = 'linear'; end
  if ~isfield(options,'num_folds'), options.num_folds = 5; end
  if ~isfield(options,'solver_options'), options.solver_options = []; end

  nargs = size(options.arg,2);
  nc = length(options.C);
  [dim,num_data] = size(arg1.X);

  solver_options = options.solver_options;
  solver_options.ker = options.ker;
  min_error = inf;
  Errors = [];
  cnt_model = 0;
  num_model = nc*nargs;
  
  % data partitioning
  [itrn,itst] = crossval(num_data,options.num_folds);
  
  for i = 1:nargs,
    arg = options.arg(:,i);
    for j = 1:nc,
      cnt_model = cnt_model + 1;
      C = options.C(j);
      
      % svm parameters
      solver_options.C = C;
      solver_options.arg = arg;
    
      % display info
      if options.verb == 1,
        fprintf('Model %d/%d: ker=%s, C=%f, arg=', ...
            cnt_model, num_model,solver_options.ker, C);
        fprintf('%f ', arg);
        fprintf('\n');
      end
      
      fold_error = zeros(options.num_folds,1);
      
      for k=1:options.num_folds,

        if options.verb == 1,
          fprintf('fold %d/%d: #trn/tst = %d/%d, training', ...
              k, options.num_folds, length(itrn{k}), length(itst{k}));
        end
        
        trn.X = arg1.X(:,itrn{k});
        trn.y = arg1.y(:,itrn{k});
        tst.X = arg1.X(:,itst{k});
        tst.y = arg1.y(:,itst{k});
        
        % run solver
        model = feval( options.solver, trn, solver_options );
        
        if options.verb == 1,
          ypred = feval( model.fun, trn.X, model );
          fprintf(' err = %.4f, testing', cerror( ypred,trn.y));
        end

        % classify validation data
        ypred = feval( model.fun, tst.X, model );
        
        fold_error(k) = cerror( ypred, tst.y );

        if options.verb == 1,
          fprintf(' err = %.4f\n', fold_error(k));
        end
        
      end
    
      % cross-validation error
      err = mean( fold_error );
      Errors(j,i) =  err;

      if min_error > err,
        min_error = err;
        best_model = model;
        
        if options.verb == 1,
          fprintf('cross-validation error = %.4f (best so far)\n\n', err);
        end
      else
        if options.verb == 1,
          fprintf('cross-validation error = %.4f\n\n', err);
        end
      end
      
    end
  end
  
  % disp info
  if options.verb == 1,
    fprintf('best model: ker=%s, C=%f, arg=', ...
        best_model.options.ker, best_model.options.C);
    fprintf('%f ', best_model.options.arg);
    fprintf('\ncross-validation error = %.4f\n', min_error);
  end
  
elseif nargin == 3, 
  % evaluation using val_data
  %=================================================================  
  options = c2s(arg3);
  if ~isfield(options,'verb'), options.verb = 0; end
  if ~isfield(options,'ker'), options.ker = 'linear'; end
  if ~isfield(options,'solver'), options.solver = 'smo'; end
  if ~isfield(options,'solver_options'), options.solver_options = []; end

  nargs = size(options.arg,2);
  nc = length(options.C);
  [dim,num_data] = size(arg1.X);

  solver_options = options.solver_options;
  solver_options.ker = options.ker;
  min_error = inf;
  Errors = [];
  cnt_model = 0;
  num_model = nc*nargs;
  
  for i = 1:nargs,
    arg = options.arg(:,i);
    for j = 1:nc,
      cnt_model = cnt_model + 1;
      C = options.C(j);
      
      % svm parameters
      solver_options.C = C;
      solver_options.arg = arg;
    
      % display info
      if options.verb == 1,
        fprintf('Model %d/%d: ker=%s, C=%f, arg=', ...
            cnt_model, num_model,solver_options.ker, C);
        fprintf('%f ', arg);
        fprintf('\n');
      end

      if options.verb == 1, 
          fprintf('#trn/tst = %d/%d, training', ...
              length(arg1.y), length(arg2.y));
      end
      
      % run SVM solver
      model = feval( options.solver, arg1, solver_options );
        
      if options.verb == 1,
        ypred = feval( model.fun, arg1.X, model );
        fprintf(' err = %.4f, testing', cerror( ypred,arg1.y));
      end
        
      % classify validation data
      ypred = feval( model.fun, arg2.X, model );

      err = cerror(ypred,arg2.y);
      Errors(j,i) =  err;

      if min_error > err,
        min_error = err;
        best_model = model;
        
        if options.verb == 1,
          fprintf(' err = %.4f (best so far)\n\n', err);
        end
      else
        if options.verb == 1,
          fprintf(' err = %.4f\n\n', err);
        end
      end
    end
  end

  % disp info
  if options.verb == 1,
    fprintf('best model: ker=%s, C=%f, arg=', ...
        best_model.options.ker, best_model.options.C);
    fprintf('%f ', best_model.options.arg);
    fprintf('\ntesting error = %.4f\n', min_error);
  end
  
end

%EOF