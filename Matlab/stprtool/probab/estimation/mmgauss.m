function model=mmgauss(X,options,init_model)
% MMGAUSS Minimax estimation of Gaussian distribution.
%
% Synopsis:
%  model = mmgauss(X)
%  model = mmgauss(X,options)
%  model = mmgauss(X,options,init_model)
% 
% Description:
%  This function computes the minimax estimation of Gaussian 
%  parameters. The minimax estimation (reffer to [SH10]) for 
%  Gaussian model is defined as:
%
%   (Mean,Cov) = argmax   min( pdfgauss(X, Mean, Cov) ).
%               Mean,Cov   
%    
%  The sample data X should be good representatives of the
%  distribution. In contrast to maximum-likelihood estimation, 
%  the data do not have to be i.i.d.
%
%  An itrative algorithm is used for estimation. It iterates
%  until 
%     upper_bound - lower_bound < eps,
%  where eps is prescribed precission and upper_bound, lower_bound
%  are bounds on the optimal solution
%   upper_bound >   max   min( pdfgauss(X, Mean, Cov) ) > lower_bound
%                 Mean,Cov   
%
% Input:
%  X [dim x num_data] Data sample.
%  
%  options [struct] Control parameters:
%   .eps [1x1] Precision of found estimate (default 0.1).
%   .tmax [1x1] Maximal number of iterations (default inf).
%   .cov_type [int] Type of estimated covariance matrix:
%     cov_type = 'full'      full covariance matrix (default)
%     cov_type = 'diag'      diagonal covarinace matrix
%     cov_type = 'spherical' spherical covariance matrix
%   .verb [int] If 1 then info is printed (default 0).
%
%  init_model [struct] Initial model:
%   .Alpha [1xnum_data] Weights of training vectors.
%   .t [1x1] (optional) Counter of iterations.
%
% Output:
%  model [struct] Gaussian distribution:
%   .Mean [dim x 1] Estimated mean vector.
%   .Cov [dim x dim] Estimated covariance matrix.
%
%   .t [1x1] Number of iterations.
%   .exitflag [1x1] 1 ... (upper_bound - lower_bound) < eps
%                   0 ... maximal number of iterations tmax exceeded.
%   .upper_bound [1x1] Upper bound on the optimized criterion.
%   .lower_bound [1x1] Lower bound on the optimized criterion.
%   .Alpha [1 x num_data] Data weights. The minimax estimate
%     is equal to maximum-likelihood estimate of weighted data.
%   .options [struct] Copy of used options.
%
% Example: 
%  X = [[0;0] [1;0] [0;1]];
%  mm_model = mmgauss(X);
%  figure; ppatterns(X);
%  pgauss(mm_model, struct('p',exp(mm_model.lower_bound')));
%
% See also 
%  PDFGAUSS, MLCGMM, EMGMM.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 26-may-2004, VF
% 30-apr-2004, VF
% 19-sep-2003, VF
% 27-feb-2003, VF
% 24. 6.00 V. Hlavac, comments polished.

[dim,num_data]=size(X);

% processing input arguments
% ------------------------------------------
if nargin < 2, options=[]; else options = c2s(options); end
if ~isfield(options,'eps'), options.eps =0.1; end
if ~isfield(options,'tmax'), options.tmax = inf; end
if ~isfield(options,'verb'), options.verb = 0; end
if ~isfield(options,'cov_type'), options.cov_type = 'full'; end

% inicialization
%---------------------------------
if nargin < 3,
  model.Alpha = ones(1,num_data);  
  model.t = 0;
  model.fun = 'pdfgauss';
  model.options = options;
else
  model = init_model;
  if ~isfield(init_model,'t'), model.t = 0; end
end

% Main loop 
% ----------------------------------------
stop = 0;
while ~stop & options.tmax > model.t,

   if options.verb,
     fprintf('iteration %d: ', model.t );
   end

   % compute ML estimate for given weights model.Alpha 
   tmp_model = melgmm( X, model.Alpha, options.cov_type);
    
   model.Mean = tmp_model.Mean;
   model.Cov = tmp_model.Cov;
   
   % find a sample with the minimal probability
   logPx = log( pdfgauss(X, model));
   [minLogPx,min_inx] = min( logPx );
   
   % compute upper bound and lower bound
   model.upper_bound=sum(model.Alpha.*logPx)/sum(model.Alpha);
   model.lower_bound=minLogPx;

   if options.verb,
    fprintf('upper_bound=%f, lower_bound=%f\n', model.upper_bound, ...
     model.lower_bound );
   end

   % check stopping condition
   if model.upper_bound - model.lower_bound < options.eps, 
     stop = 1;
     model.exitflag = 1;
   else
     % increase occurance of the 'worst' sample by 1
     model.Alpha(min_inx) = model.Alpha(min_inx) + 1;
     model.t = model.t + 1;
     model.exitflag = 0;
   end
   
end

return;
