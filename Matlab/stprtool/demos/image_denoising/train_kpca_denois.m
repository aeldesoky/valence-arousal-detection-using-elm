% TRAIN_KPCA_DENOIS Training of kernel PCA model for image denoising. 
%
% Description:
%  The kernel PCA model is trained to describe an input
%  class of images corrupted by noise [Mika99b]. The training 
%  data contains images corrupted by noise and corresponding 
%  ground truth. The free paramaters of the kernel PCA
%  are tuned by cross-validation. The objective function 
%  is a sum of squared differences between ground truth 
%  images and reconstructed images. The greedy KPCA algorithm 
%  is used to train the kernel PCA model.
%
% See also
%  GREEDYKPCA, KPCAREC, KPCA.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 07-jun-2004, VF
% 06-jun-2004, VF
% 17-mar-2004, VF

% Setting
% -------------------------------------

options.ker = 'rbf';   % kernel
options.m = 500;       % # of basis vectors
options.p = 10;        % deth of search for the best basis vector
options.verb = 1;

% # folds for cross-validation; 
% num_folds = 1 means 50/50 - training/testing part
num_folds = 1;  

% algorithm to compute kernel PCA
%KPCA_Algo = 'kpca';
KPCA_Algo = 'greedykpca';

% parameters to be evaluated by cross-validation:
%New_Dim_Range = [50 100 200 300]; % usps
%Arg_Range = [3.5 4 5 6 7 8];      % usps

New_Dim_Range = [1 2];   % noisy_circle
Arg_Range = [0.5 1 2 3]; % noisy_circle

% input/output files
input_data_file = 'noisy_circle';
output_data_file = [];
%input_data_file = '/home.dokt/xfrancv/data/usps/usps_noisy';
%output_data_file = 'USPSModelGreedyKPCA';

% Loads training and testing data.
% -------------------------------------
load(input_data_file,'trn','tst');
[dim,num_data] = size(trn.X);

% Data partitioning for cross-validation
[itrn,itst] = crossval(num_data,num_folds);

% Tuning kernel PCA model
% -------------------------------------
Mse = [];

for arg = Arg_Range,
  for new_dim = New_Dim_Range,
    
    fprintf('\nnew_dim = %d, arg = %f\n', new_dim, arg);
   
    cv_mse = 0;  
    for i=1:num_folds,
    
      fprintf('\n');

      % training and validation part of data
      trn_X = trn.gnd_X(:,itrn{i});
      val_gnd_X = trn.gnd_X(:,itst{i});
      val_corr_X = trn.X(:,itst{i});
    
      fprintf('Computing Kernel PCA...');
      options.arg = arg;
      options.new_dim = new_dim;
      kpca_model = feval( KPCA_Algo, trn_X, options);
      fprintf('done.\n');

      % data restoration
      val_reconst_X = kpcarec(val_corr_X, kpca_model);
  
      % compute error
      dummy = (val_reconst_X - val_gnd_X).^2;
      mse = sum(dummy(:))/size(val_gnd_X,2);
  
      fprintf('folder %d/%d: validation errors mse=%f\n', ...
        i, num_folds, mse);
   
      cv_mse = cv_mse + mse;
    end

    % compute cross-validation error
    cv_mse = cv_mse/num_folds;
 
    Mse(find(new_dim==New_Dim_Range),find(arg==Arg_Range)) = cv_mse;
 
    fprintf('Kernel arg = %f: mse = %f\n', arg, cv_mse);
  end
end

% take the best parameters
%----------------------------------------------
[inx1,inx2] = find(Mse==min(Mse(:)));
fprintf('\nMin(mse) = %f, dim = %f, arg = %f\n', ...
   Mse(inx1,inx2), New_Dim_Range(inx1), Arg_Range(inx2) );

% compute kernel PCA model with best parameters
% using all training data
%---------------------------------------------
fprintf('Computing optimal Kernel PCA...');
options.arg = Arg_Range(inx2);
options.new_dim = New_Dim_Range(inx1);
kpca_model = feval( KPCA_Algo, trn.X, options);
fprintf('done.\n');

if isempty(output_data_file),
  % plot results of tuning
  figure; hold on;
  xlabel('\sigma'); ylabel('mse');

  h = [];
  clear Str;
  for i=1:length(New_Dim_Range),
    h = [h, plot(Arg_Range, Mse(i,:),marker_color(i) )];
    Str{i} = sprintf('dim = %d', New_Dim_Range(i));
  end

  legend(h,Str);
else
  % save model to file
  save(output_data_file,'Arg_Range','New_Dim_Range',...
     'options','Mse','num_folds','input_data_file',...
     'output_data_file','KPCA_Algo','kpca_model');
end

% plot denosing in 2D case only
%-------------------------------------
if dim == 2 & isempty(output_data_file),

  X = kpcarec(tst.X,kpca_model);

  mse = sum(sum((X-tst.gnd_X).^2 ));
  fprintf('\ntest mse=%f\n', mse);

  figure; hold on;
  h0=ppatterns(tst.gnd_X,'r+');
  h1=ppatterns(tst.X,'gx');
  h2=ppatterns(X,'bo');
  legend([h0 h1 h2],'Ground truth','Noisy','Reconst');
end

% EOF
