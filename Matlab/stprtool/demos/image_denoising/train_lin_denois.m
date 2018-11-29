% TRAIN_LIN_DENOIS Training of linear PCA model for image denoising. 
%
% Description:
%  The linear PCA model is trained to describe an input
%  class of images corrupted by noise. The training data 
%  contains images corrupted by noise and corresponding 
%  ground truth. The output dimension of the linear PCA
%  is tuned by cross-validation. The objective function 
%  is a sum of squared differences between ground truth 
%  images and reconstructed images. 
%
% See also
%  PCA, PCAREC, LINPROJ.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 07-jun-2004, VF
% 05-may-2004, VF
% 17-mar-2004, VF

% Setting 
%-------------------------------------------------------

% # folds for cross-validation; 
% num_folds = 1 means 50/50 - training/testing part
num_folds = 1;  

% parameters to be evaluated by cross-validation:
%New_Dim_Range = [5 10 20 30 40 60 80 100]; % USPS
New_Dim_Range = [1 2]; % noisy_circle

input_data_file = 'noisy_circle';
output_data_file = [];
%input_data_file = '/home.dokt/xfrancv/data/usps/usps_noisy';
%output_data_file = 'USPSModelLinPCA';

%-------------------------------------------------------

% Loads training and testing data
load(input_data_file,'trn','tst');
[Dim,num_data] = size( trn.X );

% Data partitioning for cross-validation
[itrn,itst]=crossval(num_data,num_folds);

% Tuning linear PCA model
%-------------------------------------------------------
Mse = [];

for new_dim = New_Dim_Range,
    
 fprintf('\nnew_dim = %d\n', new_dim);
   
 cv_mse = 0;  
 for i=1:num_folds,
    
   fprintf('\n');
    
   trn_X = trn.gnd_X(:,itrn{i});
   val_gnd_X = trn.gnd_X(:,itst{i});
   val_corr_X = trn.X(:,itst{i});

   
   fprintf('Computing Linear PCA...');
   lin_model = pca(trn_X, new_dim);
   fprintf('done\n');

   fprintf('Projecting validation data...');
   val_reconst_X = pcarec( val_corr_X, lin_model );
   fprintf('done.\n');

   dummy = (val_reconst_X - val_gnd_X).^2;
  
   mse = sum(dummy(:))/size(val_gnd_X,2);
  
   fprintf('folder %d/%d: validation errors mse=%f\n', ...
      i, num_folds, mse);
   
   cv_mse = cv_mse + mse;
 end
 
 cv_mse = cv_mse/num_folds;
 
 Mse(find(new_dim==New_Dim_Range)) = cv_mse;
 
 fprintf('new_dim = %d: mse = %f\n', new_dim, cv_mse);
end

% take the best dimension
%--------------------------------------------------
[dummy,inx] = min(Mse);
fprintf('\nMin(mse) = %f, dim = %f\n', ...
   Mse(inx), New_Dim_Range(inx) );

% compute PCA with tbe best dimesion and all training data
%----------------------------------------------------------
fprintf('Computing optimal Kernel PCA...');
lpca_model = pca( trn.X, New_Dim_Range(inx) );
fprintf('done.\n');

if isempty(output_data_file),
  figure; hold on;
  xlabel('dim'); ylabel('mse');

  plot(New_Dim_Range,Mse);
else
   save(output_data_file,'New_Dim_Range',...
      'Mse','num_folds','input_data_file',...
      'output_data_file','lpca_model');
end

% plot denosing in 2D case only
%-------------------------------------
if Dim == 2 & isempty(output_data_file),
  X = pcarec(tst.X, lin_model );

  mse = sum(sum((X-tst.gnd_X).^2 ));
  fprintf('\ntest mse=%f\n', mse);
  
  figure; hold on;
  h0=ppatterns(tst.gnd_X,'r+');
  h1=ppatterns(tst.X,'gx');
  h2=ppatterns(X,'bo');
  legend([h0 h1 h2],'Ground truth','Noisy','Reconst');
end

% EOF
