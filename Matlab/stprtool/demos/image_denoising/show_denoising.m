% SHOW_DENOISING Image denosing of USPS hand-written numerals.
%
% Description:
%  The input noisy images are denoised using the 
%  kernel PCA model. The denosing based on linear PCA
%  is also computed for comparision. Imgase of
%   - Ground truth numerlas
%   - Noisy numerals
%   - Numerals denoised by kernel PCA
%   - Numerals denoised by linear PCA
%  are displayed.
%
% See also
%  KPCAREC, PCAREC.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modification:
% 07-jun-2004, VF
% 05-may-2004, VF
% 22-apr-2004, VF

help show_denoising;

% setting 
%-------------------------------------------------------
kpca_filename = 'USPSModelGreedyKPCA.mat'; % kpca model
lpca_filename = 'USPSModelLinPCA.mat';     % linear PCA model

% USPS databes with noisy images (see help make_noisy_usps).
input_data_file = '/home/xfrancv/data/usps/usps_noisy';

% loading
%----------------------------------------
load(kpca_filename,'kpca_model');
load(lpca_filename,'lpca_model');
load(input_data_file,'tst');

% get indices of examples to be denoised
% ---------------------------------------
inx = [];
for i=1:10,
  tmp = find(tst.y == i);
  inx = [inx, tmp(1) ];
end

% get noisy and ground truth numerals
%----------------------------------------
noisy_X = tst.X(:,inx);  
gnd_X = tst.gnd_X(:,inx);

% Kernel PCA and linear PCA denoising
%----------------------------------------
kpca_X = kpcarec( noisy_X, kpca_model);
lpca_X = pcarec( noisy_X, lpca_model);

% display results
%----------------------------------------
h=figure; set(h,'name','Denoised by greedy KPCA');
showim( kpca_X);

h=figure; set(h,'name','Denoised by linear PCA');
showim( lpca_X);

h=figure; set(h,'name','Ground truth');
showim( gnd_X);

h=figure; set(h,'name','Noisy');
showim( noisy_X);

%EOF