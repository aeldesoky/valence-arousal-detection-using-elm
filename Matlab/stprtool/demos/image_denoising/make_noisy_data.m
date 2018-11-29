% MAKE_NOISY_DATA Adds Gaussian noise to USPS database.
%
% Description:
%  It adds Gaussian noise to the USPS images. The input
%  file usps.mat contains training trn.X and testing 
%  tst.X part. This script generates file usps_noisy
%  which contains
%    trn.gnd_X [256x7291] Original training USPS data.
%    trn.X     [256x7291] USPS data with added Gaussian noise.
%    trn.y     [1x7291] Labels (1..10).
%
%    tst.gnd_X [256x2007] Original testing USPS data.
%    tst.X     [256x2007] USPS data with added Gaussian noise.
%    tst.y     [1x2007] Labels (1..10).
%    

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 07-jun-2004, VF

% setting
%---------------------------------------------
% signal to noise ratio
snr = 1;

input_data_file = '/home.dokt/xfrancv/data/usps/usps.mat';
output_data_file = '/home.dokt/xfrancv/data/usps/usps_noisy.mat';

% load input file
orig = load(input_data_file);

% add noise
trn.X = orig.trn.X + randn(size(orig.trn.X))*(std(orig.trn.X(:))/snr);
trn.y = orig.trn.y;
trn.gnd_X = orig.trn.X;
tst.X = orig.tst.X + randn(size(orig.tst.X))*(std(orig.tst.X(:))/snr);
tst.y = orig.tst.y;
tst.gnd_X = orig.tst.X;

% save it to file
save(output_data_file,'tst','trn');

% EOF