% SHOW_DENOIS_TUNING Plots curves of tuning stage of Kernel PCA 
%  and liner PCA trained for denosing of UPSP handwritten numerals.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 06-jun-2004, VF
% 28-apr-2004, VF

% Models with results of tuning
%-------------------------------------------------------
kpca_filename = 'USPSModelGreedyKPCA.mat';  % kpca model
lpca_filename = 'USPSModelLinPCA.mat';      % linear PCA model

% Kernel KPCA 
%-------------------------------------------------------
load( kpca_filename );

figure; hold on; title('Tuning greedy KPCA');
xlabel('\sigma'); ylabel('\epsilon_{MS}');

h = [];
clear Str;
for i=1:length(New_Dim_Range),
  h = [h, plot(Arg_Range, Mse(i,:),marker_color(i) )];
  Str{i} = sprintf('dim = %d', New_Dim_Range(i));
end
legend(h,Str);

% Linear PCA 
%-------------------------------------------------------
load( lpca_filename ); 

figure; hold on; title('Tuning linear PCA');
xlabel('dim'); ylabel('\epsilon_{MS}');

plot(New_Dim_Range,Mse);

% EOF