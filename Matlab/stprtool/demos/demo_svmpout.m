% DEMO_SVMPOUT Fitting a posteriory probability to SVM output.
%
% A posteriory probability p(y==1|f(x)) of the first class
% given SVM output f(x) is assumed to be sigmoid function. 
% Parameters A(1) and A(2) of the sigmoid function 
%   p(y==1|f(x)) = 1/(1+exp(A(1)*f(x)+A(2)) 
% are estimated using Maximum-Likelihood [Platt99a]. 
%
% The Gaussian mixture model (GMM) is fitted to the SVM output
% and the a posteriory probability are computed for 
% comparison to the ML estimate.
% 
% The ML estimation of the sigmoid function is imlemented 
% in 'mlsigmoid' (see 'help mlsigmoid' for more info).
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 03-jun-2004, VF
% 6-May-2003, VF

help demo_svmpout;
echo on;

% load training data
data = load('riply_trn');

% train SVM model
svm_model = smo(data,struct('ker','rbf','arg',1,'C',10));

% plot SVM decision boundary 
figure; ppatterns(data); psvm(svm_model);

% compute SVM output
[dummy,svm_output.X] = svmclass(data.X,svm_model);
svm_output.y = data.y;

% ML fitting of sigmod to SVM ouput
%-------------------------------------------------------
sigmoid_model = mlsigmoid(svm_output);

% plot fitted probability
fx = linspace(min(svm_output.X), max(svm_output.X), 200);
sigmoid_apost = sigmoid(fx,sigmoid_model);

figure; hold on; 
xlabel('svm output f(x)'); ylabel('p(y=1|f(x))');
hsigmoid = plot(fx,sigmoid_apost,'k');
ppatterns(svm_output);

% ML estimation of GMM model of SVM output
%-------------------------------------------------------
gmm_model = mlcgmm( svm_output );

% compute a posteriory probability
pcond = pdfgauss( fx, gmm_model);
gmm_apost = (pcond(1,:)*gmm_model.Prior(1))./...
    (pcond(1,:)*gmm_model.Prior(1)+(pcond(2,:)*gmm_model.Prior(2)));

hgmm = plot(fx,gmm_apost,'g');
hcomp = pgauss(gmm_model); 

legend([hsigmoid,hgmm,hcomp],'P(y=1|f(x)) ML-Sigmoid',...
    'P(y=1|f(x)) ML-GMM','P(f(x)|y=1) ML-GMM','P(f(x)|y=2) ML-GMM');

echo off;

% EOF