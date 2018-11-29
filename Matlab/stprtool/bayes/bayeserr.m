function [risk,eps1,eps2,inter1]=bayeserr(model)
% BAYESERR Bayesian risk for 1D Gaussians and 0/1-loss.
%
% Synopsis:
%  [risk,eps1,eps2,inter1] = bayeserr(model)
%
% Description:
%  This function computes Bayesian risk of a classifier 
%  with the following assumptions:
%   - 1/0 loss function (risk = expectation of misclassification).
%   - Binary classification.
%   - Class conditional probabilities are univariate Gaussians.
%
% Input:
%  model [struct] Mixture of two univariate Gaussians.
%   .Mean [1x2] Mean values [Mean1 Mean2].
%   .Cov [1x2] Covariances [Cov1 Cov2].
%   .Prior [1x2] A priory probabilities.
% 
% Output:
%  risk [1x1] Bayesian risk for an optimal classifier.
%  eps1 [1x1] Integral of p(x|k=1) over x in L2, where
%    L2 is the area where x is classified to the 2nd class.
%  eps2 [1x1] Integral of p(x|k=2) over x in L1, where
%    L1 is the area where x is classified to the 1st class.
%  inter1 [1x2] or [1x4] One or two intervals describing L1.
%
% Example:
%  model = struct('Mean',[0 0],'Cov',[1 0.4],'Prior',[0.4 0.6]);
%  figure; hold on; 
%  h = pgmm(model,struct('comp_color',['r' 'g'])); 
%  legend(h,'P(x)','P(x|y=1)*P(y=1)','P(x|y=2)*P(y=2)');
%  [risk,eps1,eps2,interval] = bayeserr(model)
%  a = axis;
%  plot([interval(2) interval(2)],[a(3) a(4)],'k');
%  plot([interval(3) interval(3)],[a(3) a(4)],'k');
%
% See also 
%  BAYESDF, BAYESCLS
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 22-oct-2009, VF, fixed bug on line 157; bug reported 2009-10-08 by krejci.filip@gmail.com
% 20-mar-2006, VF, A mistake in help removed; bug reported by O. Sychrovksy.
% 02-may-2004, VF
% 19-sep-2003, VF
% 27-Oct-2001, VF

% allow input to be a cell
model = c2s(model);

% univariate variances can be given as a vector
if size(model.Cov,1) ~= size(model.Cov,2), 
  model.Cov = reshape(model.Cov,1,1,2); 
end

% finds out decision function which is generaly quadratic
quad_model=bayesdf(model);
a = quad_model.A;
b = quad_model.B;
c = quad_model.C;

% get parameters
p1=model.Prior(1); p2=model.Prior(2);
c1=model.Cov(:,:,1); c2=model.Cov(:,:,2);
m1=model.Mean(1); m2=model.Mean(2);

% Split X into X1 and X2 according to the computed quadratic 
% discriminat function ax^2 + bx + c = 0 and computes 
% eps1 and eps2.

if a==0,
   % The decision function is linear, i.e. in 1D it is a
   % single threshold.
   th=-c/b;
   inter1=[th,inf];
   
   % gets label for the interval (th,inf)
   class=classify(th+1,p1,p2,m1,m2,c1,c2);
   
   if abs(c)==inf,
     risk=0;
     if class==1,,
       eps1=0;
       eps2=1;
       inter1=[-inf,inf];
       return;
     else
      eps1=1;
      eps2=0;
      inter1=[];
     end
   end   
      
   eps1=1-erfc2(th,m1,sqrt(c1));
   eps2=erfc2(th,m2,sqrt(c2));  
         
   if class==2,
      % swap eps1 and eps2
      tmp=eps2;
      eps2=eps1;
      eps1=tmp;
   end
   
else
   % The decision function is quadratic, i.e. in 2d
   % there exis two thresholds which determine three intervals.
   
   D=b^2-4*a*c;
   if D > 0,
      th1=(-b-sqrt(D))/(2*a);
      th2=(-b+sqrt(D))/(2*a);
      
      if th1 > th2,
         tmp=th1;
         th1=th2;
         th2=tmp;         
      end;
      
      % finds out label for the interval [th1,th2].
      class=classify((th2+th1)/2,p1,p2,m1,m2,c1,c2);
      
      if class==2
         % integral eps2 = int_th2^inf + int_{-inf}^th1
         eps2 = 1 + erfc2(th2,m2,sqrt(c2)) - erfc2(th1,m2,sqrt(c2));          
         % integral eps1= int_th1^th2
         eps1 = erfc2(th1,m1,sqrt(c1)) - erfc2(th2,m1,sqrt(c1));         
         
         inter1=[-inf,th1,th2,inf];
      else
         % integral eps1 = int_th2^inf + int_{-inf}^th1
         eps1 = 1 + erfc2(th2,m1,sqrt(c1)) - erfc2(th1,m1,sqrt(c1));  
         % integral eps2= int_th1^th2
         eps2=erfc2(th1,m2,sqrt(c2))-erfc2(th2,m2,sqrt(c2));
         
         inter1=[th1,th2];
      end
      
   else
      % finds out label for the interval [-inf,inf].
      class=classify(0,p1,p2,m1,m2,c1,c2);
      
      if class == 1,
         eps1=0;
         eps2=1;
         inter1=[-inf,inf];
      else
         eps1=1;
         eps2=0;
         inter1=[];
      end
%      risk=0; % bug reported 2009-10-08 by krejci.filip@gmail.com
   end
end

% computes the Bayesian risk 
risk = p1*( eps1 - eps2 ) + eps2;

return;

%-----------------------------------------------
function class = classify(x,p1,p2,m1,m2,c1,c2)
% finds out to which class the given x belongs

if p1==1,     % only the 1st class is possible
   class=1;
elseif p2==1, % only the 2nd class is possible
   class=2;
elseif pdfn(x,m1,c1)*p1 > pdfn(x,m2,c2)*p2,
   class =1;
else
   class=2;
end

function p=pdfn(x,m,c)
 p=exp(-1/2*mahalan(x,m,c))/((2*pi)^(1/2) * sqrt(det(c)));
