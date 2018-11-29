function t=gat1dsearch(MI,SG,alpha,dalpha,tmax,tdelta)
% GAT1DSEARCH 1D search along improving direction in the GAT.
%
% Synopsis:
%  t=gat1dsearch(MI,SG,alpha,dalpha,tmax,tdelta)
%
% Description:
%  Auxciliary function for the 'ganders' algorithm.
%  It implements 1D-search based on the cutting interval 
%  algorithm according to the Fibonacci series. 
%
% See also 
%  GANDERS
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 21-may-2004, VF
% 17-sep-2003, VF
% 24. 6.00 V. Hlavac, comments polished.

LO_TH=0;

% default setting
if nargin < 5,
   tmax = inf;
   delta=1e-6;
elseif nargin < 6,
   delta=0;
end

% get dimension N and the # of distributions
K = size(MI,2);
N = size(MI,1);

% compute constants
for j = 1:K,
   s(j)= alpha'*MI(:,j);
   ss(j) = dalpha'*MI(:,j);
   ds(j) = ss(j) - s(j);
   sga(j) = alpha'*SG(:,:,j)*alpha;
   sgd(j) = dalpha'*SG(:,:,j)*dalpha;
   sgad(j) = dalpha'*SG(:,:,j)*alpha;
end


% first step
F1=1;
F2=1;
tbeg=0;
tend=1;
tmid=0.5*(tend+tbeg);

fmid=max([LO_TH,min( (s+tmid*ds)./sqrt( (1-tmid)^2*sga + 2*tmid*(1-tmid)*sgad + tmid^2*sgd ) )]);
fbeg=max([LO_TH,min( (s+tbeg*ds)./sqrt( (1-tbeg)^2*sga + 2*tbeg*(1-tbeg)*sgad + tbeg^2*sgd ) )]);

if sqrt( (1-tend)^2*sga + 2*tend*(1-tend)*sgad + tend^2*sgd ) == 0,
 fend=0;
else
 fend=max([LO_TH,min( (s+tend*ds)./sqrt( (1-tend)^2*sga + 2*tend*(1-tend)*sgad + tend^2*sgd ) )]);
end


% start up
stop=0;
while stop==0 & tmax > 0,
   tmax=tmax-1;

   % store fmid
   oldfmid=fmid;

   % Fibonacci, F(k)=F(k-1)+F(k-2)
   F=F2+F1;

   % find larger interval
   if (tmid-tbeg) < (tend-tmid),
      % new bound
      t=tmid+F1*(tend-tmid)/F;

      fvalue=max([LO_TH,min( (s+t*ds)./sqrt( (1-t)^2*sga + 2*t*(1-t)*sgad + t^2*sgd ) )]);

      if fvalue < fmid,
         tend=t;
         fend=fvalue;
      else
         tbeg=tmid;
         fbeg=fmid;
         tmid=t;
         fmid=fvalue;
      end
   else
      % new bound
      t=tbeg+F1*(tmid-tbeg)/F;

      fvalue=max([LO_TH,min( (s+t*ds)./sqrt( (1-t)^2*sga + 2*t*(1-t)*sgad + t^2*sgd ) )]);

      if fvalue < fmid,
         tbeg=t;
         fbeg=fvalue;
      else
         tend=tmid;
         fend=fmid;
         tmid=t;
         fmid=fvalue;
      end
   end

   % update Fibonacci F(k-2)=F(k-1) and F(k-1)=F(k);
   F2=F1;
   F1=F;

   % stop condition
   if tend-tbeg < tdelta,
      stop=1;
   end

end

% get the bigest value
fvalues=[fbeg fmid fend];
tvalues=[tbeg tmid tend];

[fmax, imax]=max(fvalues);
tmaxim=tvalues(imax);

% compute new alpha
%alpha=alpha*(1-tmaxim)+dalpha*tmaxim;
t=tmaxim;

return;

% debugging
if 1==1,
   vals=[];
   
   for t=0:0.01:1,
   
      fvalue=min( (s+t*ds)./sqrt( (1-t)^2*sga + 2*t*(1-t)*sgad + t^2*sgd ) );

      vals=[vals,fvalue];
   end

   figure;
   hold on;
   plot(0:0.01:1,vals,'g');
   win=axis;
   line([tmid tmid],[ win(3) win(4)],'Color','k');
   line([0 1],[vals(1) vals(1)],'Color','r');
   drawnow;
end

pause;
return;

