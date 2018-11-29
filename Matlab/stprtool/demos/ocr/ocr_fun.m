function ocr_fun(data)
% OCR_FUN Calls OCR classifier and displays result.
%
% Synopsis:
%  ocr_fun(data)
% 
% Description:
%  This function classifies images of characters stored as columns 
%  of the matrix data.X. The output is displayed in a grid 5 x 10.
%
% Input:
%  data.X [dim x (5*10)] Input images store as column vectors.
%    The images are assumed to be taken from grid 5x10.
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 04-jun-2004, VF
% 09-sep-2003, VF

% loads the multi-class SVM model 
OCRFileName = 'ocrmodel.mat';
ocr = load( OCRFileName );

% make sure that only the filled subwindows will be classified
filled_inx = find( sum(data.X) ~= 0);
if isempty(filled_inx), return; end
y = -ones(1,length(data.X)); % non-filled subwindows are labeled by -1

% classify filled subwindows
%-------------------------------------------------------
[y(filled_inx),dfce] = feval( ocr.fun, data.X(:,filled_inx), ocr);

% displays results
%-------------------------------------------------------
figure(10); clf;

% plot grid 5 x 10 
axis equal;
axis([0 1 0 0.5]);
hold on;
plot( [ 0 0 1 1 0 ], [ 0 .5 .5 0 0 ] );
for i = 1:9, plot( [i/10 i/10],[0 .5] ); end
for i = 1:4, plot( [0 1],[i/10 i/10] );  end

% display recognized numerals
for i=1:5,
   for j = 1:10,
      inx = j+(i-1)*10;
      if y(inx) ~= -1,
        character = ocr.labels(y(inx)); 
      
        h=text(0.05+(j-1)/10,-0.05+(5-i+1)/10,character);
        set(h,'fontsize',25,'VerticalAlignment','middle',...
           'HorizontalAlignment','center');
      end
   end
end

return;
% EOF