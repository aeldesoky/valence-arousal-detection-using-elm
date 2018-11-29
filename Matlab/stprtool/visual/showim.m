function varargout = showim(img,img_size,placement)
% SHOWIM Displays given image(s).
%
% Synopsis:
%  h = showimg(img)
%  h = showimg(img,img_size)
%  h = showimg(img,img_size,placement)
%
% Description:
%  It displays images to the current figure. It assumes that the 
%  images are represented as column vectors stored in the 
%  img [dim x num_img] where dim=height*width. 
%
% Input:
%  img [dim x num_img] Images store as column vectors of size
%    dim = width*height.
%  img_size [2 1] Image size img_size = [height, width]. The image
%    size does have to be geving for rectangular images.
%  placement [2 1] Determines the number of rows and columns
%    in wich the images are displayed, palecement = [rows,columns].
%    if not given then the placement is determined automatically.
%
% Output:
%  h [1 x num_img] Handles of individual axes.
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modificatrions:
% 10-aug-2006, VF
% 15-jun-2004, VF
% 10-sep-2003, VF

[dim,num_data]=size(img);

if nargin < 2,
  sr=sqrt(dim);
  if fix(sr)==sr, 
    img_size = [sr,sr]; 
  else 
    error('Not enough input arguments.'); 
  end
end

%---------------------

if num_data > 1,
  
  if nargin >=3,
    row = placement(1); col = placement(2);
  else

    col=floor(sqrt(4*num_data/3));
    for i=max(1,fix(num_data/10)):fix(sqrt(num_data)),
      if num_data/i == round(num_data/i),
        col= num_data/i;
        break;
      end
    end

    row=ceil(num_data/col);
  end

  cnt=0;
  h = [];
  while cnt < num_data,
    cnt=cnt+1;
    h = [h, subplot(row,col,cnt)];
    
    showim(img(:,cnt),img_size);
  end
  
else
  h = imshow(reshape(img,img_size(1),img_size(2)),[]);
end

if nargout > 0,
  varargout{1} = h;
end

return;

