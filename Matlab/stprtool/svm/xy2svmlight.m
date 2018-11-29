function xy2svmlight(data,file_name)
% XY2SVMLIGHT Converts data set to SVM^{light} format.
%
% Synopsis:
%  xy2svmlight(data,file_name).
%
% Description:
%  This function saves training data to text file required
%  by the SVM^{light} software. 
%
% Input:
%  data.X [dim x num_data] training data stored as column vectors.
%  data.y [1 x num_data] labels of traning data; possible values 
%   are 1 (first class) and 2 (second class).
%
% Output:
%  Text file 'file_name' in SVM^{Light} format.
%
% See also SVMLIGHT.
%

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 14-Jan-2003, VF
% 30-apr-2001, V. Franc, created

fid = fopen( file_name, 'w+');

dim=size(data.X,1);
num_data=size(data.X,2);

txt = zeros(1,2*dim);
inx1 = 1:2:2*dim;
txt(inx1) = 1:dim;
inx2 = 2:2:2*dim;

for i=1:num_data,

  if data.y(i) == 1, 
    fprintf(fid,'+1 '); 
  elseif data.y(i) == 2,
    fprintf(fid,'-1 '); 
  else
    fprintf(fid,'0 '); 
  end      
  
  txt(inx2) = data.X(:,i);
  
  fprintf( fid, '%d:%f ', txt );
  fprintf(fid,'\n');
end

fclose(fid);

return;
