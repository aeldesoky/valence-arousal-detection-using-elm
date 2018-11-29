function save_chars(data);
% SAVE_CHARS Saves images to file.
%
% Synopsis:
%  save_chars(data)
%
% Description:
%  This is an auxiliary function which saves items of the 
%  structure data into the file its name is taken from the 
%  global variable FileName. The function is used to when
%  examples of numerlas are being collected.
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 09-sep-03, VF

global FileName;

X=data.X;
img_size = data.img_size;
save(FileName,'X','img_size');

figure;
showim(X);

fprintf('Characters saved to file: %s\n',FileName);

return;
% EOF