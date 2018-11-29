function collect_chars( fname )
% COLLECT_CHARS Collect training examples for OCR.
%
% Synopsis:
%  collect_chars( fname )
%  
% Description:
%  It calls function that allows a user to draw examples
%  of characters. The drawn characters are saved to the
%  file fname whenever the middle mouse button is pressed.
%
% Input:
%  fname [string] The filename to which the characters are saved.
%
% Example:
%  Proceed as follows:
%  1) run function
%    collect_chars('my_examples_1');   
%  2) draw examples of numeral '1' (fill up the whole form); 
%     left button -> draw, middle button -> save it to file, 
%     right -> erase.
%  3) repeate the points 1) and 2) for numerals from 2 to 9, 0 
%     (use label 10 for numeral 0, i.e., collect_chars('my_examples_10'))
%

% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 9-sep-03, VF

global FileName;
FileName = fname;

mpaper( {'fun','save_chars','width',16,'height',16});

%EOF