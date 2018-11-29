% Demo on Optical Character Recognition (OCR) for handwritten numerals.
%
% models (dir)  - Different SVM models trained for OCR system.
%
% collect_chars - Collecting training examples for OCR.
% mpaper        - Allows to enter handwritten characters by mouse.
% ocr_fun       - Classifies numerals and displays result.   
% save_chars    - Saves images to file.
% train_ocr     - Trains OCR classifier based on multiclass SVM.
% tune_ocr      - Tunes SVM for OCR problem.
%
% ocrmodel.mat  - SVM model used by OCR system.
%
% Description:
%  First try
%    demo_ocr
%
%  The demo_ocr allows to draw hand-written numerals from 0 to 9 and 
%  classify them by multi-class SVM. The standard computer mouse 
%  is used as the input device.
%
%  The complete set of function necessary to design the OCR are 
%  included. To train your own OCR proceed as follows:
%   - Collecting training set. Use 'collect_chars' to draw your 
%     examples of numerals or other characters. Then save all training 
%     examples to a directory DataDir.
%   - Tuning SVM parameters. Customize the file 'tune_ocr' or
%     use already tuned SVM paramater.
%   - Training SVM model: Customize the file 'train_ocr' and run it.  
%   - Run OCR demo 'demo_ocr'.
%     

% About: Statistical Pattern Recognition Toolbox
% (C) 1999-2003, Written by Vojtech Franc and Vaclav Hlavac
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>
% <a href="http://www.feld.cvut.cz">Faculty of Electrical Engineering</a>
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 05-jun-2004, VF
% 04-jun-2004, VF
% 10-sep-2003, VF
