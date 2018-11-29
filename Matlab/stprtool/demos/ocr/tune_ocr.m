% TUNE_OCR Tunes SVM classifier for OCR problem.
%
% Description:
%  The following steps are performed:
%   - Training set is created from data in directory ExamplesDir.
%   - Multi-class SVM is trained for a set of arguments and 
%     regularization constants. The best model is selected 
%     based on the cross-validation error.
%    
 
% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 04-jun-2004, VF
% 09-sep-2003, VF

cd /home.dokt/xfrancv/work/new_stprtool/;
stprpath;
cd /home.dokt/xfrancv/work/new_stprtool/demos/ocr;

% Setting 
%===================================
ExamplesDir = '../../data/ocr_numerals/';  % input directory with exmaples
OCRTuningFileName = 'ocrtuning';  % output file with result of tuning

% Model setting for multi-class SVM
options.ker = 'rbf';        % kernel type
options.arg = [1 2.5 5 7.5 10] ;  % kernel argument
options.C = [inf];          % regularization constant
options.verb = 1;           % display progress info
options.solver ='oaosvm';
options.num_fold = 5;
options.svm_options.solver = 'svmlight';

% Create training set 
%====================================
fprintf('Creating training set:\n');
TrainingDataFile = [ExamplesDir 'OcrTrndata.mat'];
mergesets( ExamplesDir, TrainingDataFile );
data = load(TrainingDataFile );

% Tuning SVM model
%====================================

fprintf('Tuning multi-class SVM classifier.\n');
[model,Error] = evalsvm(data,options);

fprintf('\nSaving results to: %s\n',OCRTuningFileName);
save(OCRTuningFileName,'model','Error','options');

% Visualization
%====================================
%figure; 
%load(OCRTuningFileName)
%figure; mesh(options.arg,options.C,Errors);
%hold on; xlabel('arg'); ylabel('C');
%figure; contour(options.arg,options.C,Errors);
%hold on; xlabel('arg'); ylabel('C');

% EOF