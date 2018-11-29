% TRAIN_OCR Training of OCR classifier based on multiclass SVM.
%
% Description:
%  The following steps are performed:
%   - Training set is created from data in directory ExamplesDir.
%   - Multi-class SVM is trained.
%   - The trained SVM model is saved.
%    
 
% (c) Statistical Pattern Recognition Toolbox, (C) 1999-2003,
% Written by Vojtech Franc and Vaclav Hlavac,
% <a href="http://www.cvut.cz">Czech Technical University Prague</a>,
% <a href="http://www.feld.cvut.cz">Faculty of Electrical engineering</a>,
% <a href="http://cmp.felk.cvut.cz">Center for Machine Perception</a>

% Modifications:
% 4-jun-2004, VF
% 9-sep-2003, VF

% Setting 
%===================================
ExamplesDir = '../../data/ocr_numerals/';  % input directory with examples
OCRFileName = 'ocrmodel';  % output SVM model

% Model setting for multi-class SVM
options.ker = 'rbf';    % kernel type
options.arg = 5;        % kernel argument
options.C = [inf];      % regularization constant
options.verb = 100;       % display progress info

%options.solver ='svmlight';  
options.solver ='smo';   % use if SVM^{light} is not installed

% Create training set 
%====================================
fprintf('Creating training set:\n');
TrainingDataFile = [ExamplesDir 'OcrTrndata.mat'];
mergesets( ExamplesDir, TrainingDataFile );
data = load(TrainingDataFile );

% Training SVM model
%====================================

fprintf('Training multi-class SVM classifier.\n');
%model = oaosvm(data,options);

% Multi-class BSVM with L2-soft margin can be asls used
options.solver = 'imdm';
model = bsvm2(data,options);

% One-Against-All decomposition can be also used 
%options.solver = 'svmlight';
%model = oaasvm(data,options);


% mapping class label y -> character
model.labels = ['1' '2' '3' '4' '5' '6' '7' '8' '9' '0'];

fprintf('Saving found classifier to file %s...', OCRFileName);
savestruct(model,OCRFileName);
fprintf('done.\n');

% EOF