%% Clear all/ close figs
close all 
clear
clc

%% Setting Parameters
PARTICIPANTS_NUM = 32;
VIDEOS_NUM       = 40;
%% Load Data
addpath(genpath('Matlab/TEAP-master'));
user = 'matt';

if(strcmp(user, 'amir'))
    data_path = '~/Desktop/DEAP/MATLAB_data_preprocessed/';
    feedback_path = '~/Desktop/DEAP/participant_ratings.csv';
else
    data_path = 'C:\Users\mjad9\Desktop\DEAP\MATLAB_data_preprocessed\';
    feedback_path = 'C:\Users\mjad9\Desktop\DEAP\participant_ratings.csv';
end

if ~exist([data_path '/s30_eeglab.mat'],'file')
    loading_DEAP(data_path);
end 

feedbacks = readtable(feedback_path);
%% Extracting Features
for participant = 1:PARTICIPANTS_NUM
    eeglab_file = sprintf('%ss%0.2d_eeglab.mat', data_path, participant);
    disp(eeglab_file)
    bulk = Bulk_load(eeglab_file);
    
    means_array = [];
    std_array   = [];
    kurtosis_array = [];
    skewness_array = [];
    entropy_array =[];
    energy_array =[];
    katz_array = [];
    
    for epoch = 1:VIDEOS_NUM
        %extracting EMG features
        [features(participant,epoch).EMG_feats, features(participant,epoch).EMG_feats_names] = EMG_feat_extr(bulk(epoch));
        %extracting EEG features
        [features(participant,epoch).EEG_feats, features(participant,epoch).EEG_feats_names] = EEG_feat_extr(bulk(epoch));
        %extracting GSR features
        [features(participant,epoch).GSR_feats, features(participant,epoch).GSR_feats_names] = GSR_feat_extr(bulk(epoch));
        %extracting BVP features
        %[features(participant,epoch).BVP_feats, features(participant,epoch).BVP_feats_names] = BVP_feat_extr(bulk(epoch));
        %extracting respiration features
        [features(participant,epoch).RES_feats, features(participant,epoch).RES_feats_names] = RES_feat_extr(bulk(epoch));
        
        feedback = feedbacks(feedbacks.Participant_id==participant & feedbacks.Experiment_id==epoch,:);
        
        features(participant,epoch).feedback.felt_arousal     = feedback.Arousal;
        features(participant,epoch).feedback.felt_valence     = feedback.Valence;
        features(participant,epoch).feedback.felt_dominance   = feedback.Dominance;
        features(participant,epoch).feedback.felt_liking      = feedback.Liking;
        features(participant,epoch).feedback.felt_familiarity = feedback.Familiarity;
        
        video_signal   = Bulk_get_signal(bulk(epoch), 'EEG');
        video_raw_data = cell2mat(struct2cell(video_signal.raw));
        [mean, std, kurtosis, skewness] = Signal_feat_stat_moments(video_raw_data, 'Skip');
        energy = Signal_feat_energy(video_raw_data, 'Skip');
        for i=1:32
            katz(i) = Katz_FD(video_raw_data(i,:));
        end
        means_array = [means_array, mean];
        std_array = [std_array, std];
        kurtosis_array = [kurtosis_array, kurtosis];
        skewness_array = [skewness_array, skewness];
        energy_array = [energy_array, energy'];
        entropy_array = [entropy_array, entropy(video_raw_data)'];
        katz_array = [katz_array, katz];
        
        fprintf('extracted all the features for subject %d epoch %d\n',participant, epoch);
    end
    
    means_matrix{participant} = means_array;
    std_matrix{participant} = std_array;
    kurtosis_matrix{participant} = kurtosis_array;
    skewness_matrix{participant} = skewness_array;
    entropy_matrix{participant} =  entropy_array;
    energy_matrix{participant} = energy_array;
    katz_matrix{participant} = katz_array;
end

%% Exploring the data
delta      = [];
theta      = [];
slow_alpha = [];
alpha      = [];
beta       = [];
gamma      = [];
valence_labels   = [];
arousal_labels   = [];
dominance_labels = [];
liking_labels    = [];
video       = 1;

for participant = 1:PARTICIPANTS_NUM
    for video = 1:VIDEOS_NUM
        delta      = [delta, features(participant,video).EEG_feats(1, :)];
        theta      = [theta, features(participant,video).EEG_feats(2, :)];
        slow_alpha = [slow_alpha, features(participant,video).EEG_feats(3, :)];
        alpha      = [alpha, features(participant,video).EEG_feats(4, :)];
        beta       = [beta, features(participant,video).EEG_feats(5, :)];
        gamma      = [gamma, features(participant,video).EEG_feats(6, :)];
        
        valence   = features(participant,video).feedback.felt_valence;
        arousal   = features(participant,video).feedback.felt_arousal;
        dominance = features(participant,video).feedback.felt_dominance;
        liking    = features(participant,video).feedback.felt_liking;

        valence_labels   = [valence_labels, repmat(valence > 5, 1, 32)];
        arousal_labels   = [arousal_labels, repmat(arousal > 5, 1, 32)];
        dominance_labels = [dominance_labels, repmat(dominance > 5, 1, 32)];
        liking_labels    = [liking_labels, repmat(liking > 5, 1, 32)];       
    end
end

means_features = cell2mat(means_matrix');
std_features = cell2mat(std_matrix');
kurtosis_features = cell2mat(kurtosis_matrix');
skewness_features = cell2mat(skewness_matrix');
entropy_features = cell2mat(entropy_matrix');
energy_features = cell2mat(energy_matrix');
katz_features = cell2mat(katz_matrix');
features_array = [delta' theta' slow_alpha' alpha' beta' gamma' means_features(:) std_features(:) kurtosis_features(:) skewness_features(:) entropy_features(:) energy_features(:) katz_features(:)];
features_array = zscore(features_array);
labels         = [valence_labels' arousal_labels' dominance_labels' liking_labels'];

for i = 1:size(valence_labels,2)
    if(valence_labels(i) == 1 && arousal_labels(i) == 1)
        labels_2(i) = 1;
    elseif (valence_labels(i) == 1 && arousal_labels(i) == 0)
        labels_2(i) = 2;
    elseif (valence_labels(i) == 0 && arousal_labels(i) == 1)
        labels_2(i) = 3;
    else
        labels_2(i) = 4;
    end        
end
%% Visualizing Features
figure();
varnames = {'? PSD' '? PSD' 'Slow ? PSD' '? PSD' '? PSD' '? PSD' 'Mean' 'STD ' 'Kurtosis' 'Skewness' 'Entropy' 'Energy' 'FD'};
gplotmatrix(features_array,[],labels(:,1)',['r' 'b'],[],[],'on','grpbars',varnames, varnames);
title('Matrix of Feature Scatterplots - Valence')
figure();
gplotmatrix(features_array,[],labels(:,2)',['r' 'b'],[],[],'on','grpbars',varnames, varnames);
title('Matrix of Feature Scatterplots - Arousal')

%% Filterling channels
Fp1=1;AF3=2;F3=3;F7=4;FC5=5;FC1=6;C3=7;T7=8;CP5=9;CP1=10;P3=11;P7=12;PO3=13;
O1=14;OZ=15;Pz=16;Fp2=17;AF4=18;Fz=19;F4=20;F8=21;FC6=22;FC2=23;Cz=24;C4=25;
T8=26;CP6=27;CP2=28;P4=29;P8=30;PO4=31;O2=32;
chosen_channels = zeros(32,1);
chosen_channels([Fp1 Fp2 F3 F4 F7 F8 FC5 FC6 FC1 FC2 AF3 AF4 C3 C4 T7 T8 Fz Cz]) = 1;
chosen_channels = repmat(chosen_channels, 32*40, 1);
filtered_features = features_array(chosen_channels==1, :);
filtered_arousal_labels = arousal_labels(1, chosen_channels==1)';
filtered_valence_labels = valence_labels(1, chosen_channels==1)';

%% Feature Selection
partition = cvpartition(filtered_labels, 'KFold', 10);
options   = statset('display', 'iter');

criterion = @(XT,yT,Xt,yt) ...
      (sum(yt ~= classify(Xt,XT,yT, 'quadratic')));

[fs,history] = sequentialfs(criterion,filtered_features,filtered_labels,'direction', 'backward', 'cv',partition,'options',options);

%% Selected Features
selected_features = filtered_features(:, [1 2 3 4 6 11]);
selected_arousal_labels = filtered_arousal_labels;
selected_valence_labels = filtered_valence_labels;

%% ELM
for fold=1:10
    training_indices = training(partition,fold);
    testing_indices  = test(partition,fold);
    elm_arousal_features = [selected_valence_labels selected_features];
    elm_valence_features = [selected_arousal_labels selected_features];

    testing_arousal_set = elm_arousal_features(testing_indices,:);
    training_arousal_set = elm_arousal_features(training_indices,:);
    testing_valence_set = elm_valence_features(testing_indices,:);
    training_valence_set = elm_valence_features(training_indices,:);

    dlmwrite('training_arousal_set.txt', training_arousal_set);
    dlmwrite('testing_arousal_set.txt', testing_arousal_set);
    dlmwrite('training_valence_set.txt', training_valence_set);
    dlmwrite('testing_valence_set.txt', testing_valence_set);

    %[ELMArousalTrainingTime(fold), ELMArousalTestingTime(fold), ELMArousalTrainingAccuracy(fold), ELMArousalTestingAccuracy(fold)] = ELM('training_arousal_set.txt', 'testing_arousal_set.txt',  1, 200, 'sig');
    [ArousalTrainingTime2(fold), ArousalTestingTime2(fold), ArousalTrainingAccuracy2(fold), ArousalTestingAccuracy2(fold)] = elm_kernel('training_arousal_set.txt', 'testing_arousal_set.txt', 1, 1, 'RBF_kernel',100);
    %[ELMValenceTrainingTime(fold), ELMValenceTestingTime(fold), ELMValenceTrainingAccuracy(fold), ELMValenceTestingAccuracy(fold)] = ELM('training_valence_set.txt', 'testing_valence_set.txt',  1, 200, 'sig');
    [ValenceTrainingTime2(fold), ValenceTestingTime2(fold), ValenceTrainingAccuracy2(fold), ValenceTestingAccuracy2(fold)] = elm_kernel('training_valence_set.txt', 'testing_valence_set.txt', 1, 1, 'RBF_kernel',100);
end
%% SVM
for fold=1:10
    training_indices = training(partition,fold);
    testing_indices  = test(partition,fold);
    svm_features = filtered_features(training_indices,:);
    v_labels   = filtered_valence_labels;
    svm_v_labels = v_labels(training_indices, :);
    a_labels   = filtered_arousal_labels(:,1);
    svm_a_labels = a_labels(training_indices, :);
    %Valence
    tic
    svm_valence_classifier = fitcsvm(svm_features, svm_v_labels, 'KernelFunction', 'rbf');
    svm_valence_time(fold) = toc;
    valence_prediction = predict(svm_valence_classifier, features_array(testing_indices,:));
    svm_valence_accuracy(fold) = sum(valence_prediction == v_labels(testing_indices,1))/size(valence_prediction,1);
    %Arousal
    tic
    svm_arousal_classifier = fitcsvm(svm_features, svm_a_labels, 'KernelFunction', 'rbf');
    svm_arousal_time(fold) = toc;
    arousal_prediction = predict(svm_arousal_classifier, features_array(testing_indices,:));
    svm_arousal_accuracy(fold) = sum(arousal_prediction == a_labels(testing_indices,1))/size(arousal_prediction,1);
end
%% KNN
for fold=1:10
    training_indices = training(partition,fold);
    testing_indices  = test(partition,fold);
    knn_features = filtered_features(training_indices,:);
    v_labels   = filtered_valence_labels(:,1);
    a_labels   = filtered_arousal_labels(:,1);
    knn_valence_labels = v_labels(training_indices,:);
    knn_arousal_labels = a_labels(training_indices,:);
    % Valence
    tic
    knn_valence_classifier = fitcknn(knn_features, knn_valence_labels, 'NumNeighbors',5000); 
    knn_valence_time(fold) = toc;
    valence_prediction = predict(knn_valence_classifier, features_array(testing_indices,:));
    knn_valence_accuracy(fold) = sum(valence_prediction == v_labels(testing_indices,1))/size(valence_prediction,1);

    % Arousal 
    tic
    knn_arousal_classifier = fitcknn(knn_features, knn_arousal_labels, 'NumNeighbors',5000); 
    knn_arousal_time(fold) = toc;
    arousal_prediction = predict(knn_arousal_classifier, features_array(testing_indices,:));
    knn_arousal_accuracy(fold) = sum(arousal_prediction == v_labels(testing_indices,1))/size(arousal_prediction,1);
end
%% KPCA
% options.ker = 'rbf'; 
% options.arg = 4; 
% options.new_dim = 20; 
% X = features_array';
% model = kpca(X,options);
% XR = kpcarec(X,model); 
% figure;
% h1 = ppatterns(X);
% h2 = ppatterns(XR, '+r');
% legend([h1 h2],'Input vectors','Reconstructed');
% 
% [TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy] = ELM('training_set.txt', 'testing_set.txt', 1, 20, 'sig')

%% PCA
pca1 = pca(features_array');