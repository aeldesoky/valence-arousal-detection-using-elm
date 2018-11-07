%% Clear all/ close figs
close all 
clear
clc

%% Load Data
addpath(genpath('C:\Users\mjad9\Desktop\DEAP\TEAP'));
addpath(genpath('C:\Users\mjad9\Desktop\DEAP\MATLAB_data_preprocessed'));

% Just testing with first subject
s1 = load('s01.mat');
data = s1.data;

%% Exploring data -- This is where I keep getting errors
% trying multiple methods
