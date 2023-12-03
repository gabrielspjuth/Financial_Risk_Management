%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: Assignment_2.m
% Author: Vilmer Gimbringer & Gabriel Spjuth
% Date: December 3, 2023
% Description: Assignment 2 in Financial Risk Management.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1: VaR and Expected Shortfall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load the data from the file 'data.txt'.
[S, headers, raw] = xlsread("timeSeries.xlsx");

portfolio_value = 10;

n = 15; % There are 15 different stocks in the portfolio
omega = 1 / n;

t = 1/52; % Time of observed data
confidence_level_95 = 0.95;
confidence_level_975 = 0.975;
confidence_level_99 = 0.99;

R = zeros(1646, 15);

prices = 1647;
stocks = 15;

%% Subtask 1 (a): Calculation of Value at Risk

for j = 1:stocks
    for i = prices:-1:2
        R(i-1, j) = ((S((i-1), j+1)) - S(i, j+1)) / S(i, j+1);
    end
end



























%% Subtask 1.2: Clean Data
% Remove any NaN values from the data.
cleaned_data = cleanData(data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 2: Model Training
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subtask 2.1: Split Data
% Split the data into training and testing sets.
[train_data, test_data] = splitData(cleaned_data);

%% Subtask 2.2: Train Model
% Train a machine learning model using the training data.
model = trainModel(train_data);

%% Subtask 2.3: Evaluate Model
% Evaluate the model using the testing data.
accuracy = evaluateModel(model, test_data);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 3: Results Analysis
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subtask 3.1: Display Results
% Display the accuracy of the trained model.
disp(['Model Accuracy: ' num2str(accuracy)]);

%% Subtask 3.2: Visualize Results
% Visualize the results using plots or other visualization techniques.
visualizeResults(model, test_data);
