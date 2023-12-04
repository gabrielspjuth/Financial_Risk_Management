%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% File: Assignment_2.m
% Author: Vilmer Gimbringer & Gabriel Spjuth
% Date: December 3, 2023
% Description: Assignment 2 in Financial Risk Management.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 1: VaR and Expected Shortfall
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

[S, headers, raw] = xlsread("timeSeries.xlsx");

portfolio_value = 10;

n = 15; % There are 15 different stocks in the portfolio
omega = 1 / n;
omega_vector = ones(15, 1) * omega;

t = 1/52; % Time of observed data. Assuming 52 trading weeks per year.
confidence_level_95 = 0.95;
confidence_level_975 = 0.975;
confidence_level_99 = 0.99;

R = zeros(1646, 15);

prices = 1647;
stocks = 15;

%% Subtask 1 (a): Calculation of Value at Risk

% Calculating returns matrix for all stocks from t = 1 ... T + 1
for j = 1:stocks
    for i = prices:-1:2
        R(i-1, j) = ((S((i-1), j+1)) - S(i, j+1)) / S(i, j+1);
    end
end

C = cov(R);

portfolio_volatiliy = sqrt(transpose(omega_vector) * C * omega_vector);

VaR_95 = norminv(confidence_level_95) * portfolio_volatiliy * portfolio_value;

Var_975 = norminv(confidence_level_975) * portfolio_volatiliy * portfolio_value;

VaR_99 = norminv(confidence_level_99) * portfolio_volatiliy * portfolio_value;

%% Subtask 1 (b): Calculation of Relative VaR

portfolio_returns = transpose(omega_vector) * transpose(R); % Index 1 ger portfolio return för T = 1646, index 1646 ger för t = 1

portfolio_returns_logarithmic = log(1+portfolio_returns); % Index 1 ger portfolio return för T = 1646, index 1646 ger för t = 1

lambda = 0.94;

portfolio_volatility_EWMA = zeros(size(portfolio_returns_logarithmic)); % Index 1 ger portfolio return för T = 1646, index 1646 ger för t = 1


% Calculating EWMA volatility for each time step
for i = length(portfolio_returns_logarithmic):-1:1
    if i == 1646 % t = 1
        portfolio_volatility_EWMA(i) = portfolio_returns_logarithmic(i)^2;
    else
        portfolio_volatility_EWMA(i) = lambda * portfolio_volatility_EWMA(i+1) + (1 - lambda) * portfolio_returns_logarithmic(i+1)^2;
    end
end

relative_VaR_95 = zeros(1144, 1);
relative_VaR_99 = zeros(1144, 1);


% Using the formula (1) on page 2, on assignment 2
for i = (1646-502):-1:1
relative_VaR_95(i) = norminv(confidence_level_95) * portfolio_volatility_EWMA(i);
end

for i = (1646-502):-1:1
relative_VaR_99(i) = norminv(confidence_level_99) * portfolio_volatility_EWMA(i);
end

% Måste justera portfoliovikterna?

updatedVaR_99 = norminv(confidence_level_99) * 0.000550568808026982 * portfolio_value;


%% Subtask 1 (c): Calculation of Relative VaR and Expected Shortfall


omega = alpha(i) / portfolio_value;

omega_vector = ones(15, 1) * omega(i);













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
