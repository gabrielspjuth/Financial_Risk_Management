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

portfolio_value = 10000000;

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

VaR_975 = norminv(confidence_level_975) * portfolio_volatiliy * portfolio_value;

VaR_99 = norminv(confidence_level_99) * portfolio_volatiliy * portfolio_value;

% Output in Command Window
fprintf('VaR at %d%% confidence level: SEK %.2f\n', confidence_level_95 * 100, VaR_95);
fprintf('VaR at %.1f%% confidence level: SEK %.2f\n', confidence_level_975 * 100, VaR_975);
fprintf('VaR at %d%% confidence level: SEK %.2f\n', confidence_level_99 * 100, VaR_99);

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

% updatedVaR_99 = norminv(confidence_level_99) * 0.000550568808026982 * portfolio_value;


%% Subtask 1 (c): Calculation of Relative VaR and Expected Shortfall


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 2: Extreme Value Theory
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Subtask 2 (a): Estimating Parameters

% Sorting the portfolio returns in ascending order
sortedPortfolioReturns = sort(transpose(portfolio_returns), 'ascend');

% Calculating number of observations which fall in the 95th percentile
n = 1646;
nu = round(0.05 * n);

% Calculating the threshold level, u (95th percentile of the loss)
u = sortedPortfolioReturns(nu, 1);

% Calculating y for the 81 biggest losses in the distribution
y = sortedPortfolioReturns(1:nu-1) - u;

% Definining the log likelihood function
f = @(params) -sum(log(1/params(1))*(1+params(2)*(y/params(1))).^((-1/params(2))-1));

% Start values for ksi and beta
beta = 0.15;
ksi = 0.4;
params = [beta; ksi];

% Optimizing the parameters
paramsOptimized = fmincon(f, params);

% Calculating VaR with 95% confidence level based on the optimized params
c = 0.99;
VaR_EVT =   u+(paramsOptimized(1)/paramsOptimized(2))*((n/(nu-1)*(1-c))^(-params(2))-1);

%% Subtask 2 (b): 
