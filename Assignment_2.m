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
