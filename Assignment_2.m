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

fprintf('Optimized beta: %f\n', paramsOptimized(1,1));
fprintf('Optimized ksi: %f\n', paramsOptimized(2,1));

%% Subtask 2 (b): Estimating Parameters for a Volatile Period (approx. 5 years of data)

% Calculating Logarithmic Returns for each Stock
% and the average for the complete portfolio
dates = 1647;
stocks = 15;

prices = S(:, 2:end);
pricesInOrder = flipud(prices); % Row 1 = first date, Row 1646 = today

logReturns = zeros(dates-1, stocks);

for j = 1:stocks
    for i = 1:(dates-1)
        logReturns(i, j) = log(pricesInOrder(i + 1, j) / pricesInOrder(i, j));
    end
end

% Check this out - the time period depends whether we flip it or not
%logReturnsPortfolio = flipud(sum(logReturns, 2) ./ stocks);

logReturnsPortfolio = sum(logReturns, 2) ./ stocks;

logAverageReturnPortfolio = mean(logReturnsPortfolio);

% Finding the 5-year period with the highest average volatility

% Calculating the number of weeks in a 5 year period
yearsInWeeks = 5 * 52;

% Creating a new vector with the sample std dev for each time window
volatilitySampleVariance = zeros((dates-1)-yearsInWeeks, 1);

% Calculating the sample standard deviation for 5 year time windows
% The outer for loop iterates over all possible 5 year periods

for i = 1:((dates-1)-yearsInWeeks)
    avgSqVol = 0;
    for j = 1:yearsInWeeks
        avgSqVol = avgSqVol + ((logReturnsPortfolio((i+j), 1) - logAverageReturnPortfolio)^2);
    end
    volatilitySampleVariance(i) = sqrt(avgSqVol / (yearsInWeeks-1));
end

% Finding the period with highest average volatility
% and the corresponding index
[maximumVolatility, timeSlot] = max(volatilitySampleVariance);

timePeriodMaxVol = [dates-timeSlot-yearsInWeeks, dates-timeSlot];

% Estimating the parameters
portfolioReturns = transpose(portfolio_returns);
portfolioReturnsWithHighestVolatility = portfolioReturns(timePeriodMaxVol(1):timePeriodMaxVol(2));
portfolioReturnsWithHighestVolatilitySorted = sort(portfolioReturnsWithHighestVolatility, 'ascend');


% Calculating the threshold level in the 95th percentile
nWithHighestVolatility = 0.05 * yearsInWeeks;
uWithHighestVolatility = portfolioReturnsWithHighestVolatilitySorted(nWithHighestVolatility, 1);

% Calculating y for the 13 biggest losses in the distribution 
% of the time window
yWithHighestVolatility = portfolioReturnsWithHighestVolatilitySorted(1:nWithHighestVolatility-1) - u;

% Definining the log likelihood function
fWithHighestVolatility = @(params) -sum(log(1/params(1))*(1+params(2)*(yWithHighestVolatility/params(1))).^((-1/params(2))-1));

% Start values for ksi and beta
beta = 0.15;
ksi = 0.4;
params = [beta; ksi];

% Optimizing the parameters
paramsOptimizedWithHighestVolatility = fmincon(fWithHighestVolatility, params);

fprintf('Optimized beta (5-year highest vol): %f\n', paramsOptimizedWithHighestVolatility(1,1));
fprintf('Optimized ksi (5-year highest vol) %f\n', paramsOptimizedWithHighestVolatility(2,1));



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Task 3: Risk Factor Mapping
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Importing risk factors from March 5th, 2008 to January 10th, 2022 
SPX = readmatrix('timeSeries.xlsx', 'Sheet', 'Problem 3', 'Range', 'C4:C3429');
VIX = readmatrix('timeSeries.xlsx', 'Sheet', 'Problem 3', 'Range', 'D4:D3429');
USD3MFSR = readmatrix('timeSeries.xlsx', 'Sheet', 'Problem 3', 'Range', 'E4:E3429');

% Calculating time to maturity for the two different option contract types
timeToMaturity_MarchContract = days252bus('2022-01-10', '2022-03-18');
timeToMaturity_AprilContract = days252bus('2022-01-10', '2022-04-15');

% Black-Scholes model requires an annualized continuous compounded rate

annualRiskFreeRate = (log(1+USD3MFSR(1)*(3/12))) / (3/12);


% Calculating the Option Prices at t = T +1 with S&P500 as underlying
% Parameters from timeSeries.xlsk. Implied volatility is mean of bid, ask
Mar22Call = BlackScholes(SPX(1), 4700, 15.77, timeToMaturity_MarchContract, 0, annualRiskFreeRate, 'call');
Mar22Put = BlackScholes(SPX(1), 4600, 18.28, timeToMaturity_MarchContract, 0, annualRiskFreeRate, 'put');
Apr22Call = BlackScholes(SPX(1), 4750, 16.25, timeToMaturity_AprilContract, 0, annualRiskFreeRate, 'call');

% Calculating the Gradient Vectors, using los greekazos
DeltaMar22Call = Delta(SPX(1), 4700, 15.77, timeToMaturity_MarchContract, 0, annualRiskFreeRate, 0.05, 'call');
DeltaMar22Put = Delta(SPX(1), 4600, 18.28, timeToMaturity_MarchContract, 0, annualRiskFreeRate, 0.05, 'put');
DeltaApr22Call = Delta(SPX(1), 4750, 16.25, timeToMaturity_AprilContract, 0, annualRiskFreeRate, 0.05, 'call');

VegaMar22Call = Vega(SPX(1), 4700, 15.77, timeToMaturity_MarchContract, 0, annualRiskFreeRate, 0.05);
VegaMar22Put = Vega(SPX(1), 4600, 18.28, timeToMaturity_MarchContract, 0, annualRiskFreeRate, 0.05);
VegaApr22Call = Vega(SPX(1), 4750, 16.25, timeToMaturity_AprilContract, 0, annualRiskFreeRate, 0.05);

RhoMar22Call = Rho(SPX(1), 4700, 15.77, timeToMaturity_MarchContract, 0, annualRiskFreeRate, 'call');
RhoMar22Put = Rho(SPX(1), 4600, 18.28, timeToMaturity_MarchContract, 0, annualRiskFreeRate, 'put');
RhoApr22Call = Rho(SPX(1), 4750, 16.25, timeToMaturity_AprilContract, 0, annualRiskFreeRate, 'call');


gradientMar22Call = [DeltaMar22Call, VegaMar22Call, RhoMar22Call];
gradientMar22Put = [DeltaMar22Put, VegaMar22Put, RhoMar22Put];
gradientApr22Call = [DeltaApr22Call, VegaApr22Call, RhoApr22Call];


% Function for Black-Scholes. 
function res = BlackScholes(S, K, volatility, T, t, rate, type)
    d = dBSM(S, K, volatility, T, t, rate);
    if (type == "call")
        res = S*normcdf(d(1)) - K*exp(-rate*(T-t))*normcdf(d(2));
    elseif (type == "put")
        res = K*exp(-rate*(T-t))*normcdf(-d(2)) - S*normcdf(-d(1));
    end
end

% Function for calculating Delta
function res = Delta(S, K, volatility, T, t, rate, q, type)
    d = dBSM(S, K, volatility, T, t, rate);
    if (type == "call") 
        res = exp(-q*(T-t))*normcdf(d(1));
    elseif (type == "put")
        res = exp(-q*(T-t))*(normcdf(d(1)) - 1); % Milken Finance Institute kommer formeln ifrån
    end
end

% Funtion for calculating vega
function res = Vega(S, K, volatility, T, t, rate, q)
    d = dBSM(S, K, volatility, T, t, rate);
    res = S*exp(-q*(T-t))*normcdf(d(1))*sqrt(T-t);
end

% Function for calculating rho
function res = Rho(S, K, volatility, T, t, rate, type)
    d = dBSM(S, K, volatility, T, t, rate);
    if (type == "call")
        res = K*(T-t)*exp(-rate*(T-t))*normcdf(d(2));
    elseif (type == "put")
        res = -K*(T-t)*exp(-rate*(T-t))*normcdf(-d(2));
    end
end


% Function for calculating d1 and d2
function res = dBSM(S, K, volatility, T, t, rate)
    d1 = (log(S/K) + (rate+((volatility)^2)/2)*(T-t)) / (volatility * sqrt(T-t));
    d2 = d1 - (volatility*sqrt(T-t));
    res = [d1, d2];
end


















