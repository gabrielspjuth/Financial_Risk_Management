% Assignment 1

%%%%%     Task 1 (a)     %%%%%
dailyDates = readmatrix('Data.xlsx', 'Sheet', 'Daily', 'Range', 'A2:A6032');
dailyOMXS30 = readmatrix('Data.xlsx', 'Sheet', 'Daily', 'Range', 'B2:B6032');
dailyUSDSEK = readmatrix('Data.xlsx', 'Sheet', 'Daily', 'Range', 'C2:C6032');

weeklyDates = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'A2:A1255');
weeklyOMXS30 = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'B2:B1255');
weeklyUSDSEK = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'C2:C1255');

figure(1);
yyaxis left
plot(weeklyDates, weeklyOMXS30)
xlabel('Weeks')
ylabel('OMXS30')

yyaxis right
plot(weeklyDates, weeklyUSDSEK)
xlabel('Weeks')
ylabel('USD/SEK')

for i = 1:(numel(weeklyOMXS30)-1)
    logReturnsOMXS30(i,1) = log(weeklyOMXS30(i+1)) - log(weeklyOMXS30(i));
end

for i = 1:(numel(weeklyUSDSEK)-1)
    logReturnsUSDSEK(i,1) = log(weeklyUSDSEK(i+1)) - log(weeklyUSDSEK(i));
end

figure(2);
yyaxis left
plot(weeklyDates(2:end), logReturnsOMXS30)
ylim([-1 1])
xlabel('Weeks')
ylabel('Log returns for OMXS30')

yyaxis right
plot(weeklyDates(2:end), logReturnsUSDSEK)
xlabel('Weeks')
ylabel('Log returns for USD/SEK')

logReturnsUSDSEKAvg = mean(logReturnsUSDSEK);
logReturnsOMXS30Avg = mean(logReturnsOMXS30);

% Calculating variance and volatility for OMXS30
for i = 1:numel(logReturnsOMXS30)
    returnsOMXS30Squared(1,i)= logReturnsOMXS30(i)^2;
end
varianceOMXS30 = mean(returnsOMXS30Squared) - mean(logReturnsOMXS30Avg)^2;
volatilityOMXS30 = sqrt(varianceOMXS30);

% Calculating variance and volatility for USDSEK
for i = 1:numel(logReturnsUSDSEK)
    returnsUSDSEKSquared(1,i)= logReturnsUSDSEK(i)^2;
end
varianceUSDSEK = mean(returnsUSDSEKSquared) - mean(logReturnsUSDSEKAvg)^2;
volatilityUSDSEK = sqrt(varianceUSDSEK);

z = 1.96;
sampleSize = (numel(logReturnsOMXS30));
% Confidence Interval for OMXS30
SEOMXS30 = volatilityOMXS30 / sqrt(sampleSize);
lowerboundOMSX30 = logReturnsOMXS30Avg - z*(SEOMXS30/sqrt(sampleSize));
upperboundOMXS30 = logReturnsOMXS30Avg + z*(SEOMXS30/sqrt(sampleSize));

% Confidence Interval for USDSEK
SEUSDSEK = volatilityUSDSEK / sqrt(sampleSize);
lowerboundUSDSEK = logReturnsUSDSEKAvg - z*(SEUSDSEK/sqrt(sampleSize));
upperboundUSDSEK = logReturnsUSDSEKAvg + z*(SEUSDSEK/sqrt(sampleSize));

% Calculate annualized average returns
averageAnnualReturnOMXS30 = logReturnsOMXS30Avg * 52 * 100;  % Assuming 52 weeks in a year
averageAnnualReturnUSDSEK = logReturnsUSDSEKAvg * 52 * 100;

% Calculate annualized volatility
annualVolatilityOMXS30 = volatilityOMXS30 * sqrt(52) * 100;
annualVolatilityUSDSEK = volatilityUSDSEK * sqrt(52) * 100;
fprintf('Numerical Results 1a %.2f%%\n')
fprintf('Annualized Average Return for OMXS30: %.2f%%\n', averageAnnualReturnOMXS30);
fprintf('Annualized Average Return for USD/SEK: %.2f%%\n', averageAnnualReturnUSDSEK);

fprintf('Annualized Volatility for OMXS30: %.2f%%\n', annualVolatilityOMXS30);
fprintf('Annualized Volatility for USD/SEK: %.2f%%\n', annualVolatilityUSDSEK);

fprintf('Confidence Interval for OMXS30: [%.8f, %.8f]\n', lowerboundOMSX30, upperboundOMXS30);
fprintf('Confidence Interval for USD/SEK: [%.8f, %.8f]\n', lowerboundUSDSEK, upperboundUSDSEK);


%%%%%     Task 1 (b)     %%%%%

for i = 1:(numel(dailyOMXS30)-1)
    dailyLogReturnsOMXS30(i,1) = log(dailyOMXS30(i+1)) - log(dailyOMXS30(i));
end

for i = 1:(numel(dailyUSDSEK)-1)
    dailyLogReturnsUSDSEK(i,1) = log(dailyUSDSEK(i+1)) - log(dailyUSDSEK(i));
end

dailyLogReturnsUSDSEKAvg = mean(dailyLogReturnsUSDSEK);
dailyLogReturnsOMXS30Avg = mean(dailyLogReturnsOMXS30);

% Calculating variance and volatility for OMXS30 daily
for i = 1:numel(dailyLogReturnsOMXS30)
    dailyreturnsOMXS30Squared(1,i)= dailyLogReturnsOMXS30(i)^2;
end

varianceOMXS30daily = mean(dailyreturnsOMXS30Squared) - mean(dailyLogReturnsOMXS30Avg)^2;
volatilityOMXS30daily = sqrt(varianceOMXS30daily);

% Calculating variance and volatility for USDSEK daily

for i = 1:numel(dailyLogReturnsOMXS30)
    dailyreturnsUSDSEKSquared(1,i)= dailyLogReturnsUSDSEK(i)^2;
end

varianceUSDSEKdaily = mean(dailyreturnsUSDSEKSquared) - mean(dailyLogReturnsUSDSEKAvg)^2;
volatilityUSDSEKdaily = sqrt(varianceUSDSEKdaily);

for i = 1:numel(logReturnsOMXS30)
    returnsOMXS30Squared(1,i)= logReturnsOMXS30(i)^2;
end
varianceOMXS30 = mean(returnsOMXS30Squared) - mean(logReturnsOMXS30Avg)^2;
volatilityOMXS30 = sqrt(varianceOMXS30);



% i. Calculating Skewness and (excess) kurtosis weekly
weeklyOMXS30Skewness = skewness(logReturnsOMXS30);
weeklyUSDSEKSkewness = skewness(logReturnsUSDSEK);

weeklyOMXS30kurtosis = kurtosis(logReturnsOMXS30);
weeklyUSDSEKSkurtosis = kurtosis(logReturnsUSDSEK);

% i. Calculating Skewness and (excess) kurtosis daily
dailyOMXS30Skewness = skewness(dailyLogReturnsOMXS30);
dailyUSDSEKSkewness = skewness(dailyLogReturnsUSDSEK);

dailyOMXS30kurtosis = kurtosis(dailyLogReturnsOMXS30);
dailyUSDSEKSkurtosis = kurtosis(dailyLogReturnsUSDSEK);

%ii. Histograms
% Assuming logReturnsOMXS30, logReturnsUSDSEK, dailyLogReturnsOMXS30, and dailyLogReturnsUSDSEK are defined

% Create a 2x2 subplot
figure(3);

% Plot the first histogram in the first subplot
subplot(2, 2, 1);
histfit(logReturnsOMXS30);
title('Weekly OMXS30');

% Plot the second histogram in the second subplot
subplot(2, 2, 2);
histfit(logReturnsUSDSEK);
title('Weekly USDSEK');

% Plot the third histogram in the third subplot
subplot(2, 2, 3);
histfit(dailyLogReturnsOMXS30);
title('Daily OMXS30');

% Plot the fourth histogram in the fourth subplot
subplot(2, 2, 4);
histfit(dailyLogReturnsUSDSEK);
title('Daily USDSEK');


% Calculate percentiles for the weekly OMXS30 log returns
percentilesOMXS30 = prctile(logReturnsOMXS30, [1, 5, 95, 99]);

% Calculate percentiles for the weekly USD/SEK log returns
percentilesUSDSEK = prctile(logReturnsUSDSEK, [1, 5, 95, 99]);

% Calculate percentiles for the daily OMXS30 log returns
percentilesDailyOMXS30 = prctile(dailyLogReturnsOMXS30, [1, 5, 95, 99]);

% Calculate percentiles for the daily USD/SEK log returns
percentilesDailyUSDSEK = prctile(dailyLogReturnsUSDSEK, [1, 5, 95, 99]);

%iii. QQ Plots

%QQ Plot weeklyOMXS30

QuantilevectorOMXS30weekly = zeros(numel(logReturnsOMXS30),1);
Q1 = linspace(0,1, numel(logReturnsOMXS30));
Q2 = mean(logReturnsOMXS30);
Q3 = volatilityOMXS30;

weeklyQQOMXS30(:, 1) = logReturnsOMXS30(:, 1);
weeklyQQOMXS30(:, 1) = sort(weeklyQQOMXS30(:, 1));
weeklyQQOMXS30(:, 2) = norminv(Q1, Q2, Q3);
figure(4);
subplot(2, 2, 1)
scatter(weeklyQQOMXS30(:, 2), weeklyQQOMXS30(:, 1), '+');
hold on;
plot([min(weeklyQQOMXS30(:, 1)), max(weeklyQQOMXS30(:, 1))], [min(weeklyQQOMXS30(:, 1)), max(weeklyQQOMXS30(:, 1))], 'r--');
xlim([-4*volatilityOMXS30 4*volatilityOMXS30]) %Fixa sen
xticks([-4*volatilityOMXS30 -2*volatilityOMXS30 0*volatilityOMXS30 2*volatilityOMXS30 4*volatilityOMXS30])
xticklabels({'x=-4', 'x=-2', 'x=0', 'x=2', 'x=4'})
xlabel('Standard Deviations');
ylabel('Log Returns for USD/SEK');
hold off;


Q1 = linspace(0, 1, numel(logReturnsUSDSEK));
Q2 = mean(logReturnsUSDSEK);
Q3 = volatilityUSDSEK;

weeklyQQUSDSEK(:, 1) = logReturnsUSDSEK(:, 1);
weeklyQQUSDSEK(:, 1) = sort(weeklyQQUSDSEK(:, 1));
weeklyQQUSDSEK(:, 2) = norminv(Q1, Q2, Q3);

% Add an offset to the y-values
yOffset = 0.1; % Adjust the offset as needed
offsetY = weeklyQQUSDSEK(:, 1) + yOffset;

subplot(2, 2, 2)
scatter(weeklyQQUSDSEK(:, 2), offsetY, '+');
hold on;
plot([min(weeklyQQUSDSEK(:, 1)), max(weeklyQQUSDSEK(:, 1))], [min(offsetY), max(offsetY)], 'r--');
xlim([-4*volatilityUSDSEK 4*volatilityUSDSEK])
xticks([-4*volatilityUSDSEK -2*volatilityUSDSEK 0 2*volatilityUSDSEK 4*volatilityUSDSEK]);
xticklabels({'x=-4', 'x=-2', 'x=0', 'x=2', 'x=4'});
xlabel('Standard Deviations');
ylabel('Log Returns for USD/SEK');
hold off;

% Daily Log Returns
for i = 1:(numel(dailyOMXS30)-1)
    dailyLogReturnsOMXS30(i, 1) = log(dailyOMXS30(i+1)) - log(dailyOMXS30(i));
end

for i = 1:(numel(dailyUSDSEK)-1)
    dailyLogReturnsUSDSEK(i, 1) = log(dailyUSDSEK(i+1)) - log(dailyUSDSEK(i));
end

% Daily Log Return Statistics
dailyLogReturnsOMXS30Avg = mean(dailyLogReturnsOMXS30);
dailyLogReturnsUSDSEKAvg = mean(dailyLogReturnsUSDSEK);

% Calculate variance and volatility for daily OMXS30
varianceOMXS30daily = var(dailyLogReturnsOMXS30);
volatilityOMXS30daily = sqrt(varianceOMXS30daily);

% Calculate variance and volatility for daily USD/SEK
varianceUSDSEKdaily = var(dailyLogReturnsUSDSEK);
volatilityUSDSEKdaily = sqrt(varianceUSDSEKdaily);

% Create QQ plots for daily returns

% QQ Plot for daily OMXS30
subplot(2, 2, 3)
QuantilevectorDailyOMXS30 = zeros(numel(dailyLogReturnsOMXS30), 1);
Q1DailyOMXS30 = linspace(0, 1, numel(dailyLogReturnsOMXS30));
Q2DailyOMXS30 = mean(dailyLogReturnsOMXS30);
Q3DailyOMXS30 = volatilityOMXS30daily;

dailyQQOMXS30(:, 1) = dailyLogReturnsOMXS30(:, 1);
dailyQQOMXS30(:, 1) = sort(dailyQQOMXS30(:, 1));
dailyQQOMXS30(:, 2) = norminv(Q1DailyOMXS30, Q2DailyOMXS30, Q3DailyOMXS30);

scatter(dailyQQOMXS30(:, 2), dailyQQOMXS30(:, 1), '+');
hold on;
plot([min(dailyQQOMXS30(:, 1)), max(dailyQQOMXS30(:, 1))], [min(dailyQQOMXS30(:, 1)), max(dailyQQOMXS30(:, 1))], 'r--');
xlim([-4 * volatilityOMXS30daily, 4 * volatilityOMXS30daily]);
xticks([-4 * volatilityOMXS30daily, -2 * volatilityOMXS30daily, 0, 2 * volatilityOMXS30daily, 4 * volatilityOMXS30daily]);
xticklabels({'x=-4', 'x=-2', 'x=0', 'x=2', 'x=4'});
xlabel('Standard Deviations');
ylabel('Log Returns for OMXS30');
title('QQ Plot for Daily OMXS30 Returns');
hold off;

% QQ Plot for daily USD/SEK
subplot(2, 2, 4)
QuantilevectorDailyUSDSEK = zeros(numel(dailyLogReturnsUSDSEK), 1);
Q1DailyUSDSEK = linspace(0, 1, numel(dailyLogReturnsUSDSEK));
Q2DailyUSDSEK = mean(dailyLogReturnsUSDSEK);
Q3DailyUSDSEK = volatilityUSDSEKdaily;

dailyQQUSDSEK(:, 1) = dailyLogReturnsUSDSEK(:, 1);
dailyQQUSDSEK(:, 1) = sort(dailyQQUSDSEK(:, 1));
dailyQQUSDSEK(:, 2) = norminv(Q1DailyUSDSEK, Q2DailyUSDSEK, Q3DailyUSDSEK);

scatter(dailyQQUSDSEK(:, 2), dailyQQUSDSEK(:, 1), '+');
hold on;
plot([min(dailyQQUSDSEK(:, 1)), max(dailyQQUSDSEK(:, 1))], [min(dailyQQUSDSEK(:, 1)), max(dailyQQUSDSEK(:, 1))], 'r--');
xlim([-4 * volatilityUSDSEKdaily, 4 * volatilityUSDSEKdaily]);
xticks([-4 * volatilityUSDSEKdaily, -2 * volatilityUSDSEKdaily, 0, 2 * volatilityUSDSEKdaily, 4 * volatilityUSDSEKdaily]);
xticklabels({'x=-4', 'x=-2', 'x=0', 'x=2', 'x=4'});
xlabel('Standard Deviations');
ylabel('Log Returns for USD/SEK');
title('QQ Plot for Daily USD/SEK Returns');
hold off;

fprintf('Numerical Results 1b\n');

% Skewness and Kurtosis for Weekly OMXS30
fprintf('Skewness for Weekly OMXS30: %.4f\n', weeklyOMXS30Skewness);
fprintf('Kurtosis for Weekly OMXS30: %.4f\n', weeklyOMXS30kurtosis);

% Skewness and Kurtosis for Weekly USD/SEK
fprintf('Skewness for Weekly USD/SEK: %.4f\n', weeklyUSDSEKSkewness);
fprintf('Kurtosis for Weekly USD/SEK: %.4f\n', weeklyUSDSEKSkurtosis);

% Skewness and Kurtosis for Daily OMXS30
fprintf('Skewness for Daily OMXS30: %.4f\n', dailyOMXS30Skewness);
fprintf('Kurtosis for Daily OMXS30: %.4f\n', dailyOMXS30kurtosis);

% Skewness and Kurtosis for Daily USD/SEK
fprintf('Skewness for Daily USD/SEK: %.4f\n', dailyUSDSEKSkewness);
fprintf('Kurtosis for Daily USD/SEK: %.4f\n', dailyUSDSEKSkurtosis);

% Define an array for percentiles and initialize it
percentilesArray = zeros(4, 4);

% Percentiles for Weekly OMXS30 log returns
percentilesArray(1, 1) = 1;
percentilesArray(2, 1) = 5;
percentilesArray(3, 1) = 95;
percentilesArray(4, 1) = 99;

for i = 1:4
    percentilesArray(i, 2) = round(percentilesOMXS30(i), 2);
end

% Percentiles for Weekly USD/SEK log returns
percentilesArray(1, 3) = 1;
percentilesArray(2, 3) = 5;
percentilesArray(3, 3) = 95;
percentilesArray(4, 3) = 99;

for i = 1:4
    percentilesArray(i, 4) = round(percentilesUSDSEK(i), 2);
end

% Percentiles for Daily OMXS30 log returns
percentilesArray(1, 5) = 1;
percentilesArray(2, 5) = 5;
percentilesArray(3, 5) = 95;
percentilesArray(4, 5) = 99;

for i = 1:4
    percentilesArray(i, 6) = round(percentilesDailyOMXS30(i), 2);
end

% Percentiles for Daily USD/SEK log returns
percentilesArray(1, 7) = 1;
percentilesArray(2, 7) = 5;
percentilesArray(3, 7) = 95;
percentilesArray(4, 7) = 99;

for i = 1:4
    percentilesArray(i, 8) = round(percentilesDailyUSDSEK(i), 2);
end

% Display the 4x4 percentiles matrix with improved alignment for headings
fprintf('Percentiles Matrix (4x4) with 2 significant digits:\n');
fprintf('  Weekly OMXS30        Weekly USD/SEK        Daily OMXS30      Daily USD/SEK\n');
disp(percentilesArray);
format short; % Reset format to default


% Create a scatter plot
figure(25)
qqplot(logReturnsOMXS30);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Load data
weeklyDates = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'A2:A1255');
weeklyOMXS30 = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'B2:B1255');
weeklyUSDSEK = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'C2:C1255');

% Task 2(a)

% Calculate log returns
logReturnsOMXS30 = diff(log(weeklyOMXS30));
logReturnsUSDSEK = diff(log(weeklyUSDSEK));

% Define rolling window sizes
window_30w = 30;
window_90w = 90;

N = 1223;
logOMXS30_30 = logReturnsOMXS30(31:end);
logUSDSEK_30 = logReturnsUSDSEK(31:end);
logUSDSEK_90 = logReturnsUSDSEK(91:end);
logOMXS30_90 = logReturnsOMXS30(91:end);

% EWMA Rolling Volatility


lambda_30 = 2 / (30 + 1);  % EWMA parameter for 30-week window
lambda_90 = 2 / (90 + 1);  % EWMA parameter for 90-week window

% Initialize arrays for EWMA volatilities
ewma_volatility_omxs30_30 = zeros(size(logReturnsOMXS30));
ewma_volatility_usdsek_30 = zeros(size(logReturnsUSDSEK));

ewma_volatility_omxs30_90 = zeros(size(logReturnsOMXS30));
ewma_volatility_usdsek_90 = zeros(size(logReturnsUSDSEK));

% Calculate EWMA volatilities
for t = 2:length(logReturnsOMXS30)
    ewma_volatility_omxs30_30(t) = sqrt(lambda_30 * (logReturnsOMXS30(t)^2) + (1 - lambda_30) * ewma_volatility_omxs30_30(t - 1)^2);
    ewma_volatility_usdsek_30(t) = sqrt(lambda_30 * (logReturnsUSDSEK(t)^2) + (1 - lambda_30) * ewma_volatility_usdsek_30(t - 1)^2);

    ewma_volatility_omxs30_90(t) = sqrt(lambda_90 * (logReturnsOMXS30(t)^2) + (1 - lambda_90) * ewma_volatility_omxs30_90(t - 1)^2);
    ewma_volatility_usdsek_90(t) = sqrt(lambda_90 * (logReturnsUSDSEK(t)^2) + (1 - lambda_90) * ewma_volatility_usdsek_90(t - 1)^2);
end



% Calculate 30-week rolling average for OMXS30
moving_average_omxs30_30w = movmean(logReturnsOMXS30, [0 29]);
moving_average_omxs30_30w = moving_average_omxs30_30w(1: end-30);
historicalVolOMXS30_30w = sqrt((logOMXS30_30 - moving_average_omxs30_30w).^2) * sqrt(252);

% Calculate 30-week rolling average for USDSEK
moving_average_usdsek_30w = movmean(logReturnsUSDSEK, [0 29]);
moving_average_usdsek_30w = moving_average_usdsek_30w(1: end-30);
historicalVolUSDSEK_30w = sqrt((logUSDSEK_30 - moving_average_usdsek_30w).^2) * sqrt(252);

% Calculate 90-week rolling average for OMXS30
moving_average_omxs30_90w = movmean(logReturnsOMXS30, [0 89]);
moving_average_omxs30_90w = moving_average_omxs30_90w(1: end-90);
historicalVolOMXS30_90w = sqrt((logOMXS30_90 - moving_average_omxs30_90w).^2) * sqrt(252);

% Calculate 30-week rolling average and 90-week rolling average for USDSEK
moving_average_usdsek_90w = movmean(logReturnsUSDSEK, [0 89]);
moving_average_usdsek_90w = moving_average_usdsek_90w(1: end-90);
historicalVolUSDSEK_90w = sqrt((logUSDSEK_90 - moving_average_usdsek_90w).^2) * sqrt(252);


% Subplot 1

figure(1);
subplot(2, 2, 1);
plot(moving_average_omxs30_30w);
title('OMXS30 30-week Rolling Average');

% Subplot 2
subplot(2, 2, 2);
plot(moving_average_usdsek_30w);
title('USDSEK 30-week Rolling Average');

% Subplot 3
subplot(2, 2, 3);
plot(moving_average_omxs30_90w);
title('OMXS30 90-week Rolling Average');

% Subplot 4
subplot(2, 2, 4);
plot(moving_average_usdsek_90w);
title('USDSEK 90-week Rolling Average');


% Task 2 (b)

lambda = 0.94;

OMXS30EWMA = zeros(numel(logReturnsOMXS30), 1);

for i = 1:(numel(logReturnsOMXS30))
    OMXS30EWMA(i) = sqrt((1 - lambda) * logReturnsOMXS30(i)^2);
end


%%%%%%%%%%%%%%%% Task 2 (d) %%%%%%%%%%%%%%

weeklyDatesGARCH = readmatrix('GARCH_Solver_Complete.xlsx', 'Sheet', 'OMXS30', 'Range', 'A3:A1256');
OMXS30GARCHVolatility = readmatrix('GARCH_Solver_Complete.xlsx', 'Sheet', 'OMXS30', 'Range', 'M4:M1256');
USDSEKGARCHVolatility = readmatrix('GARCH_Solver_Complete.xlsx', 'Sheet', 'USDSEK', 'Range', 'M4:M1256');

avgOMXS30Vol = mean(OMXS30GARCHVolatility(2, end));
avgUSDSEKVol = mean(USDSEKGARCHVolatility(2, end));

standardizedReturnsOMXS30 = logReturnsOMXS30 / (avgOMXS30Vol * sqrt(52));
standardizedReturnsUSDSEK = logReturnsUSDSEK / avgUSDSEKVol* sqrt(52);

figure(2);
qqplot(standardizedReturnsOMXS30);

figure(3);
qqplot(standardizedReturnsUSDSEK);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% TASK 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

weeklyDates = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'A2:A1255');
weeklyOMXS30 = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'B2:B1255');
weeklyUSDSEK = readmatrix('Data.xlsx', 'Sheet', 'Weekly', 'Range', 'C2:C1255');
normcdfOMXS30 = readmatrix('GARCH_Solver_Complete.xlsx', 'Sheet', 'Variance target OMXS30', 'Range', 'Q5:Q1256');
normcdfUSDSEK = readmatrix('GARCH_Solver_Complete.xlsx', 'Sheet', 'Variance target USDSEK', 'Range', 'P5:P1256');
residualUSDSEK = readmatrix('GARCH_Solver_Complete.xlsx', 'Sheet', 'Variance target USDSEK', 'Range', 'J4:J1256');
residualOMXS30 = readmatrix('GARCH_Solver_Complete.xlsx', 'Sheet', 'Variance target OMXS30', 'Range', 'K5:K1256');

%%%%%%%%%%%      Task 3 (a)        %%%%%%%%%%%%
 
 for i = 1:(numel(weeklyOMXS30)-1)
     logReturnsOMXS30(i,1) = log(weeklyOMXS30(i+1)) - log(weeklyOMXS30(i));
 end
 
 for i = 1:(numel(weeklyUSDSEK)-1)
     logReturnsUSDSEK(i,1) = log(weeklyUSDSEK(i+1)) - log(weeklyUSDSEK(i));
 end

corrOMXS30USDSEK = corr(logReturnsOMXS30, logReturnsUSDSEK);

%% Task 3 (b) %%
% Compute autocorrelation for lags 1 to 5 weeks
maxLag = 5;

autocorr_valuesOMXS30 = autocorr(logReturnsOMXS30, maxLag)';
autocorr_valuesUSDSEK = autocorr(logReturnsUSDSEK, maxLag)';

% Create a matrix with two rows displaying autocorrelations
autocorr_matrix = [autocorr_valuesOMXS30; autocorr_valuesUSDSEK];

% Display the autocorrelation matrix
disp('Autocorrelation matrix for lags 1 to 5 weeks (OMXS30 and USDSEK):');
disp(autocorr_matrix(:, 2:end));


%% Task 3 (c) %%



U = [normcdfOMXS30, normcdfUSDSEK];

% Determining the best copula
copula_types = {'t', 'gumbel', 'clayton', 'frank', 'Gaussian'};
best_log_likelihood = -Inf;
best_copula_type = '';
best_copula_params = [];

[rhot, nut] = copulafit("t", U);
 alphagumbel = copulafit("gumbel", U);
 alphaclayton = copulafit("clayton", U);
 alphafrank = copulafit("frank", U);
 rhogaussian = copulafit("Gaussian", U);

for i = 1:length(copula_types)
    copula_type = copula_types{i};
    copula_params = copulafit(copula_type, U);

       % Extracting the degrees of freedom parameter for the t copula
    if strcmp(copula_type, 't')
        log_likelihood = sum(log(copulapdf('t', U, rhot, nut)));
        ll_t = log_likelihood;
    elseif strcmp(copula_type, 'gumbel')
        log_likelihood = sum(log(copulapdf('gumbel', U, alphagumbel)));
        ll_gumbel = log_likelihood;
    elseif strcmp(copula_type, 'clayton')
        log_likelihood = sum(log(copulapdf('clayton', U, alphaclayton)));
        ll_clayton = log_likelihood;
    elseif strcmp(copula_type, 'frank')
        log_likelihood = sum(log(copulapdf('frank', U, alphafrank)));
        ll_frank = log_likelihood;
    elseif strcmp(copula_type, 'Gaussian')
        log_likelihood = sum(log(copulapdf('Gaussian', U, rhogaussian)));
        ll_gaussian = log_likelihood;
    end
    
    
    fprintf('Log Likelihood for %s: %f\n', copula_type, log_likelihood);

    if log_likelihood > best_log_likelihood
        best_log_likelihood = log_likelihood;
        best_copula_type = copula_type;
    end
end


disp(['Best Copula Type: ', best_copula_type]);


% Best Copula Type
best_copula_type = 't';

% Best Copula Parameters
[rhot, nut] = copulafit(best_copula_type, U);

fprintf('Best Copula Type: %s\n', best_copula_type);

% Generating samples from the best copula
num_samples = 1252;
copula_samples = copularnd(best_copula_type, rhot, nut, num_samples);

% Plotting the time series and copula samples
figure;

% Plotting the time series as a copula
subplot(1, 4, 1);
scatter(U(:, 1), U(:, 2), 'r.');
title('Time Series as a Copula');
xlabel('U1 (OMXS)');
ylabel('U2 (USD/SEK)');
axis square;
grid on;

% Plotting the copula samples for the best copula
subplot(1, 4, 2);
copula_samples_best = copularnd("t", rhot, nut, 1252);
scatter(copula_samples_best(:, 1), copula_samples_best(:, 2), 'b.');
title(['Best Copula Type: ', best_copula_type]);
xlabel('U1');
ylabel('U2');
axis square;
grid on;

sgtitle(['Time Series and Copula Samples - Best Copula Type: ', best_copula_type]);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% F PRINT %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

output.RIC = {'OMXS30', 'USDSEK'};

% stat
output.stat.mu = [logReturnsOMXS30Avg, logReturnsUSDSEKAvg];
output.stat.sigma = [volatilityOMXS30, volatilityUSDSEK];
output.stat.CI = [100*lowerboundOMSX30 100*upperboundOMXS30; 100*lowerboundUSDSEK 100*upperboundUSDSEK];
output.stat.skew = [weeklyOMXS30Skewness, weeklyUSDSEKSkewness, dailyOMXS30Skewness, dailyUSDSEKSkewness];
output.stat.kurt = [weeklyOMXS30kurtosis, weeklyUSDSEKSkurtosis, dailyOMXS30kurtosis, dailyUSDSEKSkurtosis];
output.stat.perc = [percentilesOMXS30; percentilesUSDSEK; percentilesDailyOMXS30; percentilesDailyUSDSEK]; %OK?
output.stat.corr = [corr(logReturnsOMXS30, logReturnsUSDSEK)];
output.stat.acorr = [transpose(autocorr_matrix(:, 2:end))]; % ok?

% EWMA
output.EWMA.obj = [7712.515833, 9069.995795];
output.EWMA.param = [0.923972232, 0.790284312];

% GARCH
output.GARCH.obj = [7775.26729813614000, 9236.685874894];
output.GARCH.param = [0.0000375833253363329, 0.146938556451542, 0.818102850826723, 0.0000182019784549581, 0.0800291461030378, 0.845082435788515];
output.GARCH.objVT = [7775.2672980646500, 9236.68587493];
output.GARCH.paramVT = [0.0000375835289169586, 0.14694347278145, 0.818099068971663, 0.0000183111873631994, 0.0800275039636921, 0.84508361338695];

output.copulaLogL = [ll_gaussian ll_t ll_gumbel ll_clayton ll_frank];

printResults(output, true);
