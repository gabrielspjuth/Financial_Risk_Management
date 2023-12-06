function printResults(output,printToFile)

%output - structure array (help struct)
%%% Example %%%
% output.RIC = {'.OMXS30', 'ERICb:ST'}
% output.stat.
%        mu: [6.12 7.14]        Mean returns (weekly data)
%     sigma: [15.3 23.1]        St.dev. (weekly data)
%        CI: [2x2 double]       Conf. int., [cl1 cu1; cl2 cu2], cl1 = lower limit RIC1
%      skew: [0 0 0 0]          Coeff. of skewness [RIC1w, RIC2w, RIC1d, RIC2d]
%      kurt: [0 2.2 1 7.2]      Coeff. of kurt [RIC1w, RIC2w, RIC1d, RIC2d]
%      perc: [4x4 double]       Percentiles [pw11 pw15 pw195 pw199; pw21 pw25 pw295 pw299;], pw11 weekly data, RIC1, first percentile etc.
%      corr: [0.5 0.45 0.4]     Correlations [Linear Spearman Kendall]
%     acorr: [5x2]              Autocorrelations, rows: lag 1-5, cols: RIC1-2
% output.EWMA.
%       obj: [100 101]          [log-L (RIC1), log-L (RIC2)]
%     param: [0.9400 0.9600]    [lambda (RIC1), lambda (RIC2)]
% output.GARCH.
%       obj: [1x2 double]       [log-L (RIC1) unconstrained, log-L (RIC2) unconstrained]
%     param: [1x6 double]       [sigma, alpha, beta (RIC1), sigma, alpha, beta (RIC2) (unconstrained MLE)] %sigma is the yearly volatility, i.e. sqrt(VL*52), from MLE
% output.GARCH.
%       objVT: [1x2 double]     [log-L (RIC1) variance targeting, log-L (RIC2) variance targeting.;]
%     paramVT: [1x6 double]     [sigma, alpha, beta (RIC1), sigma, alpha, beta (RIC2) (variance targeting)] %sigma is the yearly volatility, i.e. sqrt(VL*52), from sample estimate.
% output.copulaLogL = [1 2 3 4 5] Normal, Stud.-t, Gumbel, Clayton, Frank
% Assigning a variable to the struct type: output.stat.mu = [6.12 7.14];

if printToFile
    fileID = fopen('NumResults.txt','w+');
    fprintf(fileID,'%s\r\n', 'PROBLEM 1');
    fprintf(fileID,'%35s %44s\r\n','Weekly','Daily');
    fprintf(fileID,'%21s %22s %22s %22s\r\n', output.RIC{:}, output.RIC{:});
    fprintf(fileID,'%5s %15.2f %22.2f\r\n', 'Mean:', output.stat.mu);
    fprintf(fileID,'%5s %15.2f %22.2f\r\n', 'Std:', output.stat.sigma);
    fprintf(fileID,'%5s %7s%2.2f %2.2f] %11s%2.2f %2.2f]\r\n', 'CI:','[',output.stat.CI(1,:),'[',output.stat.CI(2,:));
    fprintf(fileID,'%5s %15.2f %22.2f %22.2f %22.2f\r\n', 'Skew:', output.stat.skew);
    fprintf(fileID,'%5s %15.2f %22.2f %22.2f %22.2f\r\n', 'Kurt:', output.stat.kurt);
    fprintf(fileID,'%5s [%2.2f %2.2f %2.2f %2.2f] [%2.2f %2.2f %2.2f %2.2f] [%2.2f %2.2f %2.2f %2.2f] [%2.2f %2.2f %2.2f %2.2f]\r\n\r\n', 'Perc:',output.stat.perc(1,:),output.stat.perc(2,:),output.stat.perc(3,:),output.stat.perc(4,:));
    
    fprintf(fileID,'Problem 2\r\n');
    fprintf(fileID,'%21s %22s\r\n', output.RIC{:});
    fprintf(fileID,'%8s %12.2f %22.2f\r\n', 'EWMA:', output.EWMA.obj);
    fprintf(fileID,'%7s: %12.2f %22.2f\r\n', 'lambda', output.EWMA.param);
    fprintf(fileID,'%8s %12.2f %22.2f\r\n', 'GARCH1:', output.GARCH.obj);
    fprintf(fileID,'%1s%1s%1s%1s]: %1s%2.2f %2.2f %2.2f] %7s%2.2f %2.2f %2.2f]\r\n', '[','s,','a,','b','[',output.GARCH.param(1:3),'[',output.GARCH.param(4:6));
    fprintf(fileID,'%8s %12.2f %22.2f\r\n', 'GARCH2:', output.GARCH.objVT);
    fprintf(fileID,'%1s%1s%1s%1s]: %1s%2.2f %2.2f %2.2f] %7s%2.2f %2.2f %2.2f]\r\n\r\n', '[','s,','a,','b','[',output.GARCH.paramVT(1:3),'[',output.GARCH.paramVT(4:6));
    
    fprintf(fileID,'Problem 3\r\n');
    fprintf(fileID,'%9s %11s %22s\r\n', 'A-corr:', output.RIC{:});
    fprintf(fileID,'%9s %11.2f %22.2f\r\n', 'Lag 1:', output.stat.acorr(1,:));
    fprintf(fileID,'%9s %11.2f %22.2f\r\n', 'Lag 2:', output.stat.acorr(2,:));
    fprintf(fileID,'%9s %11.2f %22.2f\r\n', 'Lag 3:', output.stat.acorr(3,:));
    fprintf(fileID,'%9s %11.2f %22.2f\r\n', 'Lag 4:', output.stat.acorr(4,:));
    fprintf(fileID,'%9s %11.2f %22.2f\r\n\r\n', 'Lag 5:', output.stat.acorr(5,:));
    
    fprintf(fileID,'%21s \r\n', 'Linear');
    fprintf(fileID,'%9s %11.2f \r\n\r\n', 'corr:', output.stat.corr);

    fprintf(fileID,'%21s %11s %10s %10s %10s\r\n', 'Normal', 'Stud.-t', 'Gumbel', 'Clayton', 'Frank');
    fprintf(fileID,'%9s %11.2f %11.2f %10.2f %10.2f %10.2f\r\n', 'log-L:', output.copulaLogL);
    fclose(fileID);
else
    fprintf('%s\n', 'PROBLEM 1')
    fprintf('%35s %44s\n','Weekly','Daily');
    fprintf('%21s %22s %22s %22s\n', output.RIC{:}, output.RIC{:});
    fprintf('%5s %15.2f %22.2f\n', 'Mean:', output.stat.mu);
    fprintf('%5s %15.2f %22.2f\n', 'Std:', output.stat.sigma);
    fprintf('%5s %7s%2.2f %2.2f] %11s%2.2f %2.2f]\n', 'CI:','[',output.stat.CI(1,:),'[',output.stat.CI(2,:));
    fprintf('%5s %15.2f %22.2f %22.2f %22.2f\n', 'Skew:', output.stat.skew);
    fprintf('%5s %15.2f %22.2f %22.2f %22.2f\n', 'Kurt:', output.stat.kurt);
    fprintf('%5s [%2.2f %2.2f %2.2f %2.2f] [%2.2f %2.2f %2.2f %2.2f] [%2.2f %2.2f %2.2f %2.2f] [%2.2f %2.2f %2.2f %2.2f]\n\n', 'Perc:',output.stat.perc(1,:),output.stat.perc(2,:),output.stat.perc(3,:),output.stat.perc(4,:));
    
    fprintf('Problem 2\n');
    fprintf('%21s %22s\n', output.RIC{:});
    fprintf('%8s %12.2f %22.2f\n', 'EWMA:', output.EWMA.obj);
    fprintf('%7s: %12.2f %22.2f\n', num2str(char(955)), output.EWMA.param);
    fprintf('%8s %12.2f %22.2f\n', 'GARCH1:', output.GARCH.obj);
    fprintf('%1s%1s %1s %1s]: %1s%2.2f %2.2f %2.2f] %7s%2.2f %2.2f %2.2f]\n', '[', num2str(char(963)), num2str(char(945)), num2str(char(946)),'[',output.GARCH.param(1:3),'[',output.GARCH.param(4:6));
    fprintf('%8s %12.2f %22.2f\n', 'GARCH2:', output.GARCH.objVT);
    fprintf('%1s%1s %1s %1s]: %1s%2.2f %2.2f %2.2f] %7s%2.2f %2.2f %2.2f]\n\n', '[', num2str(char(963)), num2str(char(945)), num2str(char(946)),'[',output.GARCH.paramVT(1:3),'[',output.GARCH.paramVT(4:6));
    
    fprintf('Problem 3\n');
    fprintf('%9s %11s %22s\n', 'A-corr:', output.RIC{:});
    fprintf('%9s %11.2f %22.2f\n', 'Lag 1:', output.stat.acorr(1,:));
    fprintf('%9s %11.2f %22.2f\n', 'Lag 2:', output.stat.acorr(2,:));
    fprintf('%9s %11.2f %22.2f\n', 'Lag 3:', output.stat.acorr(3,:));
    fprintf('%9s %11.2f %22.2f\n', 'Lag 4:', output.stat.acorr(4,:));
    fprintf('%9s %11.2f %22.2f\n\n', 'Lag 5:', output.stat.acorr(5,:));
    
    fprintf('%21s \n', 'Linear');
    fprintf('%9s %11.2f \n\n', 'corr:', output.stat.corr);

    fprintf('%21s %11s %10s %10s %10s\n', 'Normal', 'Stud.-t', 'Gumbel', 'Clayton', 'Frank');
    fprintf('%9s %11.2f %11.2f %10.2f %10.2f %10.2f\n', 'log-L:', output.copulaLogL);   
end
