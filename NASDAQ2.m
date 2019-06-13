%%  FINAL PROJECT - STOCK MARKET ANALYSIS
%   EE 573 - Random Signal Analysis and Kalman Filtering
%   
%   *DESCRIPTION*
%       This program...
%
%   *EXAMPLES*
%       Give examples of the code here...
%     
%   *NOTES*
%       Currently data resolution is limited to only daily values. We are 
%       in the process of finding a new way to get a better resolution.
%           Current Ideas:
%               - Make many pull requests and record each value
%               - Find a website that records hourly data (see API_Keys.m)
%
%       Anything else the user needs to be aware of...
%
%   *Acknowledgements*
%       Huge thank you to Josiah Renfree for providing the source code for
%       importing historical stock data.
%==========================================================================
%% IMPORT STOCK MARKET DATA
%==========================================================================
%Start with a blank workspace
%==========================================================================
clear all; close all; clc; 


%==========================================================================
% Import Data - See hist_stock_data for details on how data is collected
%==========================================================================
start_date = '14112012'; % Enter a start date: ddmmyyy
end_date = '14112017'; % Enter an end date: ddmmyyy
stocks = hist_stock_data(start_date,end_date,{'AAPL','^IXIC'});

%==========================================================================
%Plot Stock(s)
%==========================================================================
figure; hold on;
for i = 1:length(stocks)
    if i==length(stocks)
        yyaxis right
        ylabel('Price NASDAQ ($)');
    else
        yyaxis left
        ylabel('Price AAPL ($)');
    end
    plot(datenum(stocks(i).Date),stocks(i).Close);
end

%==========================================================================
% Format Plot
%==========================================================================
datetick('x','mmmyy');
xlabel('Date'); 
 
title('Historical Stock Prices');
legend(stocks.Ticker,'Location','southeast');
hold off; grid minor;

%==========================================================================
% FFT - get an idea of the frequency content of the data
%==========================================================================
figure;
FFT = db(abs(fft(stocks.Close)));
plot(FFT); % Frequency content of signal
xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)'); title('FFT of AAPL'); legend(stocks.Ticker);
%==========================================================================
%% STATISTICS

%==========================================================================
% Parameters to be varied
%==========================================================================
Window = 11; % lagging moving average window size

%==========================================================================
% Lagging moving average
%==========================================================================
movingMean = movmean(stocks(1).Close,[Window-1 0]);
movingMeanNASDAQ=movmean(stocks(2).Close,[Window-1 0]);

%==========================================================================
% Residual and Residual Variance ?
%==========================================================================
zeroMean = stocks(1).Close - movingMean;
zeroMeanNASDAQ = stocks(2).Close - movingMeanNASDAQ;
varMovingZeroMean= movvar(zeroMean,Window);
varMovingZeroMeanNASDAQ= movvar(zeroMeanNASDAQ, Window);

%==========================================================================
% Checking Covariance
%==========================================================================
%zeroMean(1:1258)
%stocks(1).Volume(1:1258)
A=cov(zeroMean(1:1258),zeroMeanNASDAQ(1:1258));  %to get cov(X,Y)
B1=cov(zeroMean(2:1259),zeroMean(1:1258)); %to get cov(X,Z)
B2=cov(zeroMean(2:1259),zeroMeanNASDAQ(1:1258)); %to get cov(Y,Z)
B=[B1(1,2);B2(1,2)];
M=A^-1*B;
a=M(1)
b=M(2)

%==========================================================================
% Plot STATISTICS AAPL
%==========================================================================
figure;
subplot(3,1,1)
    plot(stocks(1).Close); hold on
    plot(movingMean)
    xlabel('Time (Days)');
    ylabel('Price ($)');
    title('Moving Mean');
    secondLegend = 'Moving Mean (' + string(Window) +  ' Day Average)';
    legend('Raw Data',secondLegend,'Location','southeast');
    hold off; grid minor;

subplot(3,1,2);
    plot(zeroMean); hold on
    xlabel('Time (Days)'); 
    ylabel('Amplitude'); 
    title('Residuals');
    legend('Residual Data','Location','northeast');
    hold off; grid minor;
    
subplot(3,1,3);
    plot(varMovingZeroMean,'LineWidth',3)
    xlabel('Time (Days)'); 
    ylabel('Amplitude'); 
    title('Variance');
    legend('Variance of Zero Mean Data','Location','northeast');
    hold off;
%==========================================================================
% Plot STATISTICS NASDAQ
%==========================================================================
figure;
subplot(3,1,1)
    plot(stocks(2).Close); hold on
    plot(movingMeanNASDAQ)
    xlabel('Time (Days)');
    ylabel('Price ($)');
    title('Moving Mean NASDAQ');
    secondLegend = 'Moving Mean NASDAQ (' + string(Window) +  ' Day Average)';
    legend('Raw Data',secondLegend,'Location','southeast');
    hold off; grid minor;

subplot(3,1,2);
    plot(zeroMeanNASDAQ); hold on
    xlabel('Time (Days)'); 
    ylabel('Amplitude'); 
    title('Residuals');
    legend('Residual Data','Location','northeast');
    hold off; grid minor;
    
subplot(3,1,3);
    plot(varMovingZeroMeanNASDAQ,'LineWidth',3)
    xlabel('Time (Days)'); 
    ylabel('Amplitude'); 
    title('Variance');
    legend('Variance of Zero Mean Data NASDAQ','Location','northeast');
    hold off;    
%==========================================================================
%% Creating a Model

%==========================================================================
% Parameters to be varied
%==========================================================================
Threshold = .1; % Buy/Sell condition [$] 
AROrderZeroMean = 1; % order of AR process for zeroMean

%==========================================================================
% Allocate space
%========================================================================== 
zeroMeanEst = zeros(length(movingMean)+1,1);
movingMeanEst = zeros(length(movingMean)+1,1);
gainLossCumulative = zeros(length(movingMean),1);
buySellPrediction = zeros(length(movingMean)+1,1);
buySellRaw = zeros(length(movingMean)+1,1);
AZeroMean = zeros(AROrderZeroMean,AROrderZeroMean);
BZeroMean = zeros(AROrderZeroMean,1);
CZeroMean = zeros(AROrderZeroMean,1);
dayZeroMean = zeros(1,AROrderZeroMean);
buyAndHold = zeros(length(movingMean),1);
numTradesEst = 0;


%==========================================================================
% Define first few points to avoid zeros that distort the plots
%==========================================================================
lastKnown = AROrderZeroMean; % Date before prediction starts
zeroMeanEst(1:lastKnown+1) = zeroMean(1:lastKnown+1);
movingMeanEst(1:lastKnown+1) = movingMean(1:lastKnown+1);

%==========================================================================
% movingMeanEst CALCULATION - 2-Day Linear Prediction Process
%==========================================================================
for i=lastKnown+1:length(movingMean)
    movingMeanEst(i+1) = 2*movingMean(i)-movingMean(i-1);
end

%==========================================================================
% Residual CALCULATION from movingMeanEst
%==========================================================================
zeroMean = stocks(1).Close - movingMeanEst(1:length(movingMeanEst)-1);
% *** Run this zeroMean through SC573 to verify which AR order to use *** %

%==========================================================================
% zeroMeanEst CALCULATION - AR(n) Process
%==========================================================================
for i=lastKnown+1:length(zeroMean)

% Calculate Cross Correlation Vector - Recalculated for each point along the signal
    xcorrZeroMean = xcorr(zeroMean(1:i),AROrderZeroMean,'coeff'); 

% Create Covariance Matrix (A) - Recalculated for each point along the signal
    for j=1:AROrderZeroMean
        for k = 1:AROrderZeroMean
            AZeroMean(j,k) = xcorrZeroMean(AROrderZeroMean+j-k+1);
        end
    end
    
% Create R Matrix (B) - Recalculated for each point along the signal
    numCoeffsZeroMean = floor(length(xcorrZeroMean)/2);
    coeffZeroMean = xcorrZeroMean(ceil(length(xcorrZeroMean)/2):length(xcorrZeroMean));
    BZeroMean = coeffZeroMean(2:length(coeffZeroMean));
    
% Solve for AR coefficients Ax=B - Recalculated for each point along the signal
    ARCoeffZeroMean = AZeroMean\BZeroMean;
    
% Create a vector of days for each coefficient of the AR process - Recalculated for each point along the signal
    for kk = 1:AROrderZeroMean 
        dayZeroMean(kk) = zeroMean(i-kk+1);
    end
    
% Create the zeroMeanEst for the specified AR process
    zeroMeanEst(i+1) = dayZeroMean*ARCoeffZeroMean;
    
end

%==========================================================================
% Complete Model
%==========================================================================
randomness = normrnd(0,var(zeroMean-zeroMeanEst(1:length(zeroMean))),1260,1);

stockEstimateEst = movingMeanEst + zeroMeanEst;

%==========================================================================
%% EVALUATING MODEL

%==========================================================================
% Calculate Estimator Difference and Raw Difference of Tomorrow - Today
%==========================================================================
diffEst = stockEstimateEst(1:length(stocks(1).Close))-stocks(1).Close;
diffRaw = stocks(1).Close(2:length(stocks(1).Close))-stocks(1).Close(1:length(stocks(1).Close)-1);

%==========================================================================
% Create gainLossCumulative which tracks the profit of the model
%==========================================================================
for i=lastKnown:length(stocks(1).Close)
%Number of trades based on prediction
    if  diffEst(i) > Threshold
        buySellPrediction(i) = 1; % Buy 
        if i ~= length(zeroMean) % Exclude last day because we don't have the last day +1 true value
            gainLossCumulative(i+1) = gainLossCumulative(i)+(stocks(1).Close(i+1)-stocks(1).Close(i)); % buy today, sell tomorrow
        end
        numTradesEst = numTradesEst + 1;
    elseif diffEst(i) < -Threshold
        buySellPrediction(i) = -1; % Sell
        gainLossCumulative(i+1) = gainLossCumulative(i);
    else
        buySellPrediction(i) = 0; % Nothing
        gainLossCumulative(i+1) = gainLossCumulative(i);
    end
end

%==========================================================================
% Create the Theoretical Max/Min Profit/share if guesses were 100% accurate
% or inaccurate. Also Calculate the Profit/share if someone bought on day
% 1, didn't change anything, and sold on day k.
%==========================================================================
theoreticalMax = zeros(length(stocks(1).Close),1);
theoreticalMin = zeros(length(stocks(1).Close),1);

for k = 2:length(stocks(1).Close)
    
    % Theoretical Maximum and Minimum
    if stocks(1).Close(k-1) < stocks(1).Close(k)
        theoreticalMax(k) = theoreticalMax(k-1) + (stocks(1).Close(k) - stocks(1).Close(k-1));
        theoreticalMin(k) = theoreticalMin(k-1);
    elseif stocks(1).Close(k-1) > stocks(1).Close(k)
        theoreticalMax(k) = theoreticalMax(k-1);
        theoreticalMin(k) = theoreticalMin(k-1) + (stocks(1).Close(k) - stocks(1).Close(k-1));
    else
        theoreticalMax(k) = theoreticalMax(k-1);
        theoreticalMin(k) = theoreticalMin(k-1);
    end
    
    % Buying and Holding
    buyAndHold(k) = stocks(1).Close(k)-stocks(1).Close(1);
end

%==========================================================================
% Probably the most important plot
%==========================================================================
% figure;
%     hold on
%         plot(theoreticalMax,'r'); 
%         plot(buyAndHold);
%         plot(gainLossCumulative);
%         plot(theoreticalMin,'b');
%     hold off;
%     xlabel('Days');
%     ylabel('Efficiency [$/share]');
%     title('Efficiency');
%     legend('Theoretical Maximum','Buy and Hold','EE 573 Prediction','Theoretical Minimum','Location','NorthWest');
%     grid minor;

%==========================================================================
% Usefull Indicators of How the Model is Performing
%==========================================================================
zeroMeanErr = zeroMean - zeroMeanEst(1:length(zeroMean)); % Error in true and predicted residual - Is this white?
movingMeanErr = movingMean-movingMeanEst(1:length(movingMean)); % Error in true and predicted movingMean

% Usefull Values
errZeroMeanEst = xcorr(zeroMeanErr,AROrderZeroMean,'coeff') % Correlation of error between zeroMean and zeroMeanEst --- should be close to 0.
totalProfit = gainLossCumulative(length(gainLossCumulative))
%profitPerTrade = totalProfit/numTradesEst % Calculates profitPerTrade
buyAndHoldProfit = buyAndHold(length(buyAndHold))

%==========================================================================
% Display the Important Data We've Calculated
%==========================================================================
% Which process dominates the error  
% figure;
%     plot(zeroMeanErr); hold on
%     plot(movingMeanErr);
%     xlabel('Days');
%     ylabel('Errors');
%     title('Dominating source of Error');
%     legend('zeroMeanErr','movingMeanErr');
% 
% % True stock and Estimated stock
% figure; 
%     plot(stocks(1).Close); hold on
%     plot(stockEstimateEst);
%     xlabel('Days');
%     ylabel('Price ($)');
%     title('Raw data with Model');
%     legend(stocks.Ticker,'Prediction');

% % Difference between tomorrow estimate and today
% figure;
%     stem(diffEst);
%     xlabel('Days');
%     ylabel('Difference');
%     title('Tomorrows Prediction - Today');
%     legend('Difference');

% % Plot zeroMean and zeroMeanEst
% figure; 
%     plot(zeroMean); hold on;
%     plot(zeroMeanEst)
%     xlabel('Days');
%     ylabel('Price ($)');
%     title('Zero Mean Data Structure with Modeled');
%     legend(stocks.Ticker,'Modeled');
    
% % Plot detection of a buy/sell signal
% figure; 
%     stem(buySellPrediction)
%     xlabel('Days');
%     ylabel('Buy(1), Sell(-1), or Nothing(0)');
%     title('Buy or Sell');

%% END
%==========================================================================
%==========================================================================