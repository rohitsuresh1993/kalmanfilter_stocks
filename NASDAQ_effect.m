%%  FINAL PROJECT - STOCK MARKET ANALYSIS
%   EE 573 - Random Signal Analysis and Kalman Filtering
%   
%   *DESCRIPTION*
%       This program...
%
%   *EXAMPLES*
%       Give examples of the code here...
%
%   *DATA STRUCTURE*
%       Talk about stuctures here...
%     
%   *OUTPUT FORMAT*
%       Talk about the outputs here...
%     
%   *DATA FEED*
%       Talk about where the data comes from...
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
%==========================================================================
%% IMPORT STOCK MARKET DATA
clear all; close all; clc; %Start with a blank workspace

    % Use function 'hist_stock_data' to collect stock data from Yahoo!. For
    %   more information see the specific function. Huge thank you to Josiah
    %   Renfree for providing the source code for importing historical stock data.

start_date = '06112016'; % Enter a start date: ddmmyyy
end_date = '06112017'; % Enter an end date: ddmmyyy

    % Example/Common tickers: GOOGLE-'GOOGL',APPLE-'AAPL',S&P500-'^GSPC'

stocks = hist_stock_data(start_date,end_date,{'AAPL','GOOG','^IXIC'});


    % stocks - Data Structure - struct
    %   Date:
    %   Open:
    %   High:
    %   Low:
    %   Close:
    %   AdjClose:
    %   Volume:
    %   Ticker:
    
% DIDN'T WANT TO KEEP CLOSING THESE FIGURES
%Plot Stock
figure(1); hold on;
for i = 1:length(stocks)
    if i==length(stocks)
        yyaxis right
    else
        yyaxis left
    end
    plot(datenum(stocks(i).Date),stocks(i).Close);
    
end

% Format Plot
datetick('x','mmmyy'); 
xlabel('Date'); 
ylabel('Price ($)'); 
title('Historical Stock Prices');
legend(stocks.Ticker,'Location','southeast');
hold off;
% 
% % FFT - get an idea of the frequency content of the data
% figure(2);
% FFT = db(abs(fft(stocks.Close)));
% plot(FFT); % Frequency content of signal
% xlabel('Frequency (rad/s)'); ylabel('Magnitude (dB)'); title('FFT of AAPL'); legend(stocks.Ticker);
%==========================================================================
%% STATISTICS

% All Raw Data
Window = 10; % Window Size
movingMean = movmean(stocks(1).Close,Window); % Window-Day moving average

%Moving Mean of NASDAQ
mmNASDAQ=movmean(stocks(2).Close,Window);

figure(3);
    subplot(3,1,1)
    plot(stocks(1).Close); hold on
    plot(movingMean)
    xlabel('Time (Days)');
    ylabel('Price ($)');
    title('Moving Mean');
    secondLegend = 'Moving Mean (' + string(Window) +  ' Day Average)';
    legend('Raw Data',secondLegend,'Location','southeast');
    hold off;

zeroMean = stocks(1).Close - movingMean; % Subtracting out the mean
zmNASDAQ=stocks(2).Close - mmNASDAQ; %REsidual of NASDAQ
varMovingRaw = movvar(stocks(1).Close,2); % Tracking Var
varMovingZeroMean= movvar(zeroMean,Window); %var of zero mean

subplot(3,1,2)
    plot(zeroMean); hold on
    xlabel('Time (Days)'); 
    ylabel('Amplitude'); 
    title('Stationary Mean');
    legend('Zero-Mean Data','Location','northeast');
    plot(zmNASDAQ);
    hold off;
    
subplot(3,1,3)
    plot(varMovingRaw); hold on
    plot(varMovingZeroMean,'LineWidth',3)
    xlabel('Time (Days)'); 
    ylabel('Amplitude'); 
    title('Variance');
    legend('Variance of Raw Data','Variance of Zero Mean Data','Location','northeast');
    hold off;
    
%==========================================================================
%% Creating a Model

% for ii = 1:length(varMovingZeroMean) % Creates Noise with changing variance
%     n(ii) = normrnd(0,sqrt(varMovingZeroMean(ii)),1,1); % Generates Random noise with variance function of zeroMean
% end

lag = 2;
coeff = xcorr(zeroMean,lag,'coeff') 

zmest = zeros(length(zeroMean),1);
zmest(1) = zeroMean(1);
for k = 3:length(zeroMean)
    zmest(k) = coeff(lag)*zeroMean(k-1);
end

figure(6);
    plot(zeroMean);hold on;
    plot(zmest)
    xlabel('Days');
    ylabel('Price ($)');
    title('Zero Mean Data Structure with Modeled');
    legend(stocks.Ticker,'Modeled');

% Error in true and predicted residual
reserr = zeroMean - zmest;
figure(12);
    plot(reserr)

coff2 = xcorr(reserr,lag,'coeff')
    
reserrest = zeros(length(zeroMean),1);
reserrest(1) = reserr(1);
for k = 2:length(reserr)
    reserrest(k) = coeff(lag)*reserr(k-1);
end

figure(13);
plot(reserr); hold on
plot(reserrest)

reserr2 = reserr - reserrest;
coff3 = xcorr(reserr2,lag,'coeff')
%Check for xcorr of errors until a very small value turns up

%% Predicting data

known=2;
StockEStimateest = zeros(length(stocks(1).Close+1),1);
zeroMeanest = zeros(length(stocks(1).Close+1)+1,1);
StockEstimateest(1:known)=stocks(1).Close(1:known);

A=cov(zeroMean(1:250),zmNASDAQ(1:250));  %to get cov(X,Y)
B1=cov(zeroMean(2:251),zeroMean(1:250)); %to get  XZ
B2=cov(zeroMean(2:251),zmNASDAQ(1:250)); %YZ
B=[B1(1,2);B2(1,2)]
M=A^-1*B;
a=M(1)
b=M(2)


Threshold = 0.2;
for i=known:length(zeroMean)
    movingMeanest(i+1)=a*movingMean(i)+b*mmNASDAQ(i); %predicting moving mean
    zeroMeanest(i+1)=a*zeroMean(i)+b*zmNASDAQ(i); %predicting residual
    StockEstimateest(i+1) = movingMeanest(i+1) + zeroMeanest(i+1);
    %Number of trades based on prediction
    if StockEstimateest(i+1)-stocks(1).Close(i) > Threshold
        buysell(i) = 1;
    elseif StockEstimateest(i+1)-stocks(1).Close(i) < -Threshold
        buysell(i) = -1;
    else
        buysell(i) = 0;
    end
    %Number of trades based on true data
    if stocks(1).Close(i)-stocks(1).Close(i-1) > Threshold
        buysell1(i) = 1;
    elseif stocks(1).Close(i)-stocks(1).Close(i-1) < -Threshold
        buysell1(i) = -1;
    else
        buysell1(i) = 0;
    end
end

figure(15)
stem(buysell)
figure(16)
stem(buysell1)
numtradesest=length(find(buysell~=0))
numtrades=length(find(buysell1~=0))

figure(14)
%True stock and Estimated stock
plot(StockEstimateest);hold on
plot(stocks(1).Close);

% A=cov(zeroMean(1:250),zmNASDAQ(1:250));
% B1=cov(zeroMean(2:251),zeroMean(1:250));
% B2=cov(zeroMean(2:251),zmNASDAQ(1:250));
% B=[B1(1,2);B2(1,2)]
% M=A^-1*B;
% a=M(1)
% b=M(2)



%==========================================================================
%% DATA FILTERING



%==========================================================================
%% DATA OPTIMIZATION



%==========================================================================
%% ANYTHING ELSE?



%% END
%==========================================================================
%==========================================================================
%==========================================================================
%==========================================================================