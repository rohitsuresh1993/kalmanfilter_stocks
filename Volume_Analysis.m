%% Volume Analysis

%Trying to determine the correlation between volume and the change in AAPL
%stock

clear all; close all; clc; 

start_date = '14112012'; % Enter a start date: ddmmyyy
end_date = '14112017'; % Enter an end date: ddmmyyy
stocks = hist_stock_data(start_date,end_date,{'AAPL'});
stocks(1).Volume=stocks(1).Volume/(10^7);
figure; hold on;
for i = 1:length(stocks)
    plot(datenum(stocks(i).Date),stocks(i).Close);
end
datetick('x','mmmyy');
xlabel('Date'); 
ylabel('Price ($)'); 
title('Historical Stock Prices');
legend(stocks.Ticker,'Location','southeast');
hold off; grid minor;

figure;
plot(stocks(1).Volume)

diffRaw = stocks(1).Close(2:length(stocks(1).Close))-stocks(1).Close(1:length(stocks(1).Close)-1);


Window = 10;
movingMean = movmean(stocks(1).Volume,[Window-1 0]);
zeroMean = stocks(1).Volume - movingMean;
varMovingZeroMean= movvar(zeroMean,Window);
stdMovingZeroMean=sqrt(varMovingZeroMean);

figure;
plot(stdMovingZeroMean)
xlabel('Time (Days)'); 
ylabel('Amplitude'); 
title('STD');
legend('STD of Zero Mean Data','Location','northeast');
figure;
findpeaks(stocks(1).Volume,'MinPeakProminence',8)

[peaks,day]=findpeaks(stocks(1).Volume,'MinPeakProminence',8);

%for i=1:length(day)
diffRawChopped=diffRaw(day(1:length(day)));

figure;
yyaxis left
plot(peaks,'--')
ylabel('Volume (in 100 millions)')
hold on
yyaxis right
plot(diffRawChopped)
ylabel('Difference in Price')
hold off
grid minor

figure;
yyaxis left
plot(stocks(1).Close)
ylabel('Price of AAPL')
yyaxis right
plot(stocks(1).Volume)
ylabel('Volume (in 100 millions)')
grid minor

% j=zeros(length(stocks(1).Volume),1);
% %Capture first few high points
% for i=1:Window
%     if stocks(1).Volume(i)> (2*10^8)
%         j(i)=i
%     end
% end

% for i=(Window+1):1259
%     
%    VolChopped=stocks(1).Volume(i-10:i-1);
%    meanWindow=mean(VolChopped);
%    stdWindow=std(VolChopped);
%    if stocks(1).Volume(i)> (3*stdWindow)
%     
% end
% figure;
% subplot(3,1,1)
%     plot(stocks(1).Volume); hold on
%     plot(movingMean)
%     xlabel('Time (Days)');
%     ylabel('Price ($)');
%     title('Moving Mean');
%     secondLegend = 'Moving Mean (' + string(Window) +  ' Day Average)';
%     legend('Raw Data',secondLegend,'Location','southeast');
%     hold off; grid minor;
% 
% subplot(3,1,2);
%     plot(zeroMean); hold on
%     xlabel('Time (Days)'); 
%     ylabel('Amplitude'); 
%     title('Residuals');
%     legend('Residual Data','Location','northeast');
%     hold off; grid minor;
%     
% subplot(3,1,3);
%     plot(varMovingZeroMean)
%     xlabel('Time (Days)'); 
%     ylabel('Amplitude'); 
%     title('Variance');
%     legend('Variance of Zero Mean Data','Location','northeast');
%     hold off;