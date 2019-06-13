clear, clc, close all

%% Import data
start_date = '14112012'; % Enter a start date: ddmmyyy
end_date = '14112017'; % Enter an end date: ddmmyyy
stocks = hist_stock_data(start_date,end_date,{'AAPL'});

estimate= stocks.Close(1);
updated_estimate= estimate;
H= 1;
v= 0;
R= 0;
Phi= 1;
var_noise= 1;
Q= var_noise;
Pminus= var_noise;
I= 1;
noise= normrnd(0,var_noise,length(stocks.Close)-1,1);
mse= 0;

for k= 2:length(stocks.Close)
    z(k)= stocks.Close(k);
    x(k)= H*z(k);
%     var_noise= 1;
%     Pminus= var_noise(k);
%     noise= normrnd(0,var_noise,1,1);
    
    estimate(k)= Phi*updated_estimate + noise(k-1);
    error= x(k) - estimate(k-1);
    
    Kgain= Pminus*H*(H*Pminus*H+R);
    P= (I-Kgain*H)*Pminus;
    updated_estimate= estimate(k) + Kgain*(z(k)-H*estimate(k));
    
    Pminus= Phi*P*Phi+Q;
    mse= mse+(x(k)-estimate(k))^2;
end
disp(mse)

figure(1)
plot(stocks.Close)
hold on
plot(estimate)