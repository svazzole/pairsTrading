symbols = {'AAPL', 'MSFT', 'XOM', 'JNJ', 'MCD', 'WFC', 'GE', 'JPM', 'PG', 'PFE'};
%symbols = sp100;

database = 'WIKI/';
start_date = '2008-01-01';
freq = 'daily';

str = strcat(database, symbols(1));
data = Quandl.get(str, 'collapse', freq , 'start_date', start_date, 'type', 'fints');
adj_close = fts2mat(data.Adj_Close);
data = fints(data.dates, adj_close, symbols(1));

for i=2:length(symbols);
    
    str = strcat(database, symbols(i));
    %tmpData = Quandl.get(str, 'collapse', freq , 'start_date', start_date, 'type', 'fints');
    tmpData = Quandl.get(str, 'start_date', start_date, 'type', 'fints');
    adj_close = fts2mat(tmpData.Adj_Close);
    tmpFT = fints(tmpData.dates, adj_close, symbols(i));
    data = merge(data, tmpFT);
    %if size(data,1) == size(tmpData, 1)
        %data = horzcat(data, tmpData.Adj_Close);
    %    data = merge(data, tmpData.Adj_Close);
    %end;
    
end;

data = fillts(data);

%%

stocks = data; %fints(data.dates, d, 'd', 'Stocks');

%%

prices = fts2mat(stocks);
rets = prices(2:end,:)./prices(1:end-1,:) - 1;

%%

N = size(stocks,2);

%% BACKTEST

width = 250;

last = size(rets, 1);
forecasts = zeros(last-width, N);
weights = zeros(last-width, N);
dailyMat = zeros(N, N, last-width);

for i=width+1:last
    
    data = rets(i-width:i-1,:);
    [f,w,M] = genVARforecasts(data);
    forecasts(i-width, :) = f;
    weights(i-width, :) = w;
    dailyMat(:, :, i-width) = M;
    
end;

%%

bps = 0/10000;
[payoff, ~] = addCommissions(rets(width+1:end,:), forecasts, bps, 'method', 'standard');
cpayoff = cumprod(1+payoff);

