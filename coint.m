%% COINTEGRATION

logPrices = log(prices);
rets = logPrices(2:end,:)./logPrices(1:end-1,:) - 1;

numAssets = size(prices,2);
numDays = size(prices,1);

cMat = zeros(numAssets,numAssets);
coint = []; %zeros((numAssets^2 - numAssets)/2, 5);
k = 1;

for i=1:numAssets
    for j=i+1:numAssets
        
        tmpPrices = [logPrices(:,i) logPrices(:,j)];
        [h,~,~,~,reg] = egcitest(tmpPrices);
        if (h==1)
            cMat(i,j) = reg.coeff(2);
            cMat(j,i) = cMat(i,j);
            coint(k,:) = [i j reg.coeff(1) reg.coeff(2) 0];
            k = k + 1;
        end;
    end;
end;

k = k - 1; % Number of cointegrated pairs

%%

spreads = zeros(numDays,k);
quot = zeros(numDays,k);

for i=1:k
    
    tmpPrices = [logPrices(:,coint(i,1)) logPrices(:,coint(i,2))];
    c0 = coint(i,3);
    beta = [1; -coint(i,4)];
    spreads(:,i) = zscore(tmpPrices*beta - c0);
    quot(:,i) = prices(:,coint(i,1))./prices(:,coint(i,2));
    %spreads(:,i) = tmpPrices*beta - c0;
end;

%% Compute positions 

p = positionPair(spreads);
pq = positionPair(quot);


%% Some diagnostic plots

i = 5;
lQuot = logPrices(:,coint(i,1))./logPrices(:,coint(i,2));
h1 = subplot(2,2,1); plot(1:numDays, logPrices(:,coint(i,1)), 1:numDays, logPrices(:,coint(i,2))); axis tight; grid on; hold on;
h3 = subplot(2,2,3); plot(spreads(:,i)); axis tight; grid on; hold on;
                     plot(2*ones(numDays,1), 'r');
                     plot(-2*ones(numDays,1), 'r');
                     plot(zeros(numDays,1), 'g');
h2 = subplot(2,2,2); plot(lQuot); axis tight; grid on; hold on;
                     plot(mean(lQuot)*ones(numDays,1), 'g');
                     plot((mean(lQuot) + 2 * std(lQuot))*ones(numDays,1), 'r');
                     plot((mean(lQuot) - 2 * std(lQuot))*ones(numDays,1), 'r');
h4 = subplot(2,2,4); stairs(p(:,i), 'linewidth', 1); axis tight; grid on;
linkaxes([h1, h2, h3, h4], 'x');

%% Profit and losses

pl = profitAndLosses(p, prices, coint);
plq = profitAndLosses(pq, prices, coint);
totRets = sum(pl,2);
cumRets = cumprod(totRets + 1);
totRetsQuot = sum(plq,2);
cumRetsQuot = cumprod(totRetsQuot + 1);

%% Some other diagnostic plots

i = 3;
numDays = size(prices,1);
lQuot = prices(:,coint(i,1))./prices(:,coint(i,2));
h1 = subplot(2,2,1); plot(1:numDays, prices(:,coint(i,1)), 1:numDays, prices(:,coint(i,2))); axis tight; grid on; hold on;
h3 = subplot(2,2,3); plot(spreads(:,i)); axis tight; grid on; hold on;
                     plot(2*ones(numDays,1), 'r');
                     plot(-2*ones(numDays,1), 'r');
                     plot(zeros(numDays,1), 'g');
h2 = subplot(2,2,2); plot((cumprod(1+pl(:,i)))); axis tight; grid on;
h4 = subplot(2,2,4); stairs(positions(:,i), 'linewidth', 1); axis tight; grid on;
linkaxes([h1, h2, h3, h4], 'x');






