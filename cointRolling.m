tic;

w = 252;                                % window length
nAssets = size(prices,2);               % number of assets
totDays = size(prices,1);
nDays = totDays - w;                    % number of days out of sample
maxCoint = (nAssets^2 - nAssets)/2;     % max number of coint pairs

spreads = zeros(nDays,maxCoint);
logPrices = log(prices);
coint = zeros(maxCoint,4);
k = 0;

for d=w+1:totDays
    
    disp([num2str((d-w)/nDays * 100) '%']);
    cMat = zeros(nAssets,nAssets);
    logPricesRolling = logPrices(d-w:d,:);
    
    for i=1:nAssets
        for j=i+1:nAssets
        
            tmpPrices = [logPricesRolling(:,i) logPricesRolling(:,j)];
            [h,~,~,~,reg] = egcitest(tmpPrices, 'alpha', 0.99);
            if (h==1)
                cMat(i,j) = reg.coeff(2);
                cMat(j,i) = cMat(i,j);
                nr = find((coint(:,1) == i) & (coint(:,2) == j));
                if isempty(nr)
                    k = k + 1;
                    coint(k,:) = [i j reg.coeff(1) reg.coeff(2)];
                    c0 = coint(k,3);
                    beta = [1; -coint(k,4)];
                    s = zscore(tmpPrices*beta - c0);
                    spreads(d,k) = s(end);
                else
                    coint(nr,:) = [i j reg.coeff(1) reg.coeff(2)];
                    c0 = coint(nr,3);
                    beta = [1; -coint(nr,4)];
                    s = zscore(tmpPrices*beta - c0);
                    spreads(d,nr) = s(end);
                end;
                
            end;
        end;
    end;
    
end;

%%

positions = positionPair(spreads);
pl = profitAndLosses(positions, prices, coint);
totRets = sum(pl,2);
cumRets = cumprod(totRets + 1);

%%
cumProfits = [];

for i=1:maxCoint
    
    cumProfits(:,i) = cumprod(1+pl(:,i));
    
end;

toc;