function results = pairsTrading(prices, varargin)

    p = inputParser;
    
    defaultMethod = 'standard';
    expectedMethods = {'standard','log','quotient'};
    defaultConf = 0.95;
    defaultWindowWidth = 252;
    
    addRequired(p, 'prices', @isnumeric);
    addOptional(p, 'method', defaultMethod, @(x) any(validatestring(x, expectedMethods)));
    addOptional(p, 'confLev', defaultConf, @isnumeric)
    addOptional(p, 'window', defaultWindowWidth, @isnumeric);
    
    parse(p, prices, varargin{:});
    
    %%%%%%%%%%%%%%%%%%
    % Set parameters %
    %%%%%%%%%%%%%%%%%%
    
    confLevel = p.Results.confLev;
    
    if(confLevel >= 1 || confLevel <= 0)
        error('Confidence level must be in (0,1)'); 
    end;
    
    w = p.Results.window;                       % window length
    numAssets = size(prices,2);                 % number of assets
    n = (numAssets*numAssets-numAssets)/2;      % max number of cointegration's relations
    totDays = size(prices,1);                   % number of days
    nDays = totDays - w;                        % number of days out of sample
    
    spreads = zeros(nDays, n);                  % array for spreads
    cointRel = zeros(n,4);                      % table for cointegration
    
    warning('off');                             % STFU!
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Generate all the couples. TODO: Is there a better way? %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    couples = zeros(n,2);
    k=1;
    for i=1:numAssets
        for j=i+1:numAssets
            couples(k,:) = [i j];
            k = k + 1;
        end;
    end;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start pair trading on every couple. TODO: PARALLEL %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %matlabpool open;
    POOL = parpool('local',8);
    
    parfor_progress(n);
    
    for i=1:n         
     
        for d=w+1:totDays
            
            tmpPrices = [prices(d-w:d, couples(i,1)) prices(d-w:d, couples(i,2))];
            
            c = cointParam([tmpPrices(:,1) tmpPrices(:,2)], couples(i,:), confLevel);
            
            if ~isempty(c)
                
                cointRel(i,:) = c;
                c0 = c(3);
                beta = [1; -c(4)];
                s = zscore(tmpPrices*beta - c0);
                spreads(d,i) = s(end);

            end;
            
        end;
        
        parfor_progress;
    
    end;
    
    parfor_progress(0);
    
    %matlabpool close;
    delete(POOL);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now clean cointRel, spreads and prices %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = cointRel(:,1) ~= 0;
    cointRel = cointRel(t,:);
    spreads = spreads(w+1:end, t);
    prices = prices(w+1:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Then obtain positions and p&l %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    positions = positionPair(spreads);
    pl = profitAndLosses(positions, prices, cointRel);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Return results in struct %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    results.pl = pl;
    results.prices = prices;
    results.cointRel = cointRel;
    results.spreads = spreads;
    results.positions = positions;
    
    results.cumulativeRets = cumprod(sum(pl,2) + 1);
    results.sumRets = sum(sum(pl,2));
    results.totPL = sum(pl,1);
    
    warning('on');
    
end

function w = optimalWindowWidth(t)

    %  Really needed ? 

end

function c = cointParam(p, couple, level)

    [h,~,~,~,reg] = egcitest(p, 'alpha', level);
    
    if (h==1)
        c = [couple reg.coeff(1) reg.coeff(2)];
    else
        c = [];
    end;
    
end

function p = positionPair(spreads)

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % TODO: ROLLING POSITIONS! %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    ncol = size(spreads,2);
    nDays = size(spreads,1);
    
    p = zeros(nDays, ncol);
    

    for i=1:ncol
        
        %y = spreads(1:end-1,i);
        t = spreads(1:end,i);
        
        cb = mean(spreads(:,i));
        ub = cb + 2 * std(spreads(:,i));
        lb = cb - 2 * std(spreads(:,i));
        
        for j=3:nDays
            
            if (p(j-1,i) == 0) % No position
                
                if (t(j-2) > ub && t(j-1) < ub)
                    p(j,i) = -1;
                elseif (t(j-2) < lb && t(j-1) > lb)
                    p(j,i) = 1;
                end;
                
            elseif (p(j-1,i) == 1) % Buy greater sell the lower
                
                if (t(j-2) < cb && t(j-1) > cb)
                    p(j,i) = 0;
                else
                    p(j,i) = 1;
                end;
            
            else % p(j-1) == -1 Buy the lower sell the greater
            
                if (t(j-2) > cb && t(j-1) < cb)
                    p(j,i) = 0;
                else
                    p(j,i) = -1;
                end;
                
            end;
            
        end;
        
    end;

end

function pl = profitAndLosses(positions, prices, coint)

    k = size(coint,1);
    n = size(prices,1);
    pl = zeros(n,k);
    rets = zeros(n,2);
    
    for i=1:k
        
        if coint(i,1) ~= 0
            pos = positions(:,i);
            pr = [prices(:,coint(i,1)) prices(:,coint(i,2))];
            
            w = returnWeights(pr, [pos pos]);
            
            rets(2:end,:) = pr(2:end,:)./pr(1:end-1,:) - 1;
            pl(:,i) = (rets(:,1).*w(:,1) - rets(:,2).*w(:,2)).*pos;
            
        end;
        
    end;
    

end
