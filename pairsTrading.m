function results = pairsTrading(prices, varargin)

    p = inputParser;
    
    defaultMethod = 'standard';
    expectedMethods = {'standard','log','quotient'};
    defaultConf = 0.95;
    defaultWindowWidth = 252;
    defaultBps = 10;
    
    addRequired(p, 'prices', @isnumeric);
    addOptional(p, 'method', defaultMethod, @(x) any(validatestring(x, expectedMethods)));
    addOptional(p, 'confLev', defaultConf, @isnumeric)
    addOptional(p, 'window', defaultWindowWidth, @isnumeric);
    addOptional(p, 'bps', defaultBps, @isnumeric);
    
    parse(p, prices, varargin{:});
    
    %%%%%%%%%%%%%%%%%%
    % Set parameters %
    %%%%%%%%%%%%%%%%%%
    
    confLevel = p.Results.confLev;
    
    if(confLevel >= 1 || confLevel <= 0)
        error('Confidence level must be in (0,1)'); 
    end;
    
    bps = p.Results.bps;
    
    if (bps < 0)
        error('Trading costs must be >= 0');
    end;
    
    w = p.Results.window;                       % window length
    numAssets = size(prices,2);                 % number of assets
    n = (numAssets*numAssets-numAssets)/2;      % max number of cointegration's relations
    totDays = size(prices,1);                   % number of days

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
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Start pair trading on every couple. TODO: manage different MATLAB versions %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    nDays = totDays - w;                        % number of days out of sample
    spreads = zeros(nDays,n);                   % array for spreads
    cbA = spreads;
    ubA = spreads;
    lbA = spreads;
    cointRel = zeros(n,4);                      % table for cointegration
    
    POOL = parpool('local');

    parfor_progress(n);

    parfor i=1:n         
                
        for d=w+1:totDays
            
            tmpPrices = [prices(d-w:d, couples(i,1)) prices(d-w:d, couples(i,2))];
            
            c = cointParam([tmpPrices(:,1) tmpPrices(:,2)], couples(i,:), confLevel);
            
            if ~isempty(c)
                
                cointRel(i,:) = c;
                c0 = c(3);
                beta = [1; -c(4)];
                s = zscore(tmpPrices*beta - c0);
                
                spreads(d,i) = s(end);
                cbA(d,i) = mean(s);
                ubA(d,i) = mean(s) + 2*std(s);
                lbA(d,i) = mean(s) - 2*std(s);
            
            else
                cbA(d,i) = 0;
                ubA(d,i) = 2;
                lbA(d,i) = -2;
            end;
            
        end;
        
        parfor_progress;
    
    end;
    
    parfor_progress(0);
    
    delete(POOL);
    
    %matlabpool open;

%     switch p.Results.method
%         case 'standard'
%             [spreads, cointRel, cbA, ubA, lbA] = standardSpread(prices, couples, w, totDays, confLevel);
%         case 'log'
%             lPrices = log(prices);
%             [spreads, cointRel, cbA, ubA, lbA] = standardSpread(lPrices, couples, n, w, totDays, confLevel);
%         otherwise
%             [spreads, cointRel, cbA, ubA, lbA] = quotientSpread(prices, couples, n, w, totDays, confLevel);
%     end;
   
    %matlabpool close;

    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now clean cointRel, spreads, prices and bands %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    t = cointRel(:,1) ~= 0;
    cointRel = cointRel(t,:);
    spreads = spreads(w+1:end,t);
    cbA = cbA(w+1:end,t);
    ubA = ubA(w+1:end,t);
    lbA = lbA(w+1:end,t);
    prices = prices(w+1:end,:);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Then obtain positions and p&l %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    positions = positionPair(spreads, cbA, ubA, lbA);
    pl = profitAndLosses(positions, prices, cointRel, bps);
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Return results in struct %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    results.pl = pl;
    results.prices = prices;
    results.cointRel = cointRel;
    results.spreads = spreads;
    results.cbA = cbA;
    results.ubA = ubA;
    results.lbA = lbA;
    results.positions = positions;
    
    results.cumulativeRets = cumprod(sum(pl,2) + 1);
    results.sumRets = cumsum(sum(pl,2));
    results.totPL = sum(pl,2);
    
    warning('on');
    
end

function [spreads, cointRel, cbA, ubA, lbA] = standardSpread(prices, couples, w, totDays, confLevel)
    
    n = size(couples, 1);
    nDays = totDays - w;                        % number of days out of sample
    spreads = zeros(nDays,n);                   % array for spreads
    cbA = spreads;
    ubA = spreads;
    lbA = spreads;
    cointRel = zeros(n,4);                      % table for cointegration
    
    POOL = parpool('local');

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
                cbA(d,i) = mean(s);
                ubA(d,i) = mean(s) + 2*std(s);
                lbA(d,i) = mean(s) - 2*std(s);
            
            else
                cbA(d,i) = 0;
                ubA(d,i) = 2;
                lbA(d,i) = -2;
            end;
            
        end;
        
        parfor_progress;
    
    end;
    
    parfor_progress(0);
    
    delete(POOL);

end

function [spreads, cointRel, cbA, ubA, lbA] = quotientSpread(prices, couples, n, w, totDays, confLevel)

    nDays = totDays - w;                        % number of days out of sample
    spreads = zeros(nDays,n);                   % array for spreads
    cbA = spreads;
    ubA = spreads;
    lbA = spreads;
    cointRel = zeros(n,4);                      % table for cointegration

    POOL = parpool('local');

    parfor_progress(n);

    parfor i=1:n         
        
        coup = couples(i,:);
        
        for d=w+1:totDays
            
            tmpPrices = [prices(d-w:d, coup(1)) prices(d-w:d, coup(2))];
            
            c = cointParam([tmpPrices(:,1) tmpPrices(:,2)], coup, confLevel);
            
            if ~isempty(c)
                
                cointRel(i,:) = c;
                s = tmpPrices(:,1)./tmpPrices(:,2);
                
                spreads(d,i) = s(end);
                cbA(d,i) = mean(s);
                ubA(d,i) = mean(s) + 2*std(s);
                lbA(d,i) = mean(s) - 2*std(s);
            
            else
                cbA(d,i) = 0;
                ubA(d,i) = 2;
                lbA(d,i) = -2;
            end;
            
        end;
        
        parfor_progress;
    
    end;
    
    parfor_progress(0);
    
    delete(POOL);
    
end

function c = cointParam(p, couple, level)

    [h,~,~,~,reg] = egcitest(p, 'alpha', level);
    
    if (h==1)
        c = [couple reg.coeff(1) reg.coeff(2)];
    else
        c = [];
    end;
    
end

function p = positionPair(spreads, cbA, ubA, lbA)
    
    ncol = size(spreads,2);
    nDays = size(spreads,1);
    
    p = zeros(nDays, ncol);
    
    for i=1:ncol
        
        t = spreads(:,i);
        cb = cbA(:,i);
        ub = ubA(:,i);
        lb = lbA(:,i);

        for j=3:nDays
            
            if (p(j-1,i) == 0) % No position
                
                if (t(j-2) > ub(j-2) && t(j-1) < ub(j-1))
                    p(j,i) = -1;
                elseif (t(j-2) < lb(j-2) && t(j-1) > lb(j-2))
                    p(j,i) = 1;
                end;
                
            elseif (p(j-1,i) == 1) % Buy greater sell the lower
                
                if (t(j-2) < cb(j-2) && t(j-1) >= cb(j-1))
                    p(j,i) = 0;
                else
                    p(j,i) = 1;
                end;
            
            else % p(j-1) == -1 Buy the lower sell the greater
            
                if (t(j-2) > cb(j-2) && t(j-1) <= cb(j-1))
                    p(j,i) = 0;
                else
                    p(j,i) = -1;
                end;
                
            end;
            
        end;
        
    end;

end

function pl = profitAndLosses(positions, prices, coint, bps)
    
    k = size(coint,1);
    n = size(prices,1);
    pl = zeros(n,k);
    rets = zeros(n,2);
    
    bps = bps/10000;
    
    for i=1:k
        
        if coint(i,1) ~= 0
            
            pos = positions(:,i);
            b = bps .* commissions(pos);
            
            pr = [prices(:,coint(i,1)) prices(:,coint(i,2))];
            
            w = returnWeights(pr, [pos pos]);
            
            rets(2:end,:) = pr(2:end,:)./pr(1:end-1,:) - 1;
            pl(:,i) = ((rets(:,1)-b).*w(:,1) - (rets(:,2)+b).*w(:,2)).*pos;
            
        end;
        
    end;
    

end

function b = commissions(pos)

    % N.B.: the trading costs must change sign 
    % depending on the position (-1 or 1)!
    
    d = size(pos,1);
    b = zeros(d,1);
    
    for i=2:d
        if pos(i-1) == 0 && pos(i) == 1
            b(i) = 1;
        elseif pos(i-1) == 0 && pos(i) == -1
            b(i) = -1;
        elseif pos(i-1) == 1 && pos(i) == 0
            b(i-1) = b(i-1) + 1;
        elseif pos(i-1) == -1 && pos(i) == 0
            b(i-1) = b(i-1) - 1;
        end;
    end;

end


