function [pl, positions, spreads, cointRel] = pairsTrading(prices, varargin)

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
    numAssets = size(prices,2);
    n = (numAssets*numAssets-numAssets)/2;
    totDays = size(prices,1);
    nDays = totDays - w;                        % number of days out of sample
    spreads = zeros(nDays, n);
    cointRel = zeros(n,4);
    
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
    
    % k = 0;
    
    %matlabpool open;
    
    %parfor_progress(n);
    
    for i=1:n         
        disp(i);
        
        for d=w+1:totDays
            
            tmpPrices = [prices(d-w:d, couples(i,1)) prices(d-w:d, couples(i,2))];
            
            c = cointParam([tmpPrices(:,1) tmpPrices(:,2)], couples(i,:), confLevel);
            
            if ~isempty(c)
                
                % disp(['Cointegration found in the couple (' num2str(couples(i,1)) ', ' num2str(couples(i,2)) ')']);
                % do rolling pairs trading
                
                cointRel(i,:) = c;
                c0 = c(3);
                beta = [1; -c(4)];
                s = zscore(tmpPrices*beta - c0);
                spreads(d,i) = s(end);
                
%                nr = find((cointRel(:,1) == couples(i,1)) & (cointRel(:,2) == couples(i,2)));               
%                 if isempty(nr)
%                     k = k + 1;
%                     cointRel(k,:) = c;
%                     c0 = cointRel(k,3);
%                     beta = [1; -cointRel(k,4)];
%                     s = zscore(tmpPrices*beta - c0);
%                     spreads(d,k) = s(end);
%                 else
%                     cointRel(nr,:) = c;
%                     c0 = cointRel(nr,3);
%                     beta = [1; -cointRel(nr,4)];
%                     s = zscore(tmpPrices*beta - c0);
%                     spreads(d,nr) = s(end);
%                 end;

            end;
            
        end;
    end;
    
    %parfor_progress(0);
    
    %matlabpool close;
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Now clean spreads and cointRel %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Then obtain positions and p&l %
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    positions = positionPair(spreads);
    pl = profitAndLosses(positions, prices, cointRel);
    
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

    ncol = size(spreads,2);
    nDays = size(spreads,1);
    
    p = zeros(nDays, ncol);
    

    for i=1:ncol
        
        y = spreads(1:end-1,i);
        t = spreads(2:end,i);
        
        cb = mean(spreads(:,i));
        ub = cb + 2 * std(spreads(:,i));
        lb = cb - 2 * std(spreads(:,i));
        
        for j=2:nDays
            
            if (p(j-1,i) == 0) % No position
                
                if (y(j-1) > ub && t(j-1) < ub)
                    p(j,i) = -1;
                elseif (y(j-1) < lb && t(j-1) > lb)
                    p(j,i) = 1;
                end;
                
            elseif (p(j-1,i) == 1) % Buy greater sell the lower
                
                if (y(j-1) < cb && t(j-1) > cb)
                    p(j,i) = 0;
                else
                    p(j,i) = 1;
                end;
            
            else % p(j-1) == -1 Buy the lower sell the greater
            
                if (y(j-1) > cb && t(j-1) < cb)
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
