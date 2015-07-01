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