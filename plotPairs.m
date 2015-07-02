function plotPairs(results, n)
    
    try
        prices = results.prices;
        spreads = results.spreads;
        cointRel = results.cointRel;
        positions = results.positions;
        pl = results.pl;
    catch 
        error('Are you kidding me?');
    end;
        
    nDays = size(prices,1);
    
    h1 = subplot(2,2,1); 
        plot(1:nDays, prices(:,cointRel(n,1)), 1:nDays, prices(:,cointRel(n,2))); axis tight; grid on;
    h2 = subplot(2,2,2);
        plot(1:nDays, spreads(:,n)); axis tight; grid on; hold on;
        m = mean(spreads(:,n));
        sd = std(spreads(:,n));
        plot(1:nDays, m + 2*sd*ones(nDays,1), 'r', 1:nDays, m - 2*sd*ones(nDays,1), 'r');
        plot(1:nDays, m*ones(nDays,1), 'g'); hold off;
    h3 = subplot(2,2,3);
        plot(1:nDays, pl(:,n)); axis tight; grid on;
    h4 = subplot(2,2,4);
        stairs(1:nDays, positions(:,n)); axis tight; grid on;
        ylim([-1.5 1.5]);
    
    linkaxes([h1 h2 h3 h4], 'x');
    
end