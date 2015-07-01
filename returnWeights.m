function r = returnWeights(prices, positions)
    
    rPrices = 1./prices;
    
    rpr = rPrices.*abs(positions);
    nAssets = size(rpr,2);
    nDays = size(rpr,1);
    
    for i=1:nAssets
        
        for j=2:nDays
        
            if (rpr(j,i)~=0 && rpr(j-1,i)==0)
                
                s = j;
                
                while j <= nDays && rpr(j,i)~=0
                    
                    j = j+1;
                    
                end;
                
                rpr(s:j-1,i) = rpr(s,i) .* ones(j-s,1);
                
            end;
            
        end;
        
    end;
    
    r = rpr;

end