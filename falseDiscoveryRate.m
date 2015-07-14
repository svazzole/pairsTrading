function [fdr, logic, ix] = falseDiscoveryRate(pvalues, confidence)

    [pOrd, ix] = sort(pvalues);
    l = size(pOrd,1);
    b = zeros(l,1);
    
    confidence = 1 - confidence;
    
    for i=1:l
        
        b(i) = (i/(l*sum(1./(1:l))))*confidence;
        
    end;
    
    t = (pOrd <= b);
    k = find(t==1, 1, 'last');
    
    if k~=0
    
        t(1:k) = ones(k,1);
    
    end;
    
    fdr = ix.*t;
    fdr(fdr==0) = [];
    fdr = sort(fdr);
    
    logic = t(ix);
    
end