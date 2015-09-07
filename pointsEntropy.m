function H3 = pointsEntropy(pA, w, h)    
    m=1;
    n = 4;
    X = sort(pA(:,1));
    Y = sort(pA(:,2));
    c = 1;
    
    % Normaize X and Y values
    X = X * 100 / w;
    Y = Y * 100 / h;
    
    H = 0;H2 = 0;
    for j=1:n-m
        H = H + log(n/m*(max(X(j+m)-X(j),c))); %%%incomplete
    end 
              
    for j=1:n-m
        H2 = H2 + log(n/m*(max(Y(j+m)-Y(j),c))); %%%incomplete
    end
    H = H/n;   
    H2 = H2/n;
    H3 = min(H,H2);