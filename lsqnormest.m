function n = lsqnormest(points)
    x = points;
    k=size(points,2);
    p_bar = 1/k * sum(x,2);
    
    P = (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %spd matrix P
    %P = 2*cov(x);
    
    [V,D] = eig(P);
    
    [~, idx] = min(diag(D)); % choses the smallest eigenvalue
    
    n = V(:,idx);   % returns the corresponding eigenvector    
