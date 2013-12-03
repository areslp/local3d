function n = lsqnormest2(X,idxs)

points=X(:,idxs);
x = points;
% draw_points3d(x');
k=size(points,2);
p_bar = 1/k * sum(x,2);

P = (x - repmat(p_bar,1,k)) * transpose(x - repmat(p_bar,1,k)); %spd matrix P
%P = 2*cov(x);

% add epsilon to diag to prevent singular
% epsilon=1e-3;
% P=P+diag(repmat(epsilon,3,1));

% issymmetric(P)
[V,D] = eig(P);
% diag(D)
% V

[~, idx] = min(diag(D)); % choses the smallest eigenvalue

n = V(:,idx);   % returns the corresponding eigenvector    
% size(n) % 3x1
