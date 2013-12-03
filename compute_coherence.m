function [mu] = compute_coherence(X,r,gamma)
% compute_coherence as described in paper "Exact matrix completion via convex optimization"
[U,S,V]=svd(X);
[m,n]=size(X);
r=min(r,min(m,n))
U=U(:,1:r); % mxm
V=V(:,1:r); % nxn
Iu=eye(m); % n1
Iv=eye(n); % n2
% max U
maxu=0;
for i=1:m
    temp=norm(U'*Iu(:,i));
    temp=temp*temp;
    if maxu<temp
        maxu=temp;
    end
end
% max V
maxv=0;
for i=1:n
    temp=norm(V'*Iv(:,i));
    temp=temp*temp;
    if maxv<temp
        maxv=temp;
    end
end
% max UV'
temp=U*V';
maxuv=max(abs(temp(:)));
maxuv=maxuv*maxuv;
% mu
mumatrix=[m*maxu/r, n*maxv/r, maxuv*m*n/r]
% mu=max(mumatrix);
mu=mumatrix(2)*(1-gamma);
end

