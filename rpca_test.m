close all;
clear;

[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup, denoise] = get_parameters();

[X,mapping]=gen_toydata();

size(X)


disp('algorithm begin');

% draw_points3d(X');

num=size(X,2);

% RPCA 
[A, E] = inexact_alm_rpca(X, lambda);
r=rank(A)
Z=A; % 3xn

draw_points3dm(X',A','rpca');

% size(A)
% draw_points3d(A');

% PCA
% Z1=X;
% r1=2;
% [U1 S1 V1]=svd(Z1);
% V1=V1(:,1:r1);
% U1=U1(:,1:r1);
% S1=diag(S1);
% S1=S1(1:r1);
% X1=U1*diag(S1)*V1';

% draw_points3dm(X',X1','pca');


n=length(X);

[U S V]=svd(Z);
V=V(:,1:r);
Z=V*V';
% size(Z)

[U S V]=svd(Z);
% svp=length(find(diag(S)>0));
% fprintf(1,'svp is %d\n',svp);
% U=U(:,1:svp);
% S=S(:,1:svp);
% V=V(:,1:svp);
U_tiled=U*(S.^0.5);
% row normalized, due to the paper "AUTOMATIC DETERMINATION OF THE NUMBER OF CLUSTERS USING SPECTRAL ALGORITHMS", row normalized is not necessary
for i=1:n
    U_tiled(i,:)=U_tiled(i,:)/norm(U_tiled(i,:));
end
W=U_tiled*U_tiled';
W=W.^2;
Z=W;



% vis
h=figure('Visible', 'off');
imagesc(Z);
colormap(gray);
axis equal;
saveas(h,'Z.png');

num_clusters=subspace_num;
cluster_labels=zeros(num,1);

% Spectral Clustering by Shi
% find(isnan(W))
% find(isinf(W))

[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(Z,num_clusters);


for i=1:num
    label=find(NcutDiscrete(i,:));
    cluster_labels(i)=label;
end


% draw_points3d_labels(X',cluster_labels);

ff=fopen('result.txt','w');
for i=1:num
    fprintf(ff,'%d %d\n',mapping(i,2)-1,cluster_labels(i)); % remember to -1 for meshlab
end
fclose(ff);
