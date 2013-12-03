function [normals_new,global_flag] = cut(Z,X)

[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup] = get_parameters();

[m n]=size(Z);
assert(m==n);
num=m;

Z=0.5*(abs(Z)+abs(Z'));
% Z=normc(Z);
Z=mnormalize_col(Z);

% h=figure('Visible', 'off');
% imagesc(Z);
% colormap(gray);
% axis equal;
% saveas(h,'Z_ori.png');
figure;
set(gcf, 'Color', 'w');
% Z(find(Z<0.2))=0; % LRSR
Z(find(Z<0.05))=0; % LRR
% Z(find(Z<0.1))=0; % SR
spy(Z);
axis off;
% set(gca,'xtick',[],'ytick',[]);
% set(gca,'xlabel',[ ],'ylabel',[]);


% clustering

% Z=0.5*(abs(Z)+abs(Z'));

% skinny SVD and compute W=([U~U~^T])^2 
% U~=U^*(\Sigma^*)^{1/2} with row normalized
% skinny svd

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

% sigma=0;
% kmean_clusters=num_clusters;
% spectual clustering normalized
% Choose spectral clustering algorithm
% - 1 for unnormalized Clustering
% - 2 for normalized   Clustering (Shi and Malik)
% - 3 for normalized   Clustering (Jordan and Weiss)

% ClusteringType=3;
% L = SpectralClustering(W, num_clusters, kmean_clusters, ClusteringType);

% label=find(L(id,:));
% idxs=find(L(:,label));

% for i = 1:num
    % label=find(L(i,:)==1);
    % assert(length(label)==1);
    % cluster_labels(i)=label(1); % the ith sample's label is label(1)
% end

% if write_result
    % draw_points3d_labels(X',cluster_labels);
% end

% now I have z1 and cluster_labels
% z1

% write result

% ff=fopen('result.txt','w');
% for i=1:num
    % fprintf(ff,'%d %d\n',mapping(i,2)-1,cluster_labels(i)); % remember to -1 for meshlab
% end
% fclose(ff);
