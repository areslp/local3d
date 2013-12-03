function [normals_new,global_flag,vertex_new] = cut(Z,E,vertex,vertex_new,normals,normals_new,mapping,global_flag,id)

[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup, denoise] = get_parameters();

[m n]=size(Z);
assert(m==n);
num=m;

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

% save point's neighbor 
% for i=1:n
    % coefficients=Z(:,i);
    % idx=find(coefficients>0);
    % ff=fopen(['debug/' num2str(i) '.txt'],'w');
        % fprintf(ff,'%d %d\n',mapping(i,2)-1,0);
        % for j=1:length(idx)
            % id=idx(j);
            % fprintf(ff,'%d %d\n',mapping(id,2)-1,coefficients(id));
        % end
    % fclose(ff);
% end

if subspace_num==0


    D=diag(sum(Z,2));
    L=eye(n)-D^-0.5*Z*D^-0.5;
    [U S V]=svd(L);
    S=diag(S);

    % if write_result
        % S'
    % end

    % predict the cluster number by eigenvalues(lowrank paper)
    assert(length(S)==n);
    soft_thresholding=0;
    for i=1:n
        if S(i)>=tau
            soft_thresholding=soft_thresholding+1;
        else
            soft_thresholding=soft_thresholding+log2(1+(S(i)/tau)^2);
        end
    end
    subspace_num=n-round(soft_thresholding);
    % if write_result
        % fprintf(1,'subspace number by eigenvalues(lowrank paper) is %d\n',subspace_num);
    % end

    % predict the cluster number by eigenvalues gap
    % [Y,I]=max(abs(diff(S)));
    % subspace_num=n-I;

    % if write_result
        % fprintf(1,'subspace number by eigenvalues gap is %d\n',subspace_num);
    % end
end

% if write_result
    % vis
    % h=figure('Visible', 'off');
    % imagesc(Z);
    % colormap(gray);
    % axis equal;
    % saveas(h,'Z.png');
% end


% cut
% for i=1:size(Z,2)
    % if norm(Z(:,i))==0
        % fprintf(1,'Warning: Similar matrix has zero columns!\n');
        % continue;
    % end
    % Z(:,i)=Z(:,i)/norm(Z(:,i));
% end

num_clusters=subspace_num;
cluster_labels=zeros(num,1);

% Spectral Clustering by Shi
% find(isnan(W))
% find(isinf(W))

[NcutDiscrete,NcutEigenvectors,NcutEigenvalues] = ncutW(Z,num_clusters);


% outliers detection method1 TODO:
% Esum=sum(E.^2,1);
% threshold=1e-1;
% outliers=find(Esum>threshold);
% fprintf(1,'size of outliers: %d\n',length(outliers));

% TODO: dealing with outliers
% NcutDiscrete(outliers,:)=0;

% method2
ang=10/180*pi;

% iter though class
% fprintf(1,'num_clusters is %d\n',num_clusters);

% TODO: approximate
if speedup
    for i=1:num_clusters
        idxs=find(NcutDiscrete(:,i));
        
        numpcluster=length(idxs);
        % fprintf(1,'numpcluster is %d\n',numpcluster);
        
        % if numpcluster>3*k/4.0
        if numpcluster>=k/2.0
        % if numpcluster>=k
            % map idxs to global points index
            idxs=mapping(idxs,2);
            n=lsqnormest2(vertex,idxs);%3x1

            % is n.dot(normal) small? if small, ignore it; if large, leave it
            n=n/norm(n);
            costheta = dot(repmat(n',numpcluster,1),normals(idxs,:),2);
            theta = acos(costheta);

            % assert(length(theta)==length(idxs));
            
            flatidx=find(theta<ang);
            flatnumpcluster=length(flatidx); 

            normals_new(idxs(flatidx),:)=repmat(n',flatnumpcluster,1); 
            global_flag(idxs(flatidx))=1;

            % normals_new(idxs,:)=repmat(n',numpcluster,1); 
            % global_flag(idxs)=1;
            
            % if length(find(idxs==135871))~=0
                % fprintf(1,'@@@@@@@@@@@@@@@@@@@@@@@@@@id: %d\n',mapping(id,2));
            % end
            % fprintf(1,'ignore %d points\n',numpcluster);
        else
            continue;
        end
    end
end

% caculate the current point's normal
label=find(NcutDiscrete(id,:));
idxs=find(NcutDiscrete(:,label));
% map idxs to global points index
idxs=mapping(idxs,2);

% ff=fopen('temp.txt','w');
% for i=1:length(idxs)
    % fprintf(ff,'%d %d\n',idxs(i)-1,1); % remember to -1 for meshlab
% end
% fclose(ff);


% ff=fopen('points.txt','w');
% points=vertex(:,idxs);
% for i=1:length(idxs)
    % fprintf(ff,'%f %f %f\n',points(1,i),points(2,i),points(3,i));
% end
% fclose(ff);

% vertex(:,mapping(id,2))
% n=lsqnormest2(vertex,idxs);%3x1
% n
n=fitNormal(vertex(:,idxs)',false);
% [coeff]=princomp(vertex(:,idxs)');
% n=coeff(:,3);

% nn=length(idxs);
% normals=repmat(n',nn,1);
% draw_points_and_normals(vertex(:,idxs)',normals);
cur_id=mapping(id,2);
normals_new(mapping(id,2),:)=n; 

if denoise
    [Point, Normal] = lsqPlane(vertex(:,idxs)');
    [PB]=orthcomp(Normal'); % plane basis
    point=vertex(:,cur_id);
    vertex_new(:,cur_id)=projPointOnPlane(point', [Point PB(:,1)' PB(:,2)']);
end


% normals_new(mapping(id,2),:)

% ff=fopen('cur_point.txt','w');
% fprintf(ff,'%d %d\n',mapping(id,2)-1,1);
% fclose(ff);

global_flag(mapping(id,2))=1;


for i=1:num
    label=find(NcutDiscrete(i,:));
    cluster_labels(i)=label;
end

% save cr.mat cluster_labels;

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
