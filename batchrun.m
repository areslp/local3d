close all;
clear;


[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num] = get_parameters();

% prepare the standard_normal model
% TODO: noise model do not need this
% cd data
% system(['save_normal.exe ' fname '.obj']);
% cd ..


[point_with_normal]=load(['data/' fname '.xyzn']);
covv=load('cov.txt'); % 根据曲率值过滤不靠谱的字典法矢
size(covv)
covl=(covv<1e-2); % 靠谱点的label
% size(covl)
% pause;

points=point_with_normal(:,[1 2 3]);
vertex=points';
% save vertex.mat vertex;
% [vertex, face]=read_mesh(['data/' fname '.off']); % vertex : 3xn

index=kdtree_build(vertex');

% find(isnan(vertex))
% find(isinf(vertex))

% pause;

% build_params.algorithm='kdtree_single';
% build_params.trees=1;
% [index, parameters, speedup] = flann_build_index(vertex, build_params);

normals=point_with_normal(:,[4 5 6]);

num=length(vertex)

normals_new=zeros(size(normals));
vertex_new=zeros(size(vertex));

% mex call, parameters (vertex,normals,index,lambda), return normals_new
global_flag=zeros(num,1);

disp('algorithm begin');
tic;
skip=0;

% DEBUG
% temp=36683; % fandisk parallel planes
% temp=6747; % parallel1, close-by and parallel surfaces
temp=282; % parallel, close-by k=400
% temp=3; % vertex of icosahedron, k=500, subspace_num=5
% temp=71219; % noise_fandisk corner 
% temp=9038; % bunny ear
% temp=52255; % bunny body smooth
% temp=3555; % thinbox k=400;
temp=511;
% for i=temp:temp
for i=1:num
    if mod(i,10000)==0 
        fprintf(1,'processing %dth point\n',i);
        fprintf(1,'skip points %d\n',skip);
    end
    % fprintf(1,'processing %d points\n',i);

    % fprintf(1,'============================================================================\n');

    if global_flag(i)==1
        skip=skip+1;
        continue;
    end

    [X,mapping,idx]=genrealdata_batch(i,index,vertex,normals);

    % write neighbor index
    ff1=fopen('neighbors.txt','w');
    for tt=1:length(mapping)
        fprintf(ff1,'%d\n',mapping(tt,2)-1);
    end
    fclose(ff1);

    % draw_points3d(X');

    D=X;
    dell=zeros(length(D),1); % D中要删除点索引
    for ii=1:length(D)
        id=mapping(ii,2); % ii点在全局中的索引
        if(covl(id)==0) 
            dell(ii)=1;
        end
    end
    % size(D)
    % size(find(dell==1))
    % D(:,find(dell==1))=[]; % 从字典中删除不靠谱点的法矢

    % tic
    % [Z,E]=low_rank(X,lambda,1000); % 0.03, 0.04
    % [Z,E]=ladmp_lrr_fast(X,lambda,rho,DEBUG); % 0.01, 0.02
    [Z,E]=ladmp_lrr_fast_acc(X,lambda,rho,DEBUG); % 0.005,0.006,0.007
    % [Z,E]=ladmp_lrr_fast_acc_withD(X,D,lambda,rho,DEBUG); % 0.005,0.006,0.007 
    % [Z,E]=low_rank_acc(X,lambda,1000); % 0.005, 0.006, 0.007

    % [Z]=SSQP(X,1000); % 0.005,0.006,0.007
    % E=zeros(size(X));

    % TODO: non-negative
    % [Z,E]=nnlow_rank(X,lambda,1000); % 0.01, 0.02
    
    % toc
    % fprintf(1,'lowrank learning takes:%f\n',t);
    % Z=Z'*Z;

    % vis
    % h=figure('Visible', 'off');
    % imagesc(Z);
    % colormap(gray);
    % axis equal;
    % saveas(h,'Z.png');

    % tic;
    try
        [normals_new,global_flag,vertex_new]=cut(Z,E,vertex,vertex_new,normals,normals_new,mapping,global_flag,idx); % 谱分割聚类
    catch err
        getReport(err)
        fprintf(2,'error point idx is %d\n',i);
        draw_points3d(X');
        draw_points3d(D');
    end
    % t=toc;
    % fprintf(1,'clustering takes:%f\n',t);

    % tic;
    % normals_new(i,:)=lsqnormest2(vertex,idxs);
    % t=toc;
    % fprintf(1,'recompute normal takes:%f\n',t);

    % Xnew=D*Z; % 直接用XZ恢复正确法矢
    % normals_new(i,:)=Xnew(:,idx)'; 
end
t=toc;
fprintf(1,'process %d points takes %f\n',num,t);
fprintf(1,'point skip ratio is %f\n',skip/num);

% clear
% flann_free_index(index);
kdtree_delete(index);

% output

% TODO: Hack for noise model, we save the normals to original points
% [point_with_normal1]=load(['data/fandisk.xyzn']);
% points1=point_with_normal(:,[1 2 3]);
% vertex1=points1';
ff=fopen('out.xyzn','w');
for i=1:num
    fprintf(ff,'%f %f %f %f %f %f\n',vertex(1,i),vertex(2,i),vertex(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
    % fprintf(ff,'%f %f %f %f %f %f\n',vertex1(1,i),vertex1(2,i),vertex1(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
end
fclose(ff);


% ff=fopen('out_denoise.xyzn','w');
% for i=1:num
    % fprintf(ff,'%f %f %f %f %f %f\n',vertex_new(1,i),vertex_new(2,i),vertex_new(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
    % % fprintf(ff,'%f %f %f %f %f %f\n',vertex1(1,i),vertex1(2,i),vertex1(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
% end
% fclose(ff);

% system('normal_orientation.exe out.xyzn out_orientation.xyzn');

% system(['consistent_normal.exe out.xyzn --ori data/standard_normal.xyzn']);
% system(['consistent_normal.exe out_denoise.xyzn --ori data/standard_normal.xyzn']);
