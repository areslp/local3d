close all;
clear;


[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num] = get_parameters();

% prepare the standard_normal model
% TODO: noise model do not need this
% cd data
% system(['save_normal.exe ' fname '.obj']);
% cd ..


[point_with_normal]=load(['data/' fname '.xyzn']);
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

num=length(vertex);

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
temp=1;
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
    % if i==temp
        % fprintf(1,'%d\n',temp);
    % end
    % tic;
    % lowrank
    [X,mapping,idx]=genrealdata_batch(i,index,vertex,normals);

    % t=toc;
    % fprintf(1,'generate X takes:%f\n',t);

    % write neighbor index
    % ff1=fopen('neighbors.txt','w');
    % for tt=1:length(mapping)
        % fprintf(ff1,'%d\n',mapping(tt,2)-1);
    % end
    % fclose(ff1);



    % draw_points3d(X');

    tic
    % [Z,E]=low_rank(X,lambda,1000); % 0.03, 0.04
    % [Z,E]=ladmp_lrr_fast(X,lambda,rho,DEBUG); % 0.01, 0.02
    [Z,E]=ladmp_lrr_fast_acc(X,lambda,rho,DEBUG); % 0.005,0.006,0.007
    % [Z,E]=low_rank_acc(X,lambda,1000); % 0.005, 0.006, 0.007

    % TODO: non-negative
    % [Z,E]=nnlow_rank(X,lambda,1000); % 0.01, 0.02
    
    t=toc;
    fprintf(1,'lowrank learning takes:%f\n',t);

    % if i==64
        % % vis
        % h=figure('Visible', 'off');
        % imagesc(Z);
        % colormap(gray);
        % axis equal;
        % saveas(h,'Z.png');
    % end

    % DEBUG C++

    % tic;
    % try
        [normals_new,global_flag,vertex_new]=cut(Z,E,vertex,vertex_new,normals,normals_new,mapping,global_flag,idx);
        % return;
    % catch
        % fprintf(2,'error point idx is %d\n',i);
        % draw_points3d(X');
    % end
    % idxs=cut(X,Z,idx,mapping,true);
    % t=toc;
    % fprintf(1,'clustering takes:%f\n',t);

    % tic;
    % normals_new(i,:)=lsqnormest2(vertex,idxs);
    % t=toc;
    % fprintf(1,'recompute normal takes:%f\n',t);

    % Xnew=X*Z;
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


ff=fopen('out_denoise.xyzn','w');
for i=1:num
    fprintf(ff,'%f %f %f %f %f %f\n',vertex_new(1,i),vertex_new(2,i),vertex_new(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
    % fprintf(ff,'%f %f %f %f %f %f\n',vertex1(1,i),vertex1(2,i),vertex1(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
end
fclose(ff);

% system('normal_orientation.exe out.xyzn out_orientation.xyzn');

system(['consistent_normal.exe out.xyzn --ori data/standard_normal.xyzn']);
system(['consistent_normal.exe out_denoise.xyzn --ori data/standard_normal.xyzn']);
