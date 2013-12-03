close all;
clear;

[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num] = get_parameters();

[point_with_normal]=load(['data/' fname '.xyzn']);
points=point_with_normal(:,[1 2 3]);
vertex=points';
% [vertex, face]=read_mesh(['data/' fname '.off']); % vertex : 3xn
index=kdtree_build(vertex');
% build_params.algorithm='kdtree';
% build_params.trees=1;
% [index, parameters, speedup] = flann_build_index(vertex, build_params);

normals=point_with_normal(:,[4 5 6]);

num=length(vertex);

disp('algorithm begin');
normals_new=main_loop(vertex,normals,index,lambda);

% clear
% flann_free_index(index);
kdtree_delete(index);

% output

% TODO: Hack for noise model, we save the normals to original points
[point_with_normal1]=load(['data/fandisk.xyzn']);
points1=point_with_normal(:,[1 2 3]);
vertex1=points1';
ff=fopen('out.xyzn','w');
for i=1:num
    % fprintf(ff,'%f %f %f %f %f %f\n',vertex(1,i),vertex(2,i),vertex(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
    fprintf(ff,'%f %f %f %f %f %f\n',vertex1(1,i),vertex1(2,i),vertex1(3,i),normals_new(i,1),normals_new(i,2),normals_new(i,3));
end
fclose(ff);
