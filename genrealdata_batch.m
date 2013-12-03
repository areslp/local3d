function [X,mapping,idx] = genrealdata_batch(id,tree,vertex,normals)
[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k] = get_parameters();
query=vertex(:,id); % the 0 th point

result = kdtree_k_nearest_neighbors(tree,query',k);

% parameters.checks=32;
% [result, dists] = flann_search(tree, query, k, parameters);

num=length(result);
idx=find(result==id);
% mapping between index
mapping=zeros(num,2);
for i=1:num
    mapping(i,:)=[i result(i)];
end
X=normals(result,:)';

% V=vertex(:,result);
% draw_points3d(V');

% TODO: TESTING
% nn=size(X,2);
% co=randn(1,nn);
% X=X.*repmat(co,3,1);
