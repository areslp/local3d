function [X,mapping] = genrealdata()
[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k] = get_parameters();

[point_with_normal]=load(['data/' fname '.xyzn']);
points=point_with_normal(:,[1 2 3]);
vertex=points';
% [vertex, face]=read_mesh(['data/' fname '.off']); % vertex : 3xn
tree=kdtree_build(vertex');

query=vertex(:,query_id); % the 0 th point
result = kdtree_k_nearest_neighbors(tree,query',k);
% clear
kdtree_delete(tree);

% delete query point itself, only need to delete query point itself when we use p-q as samples

% idx=find(result==query_id);
% result(idx)=[];

num=length(result);

% mapping between index

mapping=zeros(num,2);
for i=1:num
    mapping(i,:)=[i result(i)];
end

% training data using points

% X=zeros(3,num);
% for i=1:num
    % X(:,i)=vertex(:,result(i))-query;
    % % normalize samples
    % % X(:,i)=X(:,i)/norm(X(:,i));
% end
% draw_direction3d(X);

% using coordinate directly

% X=vertex(:,result);

% using Shape Context feature instead

% X=load('sc.txt')';
% X=X(:,result);

% using normal

normals=point_with_normal(:,[4 5 6]);
X=normals(result,:)';

% make the normal is not normalized, coefficients should > 0

% coefficients=randn(1,num);
% for i=1:num
    % X(:,i)=X(:,i)*coefficients(i);
% end
