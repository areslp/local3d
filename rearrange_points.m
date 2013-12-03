% function [] = rearrange_points(ours, baoli)
function [] = rearrange_points(varargin)
% baoli's method will change the point order, so we rearrange it

for v=1:numel(varargin)
    disp(['varargin{' num2str(v) '} class is ' class(varargin{v})]);
end

n=numel(varargin);

if n==0
    ours='data/standard_normal.xyzn';
    baoli='baoli/out.xyzn';
end

if n==1
    baoli=varargin{1};
end

if n==2
    ours=varargin{1};
    baoli=varargin{2};
end


data0=load(ours);
data1=load(baoli);

num=size(data0,1);

points0=data0(:,[1 2 3]);
vertex0=points0';
index0=kdtree_build(vertex0');

points1=data1(:,[1 2 3]);
vertex1=points1';
index1=kdtree_build(vertex1');

ff=fopen(baoli,'w');
for i=1:num
    query=vertex0(:,i);
    result = kdtree_k_nearest_neighbors(index1,query',1);
    idx=result(1);
    fprintf(ff,'%f %f %f %f %f %f\n',data1(idx,1),data1(idx,2),data1(idx,3),data1(idx,4),data1(idx,5),data1(idx,6));
end
fclose(ff);


kdtree_delete(index0);
kdtree_delete(index1);
end

