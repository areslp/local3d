function [A] = construct_A(P) 
k=6;
n=size(P,1);
tree=kdtree_build(P);
A=sparse([],[],[],k*n,n,2*k*n);
for i=1:n
    query=P(i,:);
    result = kdtree_k_nearest_neighbors(tree,query,k+1);
    idx=find(result==i);
    result(idx)=[];
    assert(length(result)==k);
    for j=1:k
        A((i-1)*k+j,i)=1;
        A((i-1)*k+j,result(j))=-1;
    end
end
% clear
kdtree_delete(tree);
