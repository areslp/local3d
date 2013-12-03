function [X,X0,mapping,gamma] = gentoydata()
[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup, denoise, dim, di] = get_parameters();
% k=100; % sample number of each subspace
affine=false; % generate linear subspace
% affine=true; % generate affine subspace
S=subspace_num; % number of subspace


[X,X0,mapping]=gen_basedata(dim,S,di,k,affine);
temp=X-X0;
temp2=sum(temp.^2,1);
gamma=nnz(temp2)/length(temp);
