function [X,X0,mapping,gamma] = gentoydata()
[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup, denoise, dim, di] = get_parameters();
% k=100; % sample number of each subspace
affine=false; % generate linear subspace
% affine=true; % generate affine subspace
S=subspace_num; % number of subspace


<<<<<<< HEAD
dim=100; % dim of ambient space
di=4;

% lines in R^3

% dim=3; % dim of ambient space
% di=1;

% planes in R^3

% dim=3; % dim of ambient space
% di=2; % dim of subspace

% lines in R^2

% dim=2;
% di=1;

[X,mapping]=gen_basedata(dim,S,di,k,affine);
=======
[X,X0,mapping]=gen_basedata(dim,S,di,k,affine);
temp=X-X0;
temp2=sum(temp.^2,1);
gamma=nnz(temp2)/length(temp);
>>>>>>> 0b3ab03acaea6bd8981178e706d0776e28065fea
