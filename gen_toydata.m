function [X,mapping] = gentoydata()
[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num] = get_parameters();
k=100; % sample number of each subspace
affine=false; % generate linear subspace
% affine=true; % generate affine subspace
S=subspace_num; % number of subspace

% as in the paper

% dim=100; % dim of ambient space
% di=5;

% lines in R^3

% dim=3; % dim of ambient space
% di=1;

% planes in R^3

dim=3; % dim of ambient space
di=2; % dim of subspace

% lines in R^2

% dim=2;
% di=1;

[X,mapping]=gen_basedata(dim,S,di,k,affine);
