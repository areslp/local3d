function [query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup, denoise, dim, di] = get_parameters()
% fandisk.off
fname='test';
% fname='fandisk';
% fname='noise_fandisk2';
% fname='noise_fandisk';
% fname='noise_fandisk_7w';
% fname='armadillo';
% fname='noise_armadillo';
% fname='bunny';
% fname='noise_bunny';
% fname='double-torus1';
% fname='noise_double-torus1';
% fname='cone';
% fname='icosahedron';
% fname='noise_icosahedron';
% fname='consist_sharp_sphere';
% fname='parallel1';
% fname='parallel';
% fname='thinbox';
% fname='mechpart';

query_id=9038; % bunny ear
query_id=161913; % armadillo tiny geometric feature
% query_id=102961; % armadillo tiny geometric feature

k=5; % local neighborhood number
% lambda=1/sqrt(k);
% lambda=0.1414; % 50
lambda=0; % 100
% lambda=1;
alpha=0.2;
rho=1.9;
DEBUG=0;
% tau=6e-1; % with 0.2 noiseend
tau=8e-1; % with no noise and PCA normal
% tau=1e-1; % with no noise and correct normal
% tau=8e-2; % paper motion segmentation
subspace_num=2;
speedup=false;
% speedup=true;
denoise=false;


%%%% parameters for toy data
%%%% as in the paper

% dim=100; % dim of ambient space
% di=4;

% lines in R^3

dim=3; % dim of ambient space
di=1;

% planes in R^3

% dim=3; % dim of ambient space
% di=2; % dim of subspace

% lines in R^2

% dim=2;
% di=1;

% points in R^3

% dim=3;
% di=0;
