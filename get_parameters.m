function [query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup, denoise] = get_parameters()
% fandisk.off
% fname='fandisk';
% fname='noise_fandisk2';
% fname='noise_fandisk';
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
fname='mechpart';

query_id=9038; % bunny ear
query_id=161913; % armadillo tiny geometric feature
% query_id=102961; % armadillo tiny geometric feature

k=100; % local neighborhood number
% lambda=1/sqrt(k);
lambda=1;
alpha=0.1;
rho=1.9;
DEBUG=0;
% tau=6e-1; % with 0.2 noiseend
tau=8e-1; % with no noise and PCA normal
% tau=1e-1; % with no noise and correct normal
% tau=8e-2; % paper motion segmentation
subspace_num=0;
speedup=false;
% speedup=true;
denoise=true;
