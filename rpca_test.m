close all;
clear;

[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup, denoise] = get_parameters();

[X,mapping]=gen_toydata();

rank(X)

disp('algorithm begin');

draw_points3d(X');

num=size(X,2);

% lowrank
[A, E] = inexact_alm_rpca(X, lambda);

rank(A)

