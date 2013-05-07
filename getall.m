close all;
clear;

[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num] = get_parameters();

% [X,mapping]=genrealdata();
[X,mapping]=gen_toydata();

disp('algorithm begin');

% draw_points3d(X');

num=size(X,2);

% lowrank
[Z,E]=ladmp_lrr_fast(X,lambda,rho,DEBUG);
% [Z,E]=low_rank(X,lambda,1000);
% [Z,E]=l1_low_rank(X,lambda,alpha,1000);

cut_simple(Z,X);
