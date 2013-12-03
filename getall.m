close all;
clear;

[query_id, fname, lambda, alpha, rho, DEBUG, tau, subspace_num, k, speedup, denoise, dim, di] = get_parameters();

<<<<<<< HEAD
% [X,mapping]=genrealdata();
% [X,mapping]=gen_toydata();

% save temp.mat X;
% load temp.mat; % noise data, ppt

% save temp1.mat X;
load temp1.mat;
=======
[X,X0,mapping,gamma_real]=gen_toydata();

X

disp('how well is the data? data properties');

r0=min(di*subspace_num,min(size(X,1),size(X,2)));
[UX SX VX]=svd(X,'econ');
[U0 S0 V0]=svd(X0,'econ');
SX(SX<10*eps)=0;

norm(X-X0)
temp=SX.^-1;
temp(temp==Inf)=0;
beta=1/(norm(X)*norm(temp*VX'*V0));
beta1=1/cond(X0)
E0=X-X0;
theta=min(mPrinAngles(E0,X0));
beta2=sin(theta)/(cond(X0)*(1+norm(E0)/norm(X0)))
fprintf(2,'beta is %f\n',beta);

mu=compute_coherence(X0,di*subspace_num,gamma_real);
fprintf(2,'coherence is %f\n',mu);

ta=324*(beta^2)
tb=49*(11+4*beta)^2*mu*r0
% tb=49*(11+4*beta)^2*1*r0 % mu最小也是1，这样算出来也不对。。。。。
gamma=ta/(ta+tb);
fprintf(2,'real gamma is %f, and gamma_* is %f\n',gamma_real,gamma);

% lambda=3/(7*norm(X)*sqrt(gamma_real*length(X)))

% save vars.mat X X0;
% load vars.mat;

% add 0 to X, make X\in R^100，结果并没有不同
% X=[X;zeros(97,2*k)];

cluster_labels=zeros(length(X),1);

% draw_points3d_labels(X0',cluster_labels);
draw_points3d_labels(X',cluster_labels);
>>>>>>> 0b3ab03acaea6bd8981178e706d0776e28065fea

disp('algorithm begin');


num=size(X,2);

% 先通过RPCA把$E_\perp$去掉
% [XNEW E iter] = inexact_alm_rpca_l21(X, 0.28);
% X=XNEW;

% lowrank
% [Z,E]=ladmp_lrr_fast(X,lambda,rho,DEBUG);
% [Z,E]=low_rank(X,lambda,1000);
% [Z,E]=l1_low_rank(X,lambda,alpha,1000);
[Z]=l1_graph(X,lambda);
<<<<<<< HEAD

=======
% [Z]=SSQP(X,lambda);
XX=X*Z;

% [Z0,E0]=low_rank(X0,lambda,1000);
% norm(E0)
% norm(E)
% norm(X-X0)

draw_points3d_labels(XX',cluster_labels);

% [u,s,v]=svd(X0);
% size(v)
% v=v(:,1:3);
% % v=v(:,1:20);
% Z00=v*v';
% h=figure('Visible', 'off');
% imagesc(Z00);
% colormap(gray);
% axis equal;
% saveas(h,'Z00.png');


% fprintf(1,'error norm is %f\n',norm(XX-X0));
% fprintf(1,'Z-Z0 norm is %f\n',norm(Z-Z0));
% fprintf(1,'Z-Z00 norm is %f\n',norm(Z-Z00));

% rank(Z(1:k,1:k))

% h=figure('Visible', 'off');
% imagesc(Z);
% colormap(gray);
% axis equal;
% saveas(h,'Z.png');

if subspace_num==2

    % cluster_labels_small=zeros(length(X)/2,1);
    % X1=X(:,1:k);
    % X2=X(:,k+1:end);
    % X01=X0(:,1:k);
    % X02=X0(:,k+1:end);
    % Z1=Z(1:k,1:k);
    % Z2=Z(k+1:end,k+1:end);

    % test with RPCA
    % [ZR1 ER1 iter] = inexact_alm_rpca_l21(X1, lambda);
    % [ZR2 ER2 iter] = inexact_alm_rpca_l21(X2, lambda);
    % [ZR1 ER1 iter] = inexact_alm_rpca(X1, lambda);
    % [ZR2 ER2 iter] = inexact_alm_rpca(X2, lambda);

    % LRR
    % draw_points3d_labels((X1*Z1)',cluster_labels_small);
    % draw_points3d_labels((X2*Z2)',cluster_labels_small);
    % draw_points3d_labels([(X1*Z1)';(X2*Z2)'],cluster_labels);

    % RPCA
    % draw_points3d_labels(ZR1',cluster_labels_small);
    % draw_points3d_labels(ZR2',cluster_labels_small);
    % draw_points3d_labels([ZR1';ZR2'],cluster_labels);

end




if subspace_num==1

    % [Z E iter] = inexact_alm_rpca_l21(X, lambda);
    % [Z E iter] = inexact_alm_rpca(X, lambda);
    % XXR=Z;
    % draw_points3d_labels(XXR',cluster_labels);
end
Z=full(Z)
return;
>>>>>>> 0b3ab03acaea6bd8981178e706d0776e28065fea

cut_simple(Z,X,X0);
% cut_simple(Z00,X,X0);
