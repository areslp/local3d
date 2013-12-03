function [X,mapping] = gen_basedata(dim,S,di,k,affine)
num=S*k; % total sample number

X=cell(S,1); % array of sampled points in each subspace

% random rotation matrix

T=random_rotation(dim);
% [T, ~] = qr(randn(dim));

% random orthogonal matrices
U=cell(S,1);
U{1}=orth(randn(dim,di)); % first U_i, U_1
for i=2:S
    U{i}=T*U{i-1};
end

% random sample matrix
Q=cell(S,1);
for i=1:S
    % Q{i}=randn(di,k);
    Q{i}=rand(di,k);
    Q{i}=Q{i}-0.5;
end

% generate samples
[Xcell]=cellfun(@generate_samples,U,Q,'UniformOutput',false);

% check is affine or not
if affine
    for j=1:S
        X=Xcell{j};
        offset=randn(dim,1);
        X=X+repmat(offset,1,k);
        X=[X;ones(1,k)];
        % need sample more data from the line O-X_i
        knew=3;
        for i=1:num
            Xadd=repmat(X(:,i),1,knew);
            co=rand(1,knew);
            co=repmat(co,dim+1,1);
            Xadd=Xadd.*co;
            X=[X,Xadd];
        end
        Xcell{j}=X;
    end
end
X=cat(2,Xcell{:});


% add noise, outliers
% p=0.97;
% N=rand(size(X))>p;
% N(find(N>0))=1;
% N=N.*rand(size(X));
% X=X+N;

% should not normalize this, because the points is in X;

% for i=1:num
    % X(:,i)=X(:,i)/norm(X(:,i));
% end

mapping=zeros(num,2);
for i=1:num
    mapping(i,:)=[i i];
end

function [X] = generate_samples(U,Q)
    X=U*Q;
