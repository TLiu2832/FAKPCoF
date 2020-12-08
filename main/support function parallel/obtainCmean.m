function [C] = obtainCmean(y,h0, ksize0, ksize, gamma,K1,Srate)

% These codes are mainly based on the "Fast High-Dimensional Bilateral and
% Nonlocal Means Filtering" paper: 
% Distributed under:
% https://github.com/pravin1390/FastHDFilter.
%    December 8th, 2020.
%    Zhonggui Sun, Tingting Liu

% tic
sz = size( y );
% the regularization for the elongation parameter
lambda=1; 
[~, zx1, zx2] = ckr2(y, h0, ksize0);
% toc
tic
[N,M] = size(zx1);
C = zeros(2, 2, N, M);
if mod(ksize, 2) == 0
    ksize = ksize + 1;
end
win = floor(ksize / 2);

K = fspecial('disk', win);
K = K ./ K(win+1, win+1);

% mirroring
zx1 = EdgeMirror(zx1, [win, win]);
zx2 = EdgeMirror(zx2, [win, win]);
% toc
% tic
%% Convert image slices to vectors,merging zx,zy
patchesx = im2col(zx1,[ksize,ksize],'sliding');
patchesy = im2col(zx2,[ksize,ksize],'sliding');
patches = [patchesx;patchesy];
patches = patches';
%% PCA dimensionality reduction for vectors
pca_outdim = 6;
H = bsxfun(@minus, patches, mean(patches));
[Evec , ~] = eig( H'*H );
Evec = flip(Evec, 2);
H1 = H*Evec(:,1:pca_outdim);
% Transpose the column vector
% Eval = flipdim(diag(Eval),1)';
% Eval = Eval(1:pca_outdim);
%% The points of each vector after dimension reduction are sampled
MN = size(H1, 1);
mnpoint=randperm(MN, round(Srate * MN));
H2=H1(mnpoint,:);
% toc
% tic
%% Clustering the obtained vectors by k-means
% tic
% [~, cc] = kmeans(H2, K1);
% [~, cc] = kmeans(H2, K1,'MaxIter',500);
[cc, ~] = kmeans_recursive(H2, K1);
% toc
%% Get class center and class label
cc=sortrows(cc);
rr=reshape(H1,[N M pca_outdim]);
qim = zeros( sz(1:2) );
min_err = inf( sz(1:2) );
%% Categorize the classes of points throughout the image
for k = 1:K1
    [min_err, idx1] =min( cat(3, min_err, ...
    sum((rr - repmat(permute(cc(k,:), [3 1 2] ), sz(1:2))).^2, 3)), [], 3);
    qim(idx1 == 2) = k;
end

%% Find the location of the cluster center in the original image
NUcentre=zeros(K1,1);
for k = 1:K1
    DD=bsxfun(@minus, H2, cc(k, :) );
    [~, idx] =min(sum((DD.^2),2));
    NUcentre(k)=idx;
end
tt = mnpoint(NUcentre);
cen = patches(tt,:);
%% The class center is decomposed by SVD to obtain the C matrix
% tic
C2=zeros(2,2,K1);
Centre=cen';
len=size(Centre,1)/2;
for i = 1 : K1
    
    gx = Centre(1:len,i) ;
    gx=gx.*K(:);
    gy = Centre(len+1:end,i) ;
    gy=gy.*K(:);
    G = [gx(:), gy(:)];
    
    [~,s,v] = svd(G, 0);
    
    S(1) = (s(1,1) + lambda) / (s(2,2) + lambda);
    S(2) = (s(2,2) + lambda) / (s(1,1) + lambda);
    
    tmp = (S(1) * v(:,1) * v(:,1)' + S(2) * v(:,2) * v(:,2)')  * gamma;
    C2(1:2, 1:2, i) = tmp;
    
end
% toc
%% Points in the same class share information about the class center of the same class
% tic
for i = 1 : N
    for j = 1 : M
        ff=qim(i,j);
        C(:,:,i,j)=C2(:,:,ff);
    end
end
toc
% toc