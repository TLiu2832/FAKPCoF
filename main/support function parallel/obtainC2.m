function C = obtainC2(y, h0, ksize0, ksize, gamma)

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu 

sz = size( y );
m=sz(1);
n=sz(2);
% the regularization for the elongation parameter
lambda=1; 
[~, zx1, zx2] = ckr2(y, h0, ksize0);
% tic
C = getC1(zx1, zx2, ones(m,n), ksize, lambda, gamma);
% toc