function [z, zx1, zx2] = ckr2(y, h, ksize)
% [CKR2]
% The second order classic kernel regression function for regularly sampled
% data.
%
% [USAGE]
% [z, zx1, zx2] = ckr2(y, h,ksize)
%
% [RETURNS]
% z     : the estimated image
% zx1   : the estimated gradient image along the x1 direction (vertical
%        direction)
% zx2   : the estimated gradient image along the x2 direction (horizontal
%        direction)
%
% [PARAMETERS]
% y     : the input image
% h     : the global smoothing parameter
% ksize : the size of the kernel (ksize x ksize, and "ksize" must be
%         an odd number)

% These codes are mainly based on the "Kernel regression for image
% processing and reconstruction" paper: distributed under:
% http://alumni.soe.ucsc.edu/~htakeda/KernelToolBox.htm

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun, Tingting Liu


% Get the oritinal image size
[N, M] = size(y);

% Initialize the return parameters
z = zeros(N, M);
zx1 = zeros(N, M);
zx2 = zeros(N, M);

% Create the equivalent kernels
radius = (ksize - 1) / 2;
[x2, x1] = meshgrid(-radius : 1 : radius, -radius : 1 : radius);
A = zeros(6, ksize^2);

% The feture matrix
X = [ones(ksize^2,1), x1(:), x2(:), x1(:).^2, x1(:).*x2(:), x2(:).^2];
% The weight matrix (Gaussian kernel function)
tt = x1.^2 + x2.^2;
G = exp(-(0.5/h^2) * tt);

% Equivalent kernel
XG = [X(:,1).*G(:), X(:,2).*G(:), X(:,3).*G(:),...
    X(:,4).*G(:), X(:,5).*G(:), X(:,6).*G(:)];
A(:,:) = inv(X' * XG) * (XG');

% Mirroring the input image
y = EdgeMirror(y, [radius, radius]);

% Estimate an image and its first gradients with pixel-by-pixel
for n = 1 : N
    for m = 1 : M
        
        % Neighboring samples to be taken account into the estimation
        yp = y(n:n+ksize-1, m:m+ksize-1);
        
        % Estimate the pixel values at (nn,mm)
        z(n,m)   = A(1,:) * yp(:);
        zx1(n,m) = A(2,:) * yp(:);
        zx2(n,m) = A(3,:) * yp(:);
                
    end
end




