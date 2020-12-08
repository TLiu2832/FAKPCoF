function [zs,zw] = newgaussian1(y, GW, ksize, ys, yw, mvals, iLevel)
% [zs,zw] = newgaussian(y, h, C, ksize, ys, yw, mvals, iLevel )

% The second order steering kernel regression function for regularly sampled
% data.
%
% [RETURNS]
% z     : the estimated image
% zs   : the estimated gradient image along the x1 direction (vertical
%        direction)
% zw   : the estimated gradient image along the x2 direction (horizontal
%        direction)
%
% [PARAMETERS]
% y     : the input image with co-occurence matrix
% GW     : the global smoothing parameter
% ksize : the size of the kernel (ksize x ksize, and "ksize" must be
%         an odd number)

% These codes are mainly based on the "Kernel regression for image
% processing and reconstruction" paper: 
% Distributed under:
% http://alumni.soe.ucsc.edu/~htakeda/KernelToolBox.htm

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu 

% tic
% Get the oritinal image size
[N, M] = size(y);
if mod(ksize, 2) == 0
    ksize = ksize + 1;
end
% Pixel sampling positions
radius = (ksize -1) / 2;
zs=zeros(N, M);
zw=zeros(N, M);
W=zeros(ksize,ksize);
% Mirroring
ys = padarray(ys, [radius, radius],'symmetric');
yw = padarray(yw, [radius, radius],'symmetric');
% adaptive kernel
for i = 1 : N
%     tic
    for j = 1 : M 
        % Neighboring samples to be taken account into the estimation
        if double( mvals(i,j) == iLevel )
            ysp = ys(i:i+ksize-1, j:j+ksize-1);
            ywp = yw(i:i+ksize-1, j:j+ksize-1);
            % Estimate the pixel values at (i,j)
%             zs(i,j)   = GW(:,i,j)' * ysp(:);
%             zw(i,j)   = GW(:,i,j)' * ywp(:);
            W(:,:) = GW(:,:,i,j);
            zs(i,j)   = W(:)' * ysp(:);
            zw(i,j)   = W(:)' * ywp(:);
        end  
    end
%    toc 
end