% MAIN illustrates the use of FAKPCoF in color images.
%
%    This demo shows the typical use of the 'Patch-based Co-occurrence
%    Filter with Fast Adaptive Kernel' in grayscale images. The function 
%    takes as input a color image img.

%    December 8th, 2020.
%    Zhonggui Sun, Tingting Liu

% parpool(4)
clc;close all;
clear;
% tic

% Processing Parameters:

% number of bins (for the quantization)
nBins   = 32;   
% filter - radius
r  = 7;    
t  = 3;
K1 = 30;
sample_rate = 0.1;
gamma=0.8;
h0 = 2.4;
ksize0 = 15;
% the global smoothing parameter
h=2 * sqrt(r) + 1;  
% the kernel size
ksize = 2*r+1;  
l = 3;
sigmaI = 0.3055;
% Read Data:
img = imread('Bishapur zan.png');
lab = rgb2lab(img);
sz  = size(img);
% improve ( give weighted )
fI = fspecial( 'gaussian', [l,l], sigmaI );
Imean = imfilter(img,fI);
imgmean = round(Imean);
[idx1,cc] = quantize(imgmean, nBins);
% Collect Co-occurrence Statistics:
pab = collectPab0(idx1, ones(sz(1),sz(2)), nBins);
pmi = pab./( sum(pab).' * sum(pab) + eps );
% Smooth:
%% lab
lab1=lab(:,:,1);
lab2=lab(:,:,2);
lab3=lab(:,:,3);
spmd
    Lab1 = lab1(:,(labindex-1)*end/4+1:min(labindex*end/4+2*ksize, end));
    Lab2 = lab2(:,(labindex-1)*end/4+1:min(labindex*end/4+2*ksize, end));
    Lab3 = lab3(:,(labindex-1)*end/4+1:min(labindex*end/4+2*ksize, end));
    Idx1 = idx1(:,(labindex-1)*end/4+1:min(labindex*end/4+2*ksize, end));
    [J1,J2,J3] = speedparimproveCoF(Lab1,Lab2,Lab3,h0,ksize0,r,gamma,K1,sample_rate,h,ksize,Idx1,pmi,t);
end
J11 = J1{1};J12 = J1{2};J13 = J1{3};J14 = J1{4};
J21 = J2{1};J22 = J2{2};J23 = J2{3};J24 = J2{4};
J31 = J3{1};J32 = J3{2};J33 = J3{3};J34 = J3{4};
J(:,:,1) = [J11(:, 1:end-ksize), J12(:, 1+ksize:end-ksize), J13(:, 1+ksize:end-ksize), J14(:, 1+ksize:end)];
J(:,:,2) = [J21(:, 1:end-ksize), J22(:, 1+ksize:end-ksize), J23(:, 1+ksize:end-ksize), J24(:, 1+ksize:end)];
J(:,:,3) = [J31(:, 1:end-ksize), J32(:, 1+ksize:end-ksize), J33(:, 1+ksize:end-ksize), J34(:, 1+ksize:end)];
% toc
J = lab2rgb(J);
% Visualize:
figure; imshow(img); title('input image');
figure; imshow(J,[]); title('improve CoF smoothed image');