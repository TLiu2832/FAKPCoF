% MAIN illustrates the use of FAKPCoF in grayscale images.
%
%    This demo shows the typical use of the 'Patch-based Co-occurrence
%    Filter with Fast Adaptive Kernel' in grayscale images. The function 
%    takes as input a grayscale image img.
%
%    December 8th, 2020.
%    Zhonggui Sun, Tingting Liu

% delete(gcp)
% parpool(4)
clear;clc;
% close all;

tic
% Processing Parameters:
% filter - radius
r  = 7;   
t  = 1;
K1 = 30;
sample_rate = 0.1;
gamma=1;
l = 3;
% chek
sigmaI = 0.275;
% star
% sigmaI = 0.267;
h0 = 0.5;
ksize0 = 5;
% the global smoothing parameter
h = 2 * sqrt(r) + 1;     
% the kernel size
ksize = 2*r+1;  

% Read Data:
img = imread('test Checkerboard.bmp');
gray=double(img);
sz  = size(gray);
% improve ( give patch )
fI = fspecial( 'gaussian', [l,l],sigmaI );
Imean = imfilter(gray,fI);
idx1 = round(Imean);
norm(idx1-gray)
% Collect Co-occurrence Statistics:
pab = collectPab0(idx1, ones(sz(1:2)));
pmi = pab./( sum(pab).' * sum(pab) + eps );
% Smooth:  
spmd
  gray1 = gray(:,(labindex-1)*end/4+1:min(labindex*end/4+2*ksize, end));
  Idx1 = idx1(:,(labindex-1)*end/4+1:min(labindex*end/4+2*ksize, end));
 J = parimproveCoFgray(gray1,h0,ksize0,r,gamma,h,ksize,Idx1,pmi,t);
end
J1 = J{1};J2 = J{2};J3 = J{3};J4 = J{4};
J = [J1(:, 1:end-ksize), J2(:, 1+ksize:end-ksize), J3(:, 1+ksize:end-ksize), J4(:, 1+ksize:end)];
toc
J=uint8(J);
% Visualize: 
figure; imshow(img,[]); title('test image');
figure; imshow(J,[]); title('CoF smoothed image');