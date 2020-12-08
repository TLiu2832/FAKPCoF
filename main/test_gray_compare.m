clc;
clear;
close all

% Processing Parameters:
r  = 7;    % filter - radius
img = imread('test Checkerboard.bmp');
% img = imread('test Stars Image.bmp');
gray = double(img);
t = 1;
% t = 5;
%%
idx0 = gray;
res0 = gray;
res = gray;
sz = size(gray);
%%
% improve ( give patch )
l = 3;
sigmaI = 0.275;%chek
% sigmaI = 0.267;%star
fI = fspecial( 'gaussian', [l,l], sigmaI );
imgmean = imfilter(gray,fI);
idx1 = round(imgmean);
% idx=idx1-idx0;
% norm(idx)
gamma = 1;
h0 = 0.5;
ksize0 = 5;
wsize=r;
h = 2 * sqrt(r) + 1;     % the global smoothing parameter
ksize = 2*r+1;  % the kernel size
%%
% Collect Co-occurrence Statistics:
pab0 = collectPab0(idx0, ones(sz(1:2)));
pab = collectPab0(idx1, ones(sz(1:2)));
pmi0 = pab0./( sum(pab0).' * sum(pab0) + eps );
pmi = pab./( sum(pab).' * sum(pab) + eps );

% Smooth:  
f = fspecial('gaussian', 2 * r, 2 * sqrt(r) + 1);
C = obtainC2(res, h0, ksize0, wsize,gamma);
GW = getnewgw1(h, C, ksize,sz(1),sz(2));
J=res;
for i=1:t
res0 = coFilter0(res0, idx0, pmi0, f ); 
tic

J = improvecoF2( J, idx1, pmi, ksize, GW);
toc
end
res0=uint8(res0);
J=uint8(J);
figure; imshow(img,[]); title('test image');
figure; imshow(res0,[]); title(' CoF0 smoothed image');
figure; imshow(J,[]); title('CoF smoothed image');
