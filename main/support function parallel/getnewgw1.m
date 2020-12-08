function GW = getnewgw1(h, C, ksize, N, M )

% These codes are mainly based on the "Kernel regression for image
% processing and reconstruction" paper: 
% Distributed under:
% http://alumni.soe.ucsc.edu/~htakeda/KernelToolBox.htm

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu 

% Get the oritinal image size
if mod(ksize, 2) == 0
    ksize = ksize + 1;
end
% Pixel sampling positions
radius = (ksize -1) / 2;
[x2, x1] = meshgrid(-radius : radius, -radius : radius);
C11 = zeros(N, M);
C12 = zeros(N, M);
C22 = zeros(N, M);
GW=zeros(ksize,ksize,N,M);

for n = 1 : N
    for m = 1 : M
C11(n,m) = C(1,1,n,m);
C12(n,m) = C(1,2,n,m);
C22(n,m) = C(2,2,n,m);
    end
end

% Mirroring
C11 = padarray(C11, [radius, radius],'symmetric');
C12 = padarray(C12, [radius, radius],'symmetric');
C22 = padarray(C22, [radius, radius],'symmetric');

% adaptive kernel
for i = 1 : N
    for j = 1 : M

        % compute the weight matrix
        tt = x1 .* (C11(i:i+ksize-1, j:j+ksize-1) .* x1...
            + C12(i:i+ksize-1, j:j+ksize-1) .* x2)...
            + x2 .* (C12(i:i+ksize-1, j:j+ksize-1) .* x1...
            + C22(i:i+ksize-1, j:j+ksize-1) .* x2);

%         tt = x1 .* (C11(i+radius, j+radius) * x1...
%             + C12(i+radius, j+radius) * x2)...
%             + x2 .* (C12(i+radius, j+radius) * x1...
%             + C22(i+radius, j+radius) * x2);  

        
        W = exp(-(0.5/h^2) * tt) ;
        W(W<eps*max(W(:))) = 0;
        sumW = sum(W(:));
        if sumW ~= 0
            W  = W/sumW;
        end
      
        GW(:,:,i,j) = W;   
    end
end