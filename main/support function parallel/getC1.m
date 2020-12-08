function C = getC1(zx, zy, I, wsize, lambda, gamma)

% These codes are mainly based on the "Kernel regression for image
% processing and reconstruction" paper: 
% Distributed under:
% http://alumni.soe.ucsc.edu/~htakeda/KernelToolBox.htm

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu

[N,M] = size(zx);
C = zeros(2, 2, N, M);
if mod(wsize, 2) == 0
    wsize = wsize + 1;
end
win = floor(wsize / 2);

K = fspecial('disk', win);
K = K ./ K(win+1, win+1);

% mirroring
zx = EdgeMirror(zx, [win, win]);
zy = EdgeMirror(zy, [win, win]);

for i = 1 : N
    for j = 1 : M

        if I(i,j) == 0
            continue;
        end   

        gx = zx(i:i+wsize-1, j:j+wsize-1) .* K;
        gy = zy(i:i+wsize-1, j:j+wsize-1) .* K;

        G = [gx(:), gy(:)];  
        [~,s,v] = svd(G, 0);

        S1 = (s(1,1) + lambda) / (s(2,2) + lambda);
        S2 = (s(2,2) + lambda) / (s(1,1) + lambda);

        tmp = (S1 * v(:,1) * v(:,1)' + S2 * v(:,2) * v(:,2)')  * gamma;
        C(1:2, 1:2, i, j) = tmp;

    end
end


