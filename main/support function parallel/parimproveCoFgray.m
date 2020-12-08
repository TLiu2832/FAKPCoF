function J = parimproveCoFgray(lab,h0,ksize0,wsize,gamma,h,ksize,idx1,pab,t)

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu 

C = obtainC2(lab, h0, ksize0, wsize, gamma);
szl  = size(lab);
GW = getnewgw1(h, C, ksize, szl(1), szl(2));
J = lab;
if( isa( J , 'uint8' ) ), J = double( J ); end
nBins = size( pab, 1 );
for i=1:t
    J = picfilter( J, idx1, pab, GW, nBins, ksize);
end
