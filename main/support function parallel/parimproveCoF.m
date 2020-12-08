function [J1,J2, J3] = parimproveCoF(lab1,lab2,lab3,h0,ksize0,wsize,gamma,h,ksize,idx1,pab,t)

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu 

C = obtainC2(lab1, h0, ksize0, wsize, gamma);
szl  = size(lab1);
GW = getnewgw1(h, C, ksize, szl(1), szl(2));
J1 = lab1;
J2 = lab2;
J3 = lab3;
if( isa( J1 , 'uint8' ) ), J1 = double( J1 ); end
if( isa( J2 , 'uint8' ) ), J2 = double( J2 ); end
if( isa( J3 , 'uint8' ) ), J3 = double( J3 ); end
nBins = size( pab, 1 );
for i=1:t
    J1 = picfilter( J1, idx1, pab, GW, nBins, ksize);
    J2 = picfilter( J2, idx1, pab, GW, nBins, ksize);
    J3 = picfilter( J3, idx1, pab, GW, nBins, ksize);
end