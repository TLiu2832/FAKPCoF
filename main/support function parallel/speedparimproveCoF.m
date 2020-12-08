function [J1,J2, J3] = speedparimproveCoF(lab1,lab2,lab3,h0,ksize0,wsize,gamma,K1,sample_rate,h,ksize,idx1,pab,t)

%    December 8th, 2020.
%    Zhonggui Sun, Tingting Liu

C = obtainCmean(lab1, h0, ksize0, wsize, gamma, K1,sample_rate);
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
%% Smooth Per Channel:
J1 = Currentchannel( J1, idx1, pab, GW, nBins, ksize);
J2 = Currentchannel( J2, idx1, pab, GW, nBins, ksize);
J3 = Currentchannel( J3, idx1, pab, GW, nBins, ksize);
 
end