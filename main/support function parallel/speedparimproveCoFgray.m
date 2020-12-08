function J = speedparimproveCoFgray(lab,h0,ksize0,wsize,gamma,K1,sample_rate,h,ksize,idx1,pab,t)

%    December 8th, 2020.
%    Zhonggui Sun, Tingting Liu

C = obtainCmean(lab, h0, ksize0, wsize, gamma, K1,sample_rate);
szl  = size(lab);
GW = getnewgw1(h, C, ksize, szl(1), szl(2));
J = lab;
if( isa( J , 'uint8' ) ), J = double( J ); end
nBins = size( pab, 1 );
for i=1:t
    J = picfilter( J, idx1, pab, GW, nBins, ksize);
end