function fim = Currentchannel(im2filter,im2collect,pab,GW,nBins,ksize)
    
%    December 8th, 2020.
%    Zhonggui Sun, Tingting Liu

cim = im2filter;
mvals = double( im2collect );
szfi=size(im2filter);
% Per level results:
W = zeros( szfi(1), szfi(2), nBins );
S = zeros( szfi(1), szfi(2), nBins );
FW=W;
FS=S;
% Per level smoothing:
m = mat2cell( pab,ones( nBins, 1 ));
for iLevel = 0:( nBins - 1 )
    % weight per level
    wpl = reshape( m{ iLevel + 1 }( mvals + 1 ), szfi(1:2) );
    
    yw = wpl;
    ys = wpl .* cim;
    %  tic
    [FS( :, :, iLevel + 1 ),FW( :, :, iLevel + 1 )] = ...
        newgaussian1( cim, GW, ksize, ys, yw, mvals, iLevel );
    %  toc
    W( :, :, iLevel + 1 ) = FW( :, :, iLevel + 1 );
    S( :, :, iLevel + 1 ) = FS( :, :, iLevel + 1 );
    %         waitbar(iLevel/( nBins - 1 ));
end
fim = sum( S, 3 )./ (  sum( W, 3 ) + eps );