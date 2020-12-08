function fim = improvecoF2( im2filter, im2collect, pab, ksize, GW)
% These codes are mainly based on the "Co-occurrence Filter" paper:
% Distributed under the GNU GPL license:
% http://www.gnu.org/copyleft/gpl.html

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu 

%% Default Parameters:
if( isa( im2filter , 'uint8' ) ), im2filter = double( im2filter ); end
if( size( im2filter, 3 ) == 3 * size( pab, 3 ) ); pab = repmat( pab, [1, 1, 3] ); end
if( size( im2filter, 3 ) == 3 * size( im2collect, 3 ) ); im2collect = repmat( im2collect, [1, 1, 3] ); end
%% Processing Parameters:
sz = size( im2filter ); if( length(sz) == 2 ), sz(3) = 1; end
nBins = size( pab, 1 );
fim = zeros( sz );

%% Smooth Per Channel:
for iChannel = 1:sz(3)
    % Current channel:
    cim = im2filter( :, :, iChannel );
    mvals = double( im2collect( :, :, iChannel ) );
    
    % Per level results:
    W = zeros( sz(1), sz(2), nBins );
    S = zeros( sz(1), sz(2), nBins );
    FW=W;
    FS=S;
    % Per level smoothing:
    m = mat2cell( pab( :, :, iChannel ), ones( nBins, 1 )  );
    for iLevel = 0:( nBins - 1 )
        % weight per level
        wpl = reshape( m{ iLevel + 1 }( mvals + 1 ), sz(1:2) );            	
        yw = wpl;
        ys = wpl .* cim;
        %  tic
        [FS( :, :, iLevel + 1 ),FW( :, :, iLevel + 1 )] = newgaussian1( cim, GW, ksize, ys, yw, mvals, iLevel );
        %  toc
        W( :, :, iLevel + 1 ) = FW( :, :, iLevel + 1 );
        S( :, :, iLevel + 1 ) = FS( :, :, iLevel + 1 );
%         waitbar(iLevel/( nBins - 1 ));
    end
    fim( :, :, iChannel ) = sum( S, 3 )./ (  sum( W, 3 ) + eps );
end
end