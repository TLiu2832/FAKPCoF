function [pab, s] = collectPab0( im2collect, mask2collect, nBins, f_cooc )
%COLLECTPMI collect co-occurrence statistics
%   
% Syntax:
%   pab = collectPab( im2collect, mask2collect, nBins, f_cooc );
%   [pab, s] = collectPab( im2collect, mask2collect, nBins, f_cooc );
%
% Input: 
%   im2collect: collect co-occurrence statistics from this image        
%   mask2collect: binary mask, 0 - skip this pixel; 1 - keep this pixel     [defult: all ones]
%   nBins: im2collect can have values between 0 - (nBins - 1)               [default: 256]
%   f_cooc: spatial filter for collecting co-occurrence                     [default: Gaussian]
%
% Output: 
%   pab: co-occurrence statistics
%   s: the total weight of p(a, b) before we normalize it to distribution 
%
% These codes are mainly based on the "Co-occurrence Filter" paper:
% Distributed under the GNU GPL license:
% http://www.gnu.org/copyleft/gpl.html

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu 

%% Default Parameters:
if( isa( im2collect, 'uint8' ) ), im2collect = double( im2collect ); end
if( nargin < 2 ), mask2collect = ones( size( im2collect ) ); end
if( nargin < 3 ), nBins = 256; end
if( nargin < 4 ), f_cooc = fspecial( 'gaussian', 15, sqrt(15) ); end

%% Slice image to levels, collect co-occurrence per level:
cref = im2collect;
pab = zeros( nBins, nBins );

for iLevel = 0:( nBins - 1 )
    
    w = conv2( double( cref == iLevel ) .* mask2collect, f_cooc, 'same' ).* mask2collect;
    pab( iLevel + 1, : ) = accumarray( cref(:) + 1 , w(:), [ nBins, 1 ] );
    
end

%% Normalize: 
pab = pab + pab.';
s = sum( pab(:));
pab = pab ./ s;


end

