function [qim,cc] = quantize(input_im, K, sample_rate)
% QUANTIZE k-means quantization 
%   
% Syntax:
%   qim = quantize(input_im, K, sample_rate)
%
% Input: 
%   input_im: image to be quantized
%   K:  Number of clusters for K-Means                                     [default:256] 
%   sample_rate: rate at which the image is sampled prior to K-Means       [default: 0.1]
%
% Output: 
%   qim: quantized iamge 
%
% These codes are mainly based on the "Co-occurrence Filter" paper:
% Distributed under the GNU GPL license:
% http://www.gnu.org/copyleft/gpl.html

%   Version: 1.0
%   Date: December 8th, 2020
%   Zhonggui Sun,Tingting Liu 

% Default Arguments:
if( nargin < 2 ), K = 256; end
if( nargin < 3 ), sample_rate = 0.1; end

% Reformat the Image 
input_im = rgb2lab(input_im);
lab = reshape(input_im, [], 3); 
nPoints = size(lab, 1);
lab = lab(randperm(nPoints, round(sample_rate * nPoints)), : );

% K-Means:
warning('off' ); % diable l means warnings:

%%
% K cluster centroid locations are stored
[~, cc] = kmeans(lab, K); 
warning('on');

% Qunatized Image
sz = size( input_im ); 
qim = zeros( sz(1:2) );
min_err = inf( sz(1:2) );
for k = 1:K
    
    [min_err, idx] = ...
        min( cat(3, min_err, ...
        sum((input_im - repmat(permute(cc(k, :), [3, 1, 2] ), sz(1:2))).^2, 3) ), [], 3 ); 
    qim(idx == 2) = k-1;
    
end


