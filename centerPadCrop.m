function [outputImage] = centerPadCrop(inputImage,xSize,ySize)
% centerPadCrop
% [outputImage] = centerPadCrop(inputImage,xSize,ySize)
% Performs a central crop on an image.
%
% INPUTS:
% inputImage - Image to be padded/cropped
% xSize, ySize - Size of output image
%
% OUTPUTS:
% outputImage - Cropped image from center of inputImage
%
% Author: Baris N. Ozbay, University of Colorado Denver
% Version: 1.0

inputXSize = size(inputImage,2);
inputYSize = size(inputImage,1);
inputZSize = size(inputImage,3);

bgPad = -[inputYSize-ySize, inputXSize-xSize]/2;
bgPad(bgPad<0) = 0;
bgCrop = [inputYSize-ySize, inputXSize-xSize]/2;
bgCrop(bgCrop<0) = 0;

padStart = floor(bgPad)+1;
padEnd = padStart + [inputYSize,inputXSize]-1;
imPadded = zeros(ySize,xSize,inputZSize);
imPadded(padStart(1):padEnd(1),padStart(2):padEnd(2),:) = inputImage;

cropStart = floor(bgCrop)+1;
cropEnd = cropStart + [ySize,xSize]-1;
imCropped = imPadded(cropStart(1):cropEnd(1),cropStart(2):cropEnd(2),:);

outputImage = imCropped;
