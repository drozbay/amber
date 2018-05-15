function [correctorVal] = getCorrectorValues(imFlatVal, centroids)
% getCorrectorVals
% [correctorVal] = getCorrectorValues(imFlatVal, centroids, cropIn)
% Gets correction values from imaging a flat fluorescence sample to correct
% for inhomogeneity in fiber images
%
% INPUTS:
% imFlatVal -   Image of flat values registered to input image
% centroids -   Nx2 array of N centroid pixel locations, aligned with input
%               image
% OUTPUTS:
% correctorVal - Nx1 array of values corresponding to N fiber cores
%
% Author: Baris N. Ozbay, University of Colorado Denver
% Version: 1.0

% Create shading filter
intShadFilter = imgaussfilt(imFlatVal*1,30);
imFlatValShad = imFlatVal./intShadFilter;
imFlatValShad = imFlatValShad/max(imFlatValShad(:));
% Create corrector image with inverse flat
imCorrector = mean(imFlatVal(:))*(imFlatValShad.^-1);
% Sample corrector values back into vector
correctorVal = getCentroidValues(imCorrector,centroids,1);
% Remove erroneous values
correctorVal(isinf(correctorVal))=0;
correctorVal(isnan(correctorVal))=0;
% Limit maximum corrector value to prevent spiking signal
% For max value: Set as median value plus factor that depends on 
% median absolute deviation (mad) value
madCorrectorVal = median(abs(correctorVal-median(correctorVal)));
maxCorrectorVal = median(correctorVal)+madCorrectorVal;
correctorVal(correctorVal>maxCorrectorVal) = maxCorrectorVal;

end

