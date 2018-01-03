function [correctorVal] = getCorrectorValues(imFlatVal, centroids, cropIn)
% B Ozbay 11/15/2017
% getCorrectorVals
% Inputs:
% imFlatVal - Image of flat values translated to input image
% centroids - Array of centroids from input image
% cropIn - Struct of cropping positions
% Outputs:
% correctorVal - Array of values corresponding to centroids

% Create shading filter
intShadFilter = imgaussfilt(imFlatVal*1,30);
imFlatValShad = imFlatVal./intShadFilter;
imFlatValShad = imFlatValShad/max(imFlatValShad(:));
% Create corrector image with inverse flat
imCorrector = mean(imFlatVal(:))*(imFlatValShad.^-1);
% Crop imCorrector to size of input image
% imCorrector = imCorrector(cropIn.StartY:cropIn.EndY,cropIn.StartX:cropIn.EndX);
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
