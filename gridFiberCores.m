function [imFVal,imFValLarge] = gridFiberCores(centroids,coreVal,numRows,numCols,subFactor,method)
% gridFiberCores - B. Ozbay (5/11/2017)
% [imFVal,imFValLarge] = gridFiberCores(centroids,coreVal,numRows,numCols,subFactor,method)
% Reorganize fiber structure into uniform grid with specified interpolation
% method.
%
% INPUTS:
% centroids - Nx2 array of input image centroids pixel coordinates for N
% fiber cores
% coreVal - NxM matrix of values corresponding to N centroids for an image
% series of M frames
% numRows - Number of rows for output image
% numCols - Number of columns for output image
% subFactor - Sub-sampling factor for output image (1 is ~1 pixel for one
% core, 2 is ~2 pixels for one core, etc.)
% method - String that specifies interpolation method to use for
% scatteredInterpolation function:
%           'linear' - linear interpolation
%           'nearest' - nearest neighbor interpolation
%           'natural' - natural nearest neighbor interpolation
%
% OUTPUTS:
% imFVal - Image of gridded pixels
% imFValLarge - Image of gridded pixels dilated by mean core distance

numImages = size(coreVal,2);

% Get average distance between cores using nearest neighbor finding
% Find index of nearest neighbors for each core
nIdx = knnsearch(centroids,[centroids(:,1),centroids(:,2)],'K',2);
neighbors = centroids(nIdx(:,2),:);
% Calculate euclidian distance between nearest cores
nDistance = sqrt((centroids(:,1) - neighbors(:,1)).^2 + (centroids(:,2) - neighbors(:,2)).^2);
meanDist = mean(nDistance); % Average over all core distances
sampleDist = round(meanDist/abs(subFactor)); % Divide distance by sub-sampling factor

% Set up uniform matrices to store data
% Use mean distance between nearest cores to create mesh
xInterp = numCols-mod(numCols,sampleDist);
yInterp = numRows-mod(numRows,sampleDist);
FxLin = sampleDist-ceil(sampleDist/2):sampleDist:xInterp-ceil(sampleDist/2);
FyLin = sampleDist-ceil(sampleDist/2):sampleDist:yInterp;
[Fx,Fy] = meshgrid(FxLin,FyLin);

% Interpolate all data to uniform grid
centroidsF = [reshape(Fx,[],1),reshape(Fy,[],1)];
FVal = zeros(size(centroidsF,1),numImages);
imFOutputTemp = zeros(numRows,numCols);
FOutput = scatteredInterpolant(centroids,coreVal(:,1));
FOutput.Method = method;
FOutput.ExtrapolationMethod = 'none';
hWait = waitbar(0, sprintf('Interpolating %s with subfactor %.1f...',inputname(2),subFactor));
for ii = 1:numImages
    FOutput.Values = coreVal(:,ii);
    FOutTemp = FOutput(Fx,Fy);
    imFOutputTemp(Fy(:,1),Fx(1,:),ii) = FOutTemp;
    FVal(:,ii) = reshape(FOutTemp,[],1);
    waitbar(ii/numImages,hWait);
end
close(hWait);

FVal(isnan(FVal)) = 0;
imFOutputTemp(isnan(imFOutputTemp)) = 0;
imFVal = reshape(FVal,length(FyLin),length(FxLin),[]);
if nargout>1
    imFValLarge = imdilate(imFOutputTemp,strel('square',sampleDist));
end

