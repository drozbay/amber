function [imFiberVal, imFiberBW] = makeFiberImage(centroids,coreValues,numRows,numCols,SE)
%MAKEFIBERIMAGE 

numImages = size(coreValues,2);
numCores = size(centroids,1);

imFiberVal = zeros(numRows,numCols,numImages);

if numImages>1
    for ii = 1:numImages
        for jj = 1:numCores
            imFiberVal(centroids(jj,2),centroids(jj,1),ii)=coreValues(jj,ii);
        end
    end
else
    for jj = 1:numCores
        imFiberVal(centroids(jj,2),centroids(jj,1))=coreValues(jj,1);
    end
end

imFiberVal = imdilate(imFiberVal,SE);
imFiberBW = imFiberVal;
imFiberBW(imFiberBW>0)=1;
% imFiberBW = imbinarize(imFiberBW);