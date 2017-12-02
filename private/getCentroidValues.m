function [coreVal] = getCentroidValues(imFiber,centroids,sampleRad)
%GETCENTROIDVALUES Get values of the fibers at each centroid location

numRows = size(imFiber,1);
numCols = size(imFiber,2);
numImages = size(imFiber,3);
numCores = size(centroids,1);

% Create mask for sampling cores
[rr,cc] = meshgrid(1:sampleRad*2+1);
circMask = sqrt((rr-sampleRad-1).^2+(cc-sampleRad-1).^2)<=sampleRad-0.5;
[mRow,mCol] = find(circMask);
mRow = mRow - (sampleRad + 1);
mCol = mCol - (sampleRad + 1);

% Pad centroid array
centroidPad = centroids;
centroidPad(centroidPad<sampleRad+1) = sampleRad+1;
centroidPad(centroidPad(:,1)>(numRows-sampleRad),1) = numRows-sampleRad;
centroidPad(centroidPad(:,2)>(numCols-sampleRad),2) = numCols-sampleRad;


% % - TEST - % Check sampling of original cores
% imTester = zeros(numRows,numCols);
% imTestInput = max(imFiber,[],3);
% % figure(1);
% % imshow(circMask,'initialmagnification',5000)
% %
% for jj = 1:numCores
%     cRow = mRow + centroidPad(jj,2);
%     cCol = mCol + centroidPad(jj,1);
%     imTester(sub2ind(size(imTester),cRow,cCol)) = imTester(sub2ind(size(imTester),cRow,cCol)) + 1;
% end
% figure(107); clf; imshowpair(imTester/2,imTestInput./max(imTestInput(:)),'scaling','none','colorchannels',[2 1 2]);
% % pause;
% % - TEST - %

coreVal = zeros(numCores,numImages);
topPixels = ceil(length(circMask(circMask))); % Number of pixels from region to sample
hWait = waitbar(0, sprintf('Processing stack...'));
for ii = 1:numImages
    imCurrent = imFiber(:,:,ii);
    for jj = 1:numCores
        cRow = mRow + centroidPad(jj,2);
        cCol = mCol + centroidPad(jj,1);
        tempMaskVals = imCurrent(sub2ind([numRows,numCols],cRow,cCol));
        tempMaskVals = sort(tempMaskVals(:),'descend');
        coreVal(jj,ii) = mean(tempMaskVals(1:topPixels));
    end
    waitbar(ii/numImages,hWait);
end
close(hWait);



