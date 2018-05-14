function [imFVal,imFValLarge] = gridFiberCores2(centroids,coreVal,numRows,numCols,subFactor,method)
% B Ozbay
% gridFiberCores
% Reorganize fiber structure into uniform grid
% INPUTS:
% centroids - Array of input image centroids
% coreVal - Array of values corresponding to centroids
% numRows - Number of rows for output image
% numCols - Number of columns for output image
% subFactor - Sub-sampling factor for output image (1 is ~1 pixel for one
% core, 2 is ~2 pixels for one core, etc.)
% method - String that specifies interpolation method to use for
% scatteredInterpolation function:
%           'linear' - linear interpolation
%           'nearest' - nearest neighbor interpolation
%           'natural' - natural nearest neighbor interpolation
%           'coreset' - Set each core location to the nearest pixel
% OUTPUTS:
% imFVal - Image of gridded pixels
% imFValLarge - Image of gridded pixels dilated by mean core distance
subFactor = 1;
numImages = size(coreVal,2);
% Get average distance between cores using nearest neighbor finding
% Find index of nearest neighbors for each core
nIdx = knnsearch(centroids,[centroids(:,1),centroids(:,2)],'K',2);
neighbors = centroids(nIdx(:,2),:);
% Calculate euclidian distance between nearest cores
nDistance = sqrt((centroids(:,1) - neighbors(:,1)).^2 + (centroids(:,2) - neighbors(:,2)).^2);
meanDist = mean(nDistance); % Average over all core distances
sampleDist = round(meanDist/abs(subFactor)); % Divide distance by sub-sampling factor
%%
% numCols = 1000;
% numRows = numCols;
sampleDist = 6;
% Set up uniform matrices to store data
% Use mean distance between nearest cores to create mesh
xInterp = numCols-mod(numCols,sampleDist)+sampleDist;
yInterp = numRows-mod(numRows,sampleDist)+sampleDist;
FxLin = sampleDist-ceil(sampleDist/2):sampleDist:xInterp-ceil(sampleDist/2);
FyLin = sampleDist-ceil(sampleDist/2):sampleDist:yInterp-ceil(sampleDist/2);
% FxLin = 1:sampleDist:numCols;
% FyLin = 1:sampleDist:numRows;
[Fx,Fy] = meshgrid(FxLin,FyLin);
switch lower(method)
    case {'linear','nearest','natural'}
        % Interpolate all data to uniform grid
        centroidsF = [reshape(Fx,[],1),reshape(Fy,[],1)];
        FVal = zeros(size(centroidsF,1),numImages);
        imFOutputTemp = zeros(numRows,numCols);
        FOutput = scatteredInterpolant(centroids,coreVal(:,1));
        FOutput.Method = method;
        FOutput.ExtrapolationMethod = 'none';
        hWait = waitbar(0,...
            sprintf('Interpolating %s with subfactor %.1f...',...
            inputname(2),subFactor));
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
    case 'coreset'
        % For each core, find coordinates in grid that are nearest
        centroidsF = [reshape(Fx,[],1),reshape(Fy,[],1)];
        centroidsFNaN = [reshape(Fx,[],1),reshape(Fy,[],1)];
        numGrid = size(centroidsF,1);
        numCores = size(centroids,1);
        centroidsG = zeros(numCores,2);
        transDist = zeros(numCores,1);
        distances = zeros(numCores,numGrid);
        nearestIdx = zeros(1,numCores);
        numIterations = 2;
        nearbyCdx = zeros(numCores,numIterations);
        idxArray = 1:numCores;
        % Find 2 nearest grid coordinates to each core
        for ii = 1:numCores
            % Get core centroid
            thisCentroid = centroids(ii,:);
            % Find Euclidean distance to all cores
            distances = sqrt(sum(minus(centroidsFNaN,thisCentroid).^2,2 ));
            % Retrieve 3 closest grid locations to the core
            for jj = 1:numIterations
                nearbyCdx(ii,jj) = find(distances==min(distances),1);
                distances(nearbyCdx(ii,jj)) = max(distances);            
            end
        end
 
        for kk=1:numIterations
             % Get the indices of cores that share a nearest grid coordinate
            [~, cdxA] = unique(nearbyCdx(:,1), 'stable');
            repIdx = find(~ismember(1:numCores,cdxA));
            % Assign a next-nearest grid coordinate for duplicate cores
            for ii = 1:length(repIdx)
                for jj = 1:numIterations
                    % Get indices of repeated cores
                    theseRepIdx = find(nearbyCdx(:,1)==nearbyCdx(repIdx(ii),1));
                    if length(theseRepIdx)>1
                        % Find centroid indices of second nearest
                        theseCentroids = centroids(theseRepIdx,:);
                        nextNearestCdx = nearbyCdx(theseRepIdx,2);
                        theseCentroidsF = centroidsF(nextNearestCdx,:);
                        theseDistances = sqrt(sum(minus(theseCentroidsF,theseCentroids).^2,2 ));
                        [~, sortIdx] = sort(theseDistances,'descend');
                        theseRepIdxSorted = theseRepIdx(sortIdx);
                        nearbyCdx(theseRepIdxSorted(2),1) = nearbyCdx(theseRepIdxSorted(2),2);
                    end
                end
            end
        end
        
        centroidsG = centroidsF(nearbyCdx(:,1),:);
        transDist = sqrt(sum(minus(centroidsF(nearbyCdx(:,1),:),centroids).^2,2 ));
        
        %% Compress rows
        centroidsGR = sortrows(centroidsG,1);
        midRow = numCols/2 - mod((numCols/2),sampleDist);
        ii = 1;
        currentRow = 1;
        while ii<=size(centroidsGR,1)
            currentRow = centroidsGR(ii,1)
            rowVector = centroidsGR(ii,2);
            rowStartIdx = ii;
            jj = 1;
            while centroidsGR(ii+1,1) == currentRow
                ii = ii + 1;
                jj = jj + 1;
                rowVector(jj) = centroidsGR(ii,2);
                if ii+1>size(centroidsGR,1)
                    ii=ii-1;
                    currentRow = 0;
                end
            end
%             [~, midIdx] = min(abs(rowVector-midRow));
            midIdx = round(length(rowVector)/2);
            rowVector(midIdx) = midRow;
            leftLength = length(rowVector(1:midIdx-1));
            rightLength = length(rowVector(midIdx+1:end));
            for jj = 1:leftLength
                rowVector(midIdx-jj) = midRow-sampleDist*jj;
            end
            for jj = 1:rightLength
                rowVector(midIdx+jj) = midRow+sampleDist*jj;
            end
            if ii+2>size(centroidsGR,1)
                ii=ii+1;
            end
            centroidsGR(rowStartIdx:ii,2) = rowVector;
            ii=ii+1;
        end
                
                
                
        %%
        figure(10); clf; hold on;
        colorbar; colormap(parula);
        set(gcf,'Units','Pixels','Position',[750 100 800 700]);
        set(gca,'Units','Pixels','color','k');
        gcaPos = get(gca,'Position');
        % scatter(centroids(:,2),centroids(:,1),'.b');
        maxMov = 8;
        transDistPlot = transDist;
        % transDistPlot(transDist>maxMov) = maxMov;
        scatter(centroidsF(:,2),centroidsF(:,1),5,'r','s','filled');
        scatter(centroidsGR(:,2),centroidsGR(:,1),10,transDistPlot,'s','filled');
        xlim([0 numCols]); ylim([0 numRows]);
        dispSize = sqrt(numGrid)*3.8;
        set(gca,'Position',[gcaPos(1)/2, gcaPos(2)/2, dispSize, dispSize/1.5]);

        meanTrans = mean(transDist)
    otherwise
        error('Invalid method');
end
%%


% 
% % Find the nearest grid coordinate
% nearestIdx(ii) = find(distances==min(distances),1);
% nearestIdx(ii)
% % Record the distance the core had to move
% transDist(jj) = distances(nearestIdx(ii));
% % Assign this core to the chosen grid position
% centroidsG(jj,:) = centroidsFNaN(nearestIdx(ii),:);
% % Remove this grid position from the list
% centroidsFNaN(nearestIdx(ii),:) = NaN(1,2);

