clear all;
close all;
% Open previous configuration settings, if they exist. Otherwise create new
% configuration file with default settings.
if exist('fcp_config.mat', 'file') == 2
    load('fcp_config.mat');
else
    pnameStartI = 'C:\';
    pnameStartF = 'C:\';
    pnameStartBG = 'C:\';
    threshInput = 0.5;
    threshFlat = 0.5;
    startX = 0;
    startY = 0;
    save('fcp_config.mat','pnameStartI','pnameStartF','threshInput',...
        'threshFlat','startX','startY');
end
% Get input raw fiber image and flat file from user
[fnameInput,pnameI] = uigetfile({[pnameStartI,'\*.tif;*.tiff']},...
    'Select the INPUT file...');
imFiberImagePath = [pnameI,fnameInput];
[fnameFlat,pnameF] = uigetfile({[pnameStartF,'\*.tif;*.tiff']},...
    'Select the FLAT file...');

pnameStartI = pnameI;
pnameStartF = pnameF;

% Get background file, if it exists
imFiberFlatPath = [pnameF,fnameFlat];
[fnameBG,pnameBG,nCancel] = uigetfile({[pnameStartBG,...
    '\*.tif;*.tiff']},'Select the BACKGROUND file...');
if ~nCancel
    imFiberBGPath = NaN(1);
else
    imFiberBGPath = [pnameBG,fnameBG];
    pnameStartBG = pnameBG;
end

% Save configuration
save('fcp_config.mat','pnameStartI','pnameStartF',...
    'pnameStartBG','-append');

%% Get file parameters
frameTime_ms = getParameter([pnameI,fnameInput],'Frame length (ms):');
inputVoltX = getParameter([pnameI,fnameInput],'Scan Voltage (X):');
inputVoltY = getParameter([pnameI,fnameInput],'Scan Voltage (Y):');
if isnan(inputVoltX)
    inputVoltX = getParameter([pnameI,fnameInput],...
        'Scanning Voltage (V):');
    inputVoltY = inputVoltX;
end
if ~isnan(frameTime_ms)
    tSeries = 1;
else
    tSeries = 0;
end

flatVoltX = 3;
flatVoltY = 3;

%% Open files
infoFlat = imfinfo(imFiberFlatPath);
% Check if flat file is multi-page file
if numel(infoFlat)>1
    error('Flat file should only have one image');
end
% Read in flat file
imFlat = double(imread(imFiberFlatPath,'tif',1));
% (Optional) Resize flat file
imFlat = imresize(imFlat,[1000,1000]);
infoFlat(1).Width = size(imFlat,2);
infoFlat(1).Height = size(imFlat,1);
% Gaussian denoising filter
imFlat = imgaussfilt(imFlat,1);
% Get processing size from flat file
xFlat = infoFlat(1).Width;
yFlat = infoFlat(1).Height;
% Read input images
infoInput = imfinfo(imFiberImagePath);
% Get original size from input file
xInput = infoInput(1).Width;
yInput = infoInput(1).Height;

% Scale up according to zoom level
inputVoltPerPixel = inputVoltX/xInput;
flatVoltPerPixel = flatVoltX/xFlat;
zoomScale = inputVoltPerPixel/flatVoltPerPixel;
xInputSc = round(xInput*zoomScale);
yInputSc = round(yInput*zoomScale);
% Establish cropping values
cropIn.StartX = floor((xFlat-xInputSc)/2)+1;
cropIn.EndX = cropIn.StartX + xInputSc-1;
cropIn.StartY = floor((yFlat-yInputSc)/2)+1;
cropIn.EndY = cropIn.StartY + yInputSc-1;

numImages = numel(infoInput);
imInputPad = double(zeros(yFlat,xFlat,numImages));
imInput = double(zeros(yInputSc,xInputSc,numImages));
% Read images into stack
hWait = waitbar(0, sprintf('Importing %d images...',numImages));
for ii = 1:numImages
    imInputTemp = imread(imFiberImagePath,'tif',ii);
    imInputTemp = imresize(imInputTemp,[yInputSc,xInputSc],'bicubic');
    imInput(:,:,ii) = imInputTemp;
    imInputPad(:,:,ii) = centerPadCrop(imInputTemp,xFlat,yFlat);
    waitbar(ii/numImages,hWait);
end
% imInputPadMax = imgaussfilt(max(imInputPad,[],3),1);
imInputPadGauss = imgaussfilt(mean(imInputPad(:,:,:),3),1);
close(hWait);

% Load background file if it exists
if ~isnan(imFiberBGPath)
    imBG = imread(imFiberBGPath,'tif',1);
    imBG = imresize(imBG,[yInputSc,xInputSc],'bicubic');
    imBG = centerPadCrop(imBG,xFlat,yFlat);
    imBG = imgaussfilt(mean(imBG,3),1);
    imInputPadGauss = imInputPadGauss - imBG*0.8;
    imInputPadGauss(imInputPadGauss<0) = 0;
end

% Run findFiberCores function on input image
SE = strel('disk',3);
% Set threshold for input pixels
cont = 0;
while ~cont
    [inputCentroids,imInputBW] = getFiberCentroids(imInputPadGauss,...
        SE, threshInput);
    figChk = figure; 
    pos = get(figChk,'Position');
    set(figChk,'Position',[0.5*pos(1),0.5*pos(2),pos(3)*2,pos(4)*1]);
    subplot(1,2,1); imshow(imInputPadGauss./max(imInputPadGauss(:)));
    subplot(1,2,2); imshow(imInputBW);
    text(10,10,sprintf('threshInput = %d',threshInput),'Color','w');
    threshUser = input('Enter new threshold value or C to continue [C]: ');
    if ischar(threshUser) || isempty(threshUser)
        cont = 1;
    else
        threshInput = threshUser;
    end
    close(figChk);
end
save('fcp_config.mat','threshInput','-append');

% Run findFiberCores function on flat image
SE = strel('disk',3);
% Set threshold for flat image pixels
cont = 0;
while ~cont
    [flatCentroids,imFlatBW] = getFiberCentroids(imFlat,...
        SE, threshFlat);
    figChk = figure;
    pos = get(figChk,'Position');
    set(figChk,'Position',[0.5*pos(1),0.5*pos(2),pos(3)*2,pos(4)*1]);
    subplot(1,2,1); imshow(imFlat./max(imFlat(:)));
    subplot(1,2,2); imshow(imFlatBW);
    text(10,10,sprintf('threshFlat = %d',threshFlat),'Color','w');
    threshUser = input('Enter new threshold value or C to continue [C]: ');
    if ischar(threshUser) || isempty(threshUser)
        cont = 1;
    else
        threshFlat = threshUser;
    end

    close(figChk);
end
save('fcp_config.mat','threshFlat','-append');

%% Registration
% Select starting position
startX = 0;
startY = 0;
% Run translation registration code
[xTrans, yTrans] = rigidAlignFiber(imInputBW,...
    imFlatBW,startX,startY,100,20);
% + TEST + % Check quality of image registration
% Apply transformation to binary flat image
imFlatBWReg = imtranslate(imFlatBW,[xTrans, yTrans]);
% Check accuracy of alignment
overlap = mean(mean(imInputBW.*imFlatBWReg));
fprintf('\nOverlap number: %.3f\n',overlap);
% Show alignment overlap images
figChk = figure;
imshowpair(imInputBW, imFlatBWReg,'Scaling','joint');
pause;
close(figChk);

%% Generate core structure with registered flat core coordinates
% Translate centroids according to registration
centroids(:,1) = flatCentroids(:,1) + xTrans;
centroids(:,2) = flatCentroids(:,2) + yTrans;

% Crop centroid positions to size of input image
centroids = centroids((centroids(:,1)<=cropIn.EndX &...
    centroids(:,1)>cropIn.StartX) & (centroids(:,2)<=cropIn.EndY...
    & centroids(:,2)>cropIn.StartY),:);
centroids(:,1) = round(centroids(:,1))-cropIn.StartX+1;
centroids(:,2) = round(centroids(:,2))-cropIn.StartY+1;
numCores = size(centroids,1);

%% Corrector matrix creation
% Get values from flat image
sampleRad = 4;
[flatVal] = getCentroidValues(imFlat,flatCentroids,sampleRad);
% % Create flat image with values
imFlatVal = makeFiberImage(flatCentroids,flatVal,...
    yFlat,xFlat,strel('disk',sampleRad));
% Apply transformation to flat image with values
imFlatVal = imtranslate(imFlatVal,[xTrans, yTrans]);
[correctorVal] = getCorrectorValues(imFlatVal, centroids, cropIn);

%% Acquire core values from each image
inputVal = getCentroidValues(imInput,centroids,sampleRad);
% Apply corrector values to each pixel
outputVal = inputVal.*correctorVal;

%% Remove static values
% if (numImages>1)
%     cont = 0;
%     outputValMean = zeros(size(outputVal,1),1);
%     outputValMean = mean(outputVal(:,userBGRange(1):userBGRange(2)),2);
%     imOutputValMean = makeFiberImage(centroids,outputValMean,yInputSc,xInputSc,strel('disk',sampleRad));
%     outputVal = outputVal - outputValMean;
%     outputVal(outputVal<0) = 0;
% end
% meanOutputVal = mean(outputVal,2);
% madMeanOutputVal = median(abs(meanOutputVal-median(meanOutputVal)));
% 
% figure(8); clf; hold on;
% plot(meanOutputVal,'.b');
% meanOutputVal(meanOutputVal<(1*mean(meanOutputVal)+5*madMeanOutputVal))=0;
% plot(meanOutputVal,'ok');
% fixedOutputVal = outputVal - meanOutputVal;
% fixedOutputVal(fixedOutputVal<0) = 0;
% % figure(9);
% % imagesc(makeFiberImage(flatCentroids,fixedOutputVal(:,9),yFlat,xFlat,strel('disk',sampleRad)));
% outputVal = fixedOutputVal;
%% Filter cores independently
outputValFilt = zeros(size(outputVal));
% If not a time series, filter by simple 3D Gaussian window
if (numImages>1) && ~(tSeries)
    gWin = gausswin(3);
    for ii = 1:numCores
        outputValFilt(ii,:) = nanconv(outputVal(ii,:),...
            gWin,'edge','1d');
    end
% If time series, filter with Savitzky-Golay temporal filter
% based on fiber-core sensitivity weight
elseif (numImages>1) && (tSeries)
    frameLength = 9;
    sgOrder = 3;
    outputValFilt = zeros(size(outputVal));
    correctorValHigh = median(correctorVal)+4*std(correctorVal);
    for ii = 1:numCores
        filtWeight = correctorVal(ii,1)/correctorValHigh;
        filtWeight(filtWeight>1) = 1;
        tempValFilt = sgolayfilt(outputVal(ii,:),...
            sgOrder,frameLength);
        outputValFilt(ii,:) = filtWeight*tempValFilt +...
            (1-filtWeight)*outputVal(ii,:);
    end
else
    outputValFilt = outputVal;
end

%% Equalize means of output and input images
inputVal = inputVal./mean(inputVal(:));
outputVal = outputVal./mean(outputVal(:));
outputValFilt = outputValFilt./mean(outputValFilt(:));

%% Generate output images using core values
seIm = strel('disk',4);
imInputVal = makeFiberImage(centroids,inputVal,...
    yInputSc,xInputSc,seIm);
imOutputValFilt = makeFiberImage(centroids,outputValFilt,...
    yInputSc,xInputSc,seIm);
% Rescale core images
imInputVal = imInputVal./max(imInputVal(:));
imOutputValFilt = imOutputValFilt./max(imOutputValFilt(:));

%% Grid and write data
% Grid data with nearest interp
subfactor = 1;
method = 'nearest';
imOutputNearest = gridFiberCores(centroids,outputValFilt,...
    yInputSc,xInputSc,subfactor,method);
imOutputNearest = imOutputNearest./max(imOutputNearest(:));
subfactor = 2;
method = 'natural';
imOutputNatural2x = gridFiberCores(centroids,outputValFilt,...
    yInputSc,xInputSc,subfactor,method);
imOutputNatural2x = imOutputNatural2x./max(imOutputNatural2x(:));

chanNum = str2double(fnameInput(strfind(lower(fnameInput),'.tif')-1));
if (chanNum~=1) && (chanNum~=2)
    chanNum=0;
end
fileNum = sscanf(fnameInput,'%d_',1);
saveFileName = sprintf('%d_CH%d',fileNum,chanNum);

saveFolderName = [pnameI,'ProcessedSimple\',sprintf('%d_Output',fileNum)];
status = mkdir(saveFolderName);
bitDepth = 16;

writeFiberImages(saveFolderName,saveFileName,bitDepth,...
        imInputVal,imOutputValFilt,imOutputNearest,imOutputNatural2x);