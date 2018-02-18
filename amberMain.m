%% AMBER - Artefactual Multiphoton Bundle Effect Removal
% Algorithm for removing the multiphoton induced artifact when imaging
% through a coherent fiber-bundle
clear all;
close all;

javaaddpath('C:\Dropbox\Miscellaneous\Matlab\DropboxMatlabPath\bfmatlab\bioformats_package.jar');

% Open previous configuration settings, if they exist. Otherwise create new
% configuration file with default settings.
if exist('fcp_config.mat', 'file') == 2
    load('fcp_config.mat');
else
    pnameInputStart = 'C:\';
    pnameFlatStart = 'C:\';
    pnameBGStart = 'C:\';
    threshInput = 0.5;
    threshFlat = 0.5;
    threshBG = 0.5;
    startX = 0;
    startY = 0;
    save('fcp_config.mat','pnameInputStart','pnameFlatStart','pnameBGStart',...
        'threshInput','threshFlat','threshBG','startX','startY');
end
% Get input raw fiber image and flat file from user
[fnameInput,pnameInput,nCancel] = uigetfile({[pnameInputStart,'\*.tif;*.tiff']},...
    'Select the INPUT file...');
if ~nCancel
    error('User cancelled');
end
imFiberImagePath = [pnameInput,fnameInput];
[fnameFlat,pnameFlat,nCancel] = uigetfile({[pnameFlatStart,'\*.tif;*.tiff']},...
    'Select the FLAT file...');
if ~nCancel
    error('User cancelled');
end
imFiberFlatPath = [pnameFlat,fnameFlat];
pnameInputStart = pnameInput;
pnameFlatStart = pnameFlat;

% Ask for background file
[fnameBG,pnameBG,nCancel] = uigetfile({[pnameBGStart,...
    '\*.tif;*.tiff']},'Select the BACKGROUND file or press Cancel');
if ~nCancel
    imFiberBGPath = NaN(1);
else
    imFiberBGPath = [pnameBG,fnameBG];
    pnameBGStart = pnameBG;
end

% Save configuration
save('fcp_config.mat','pnameInputStart','pnameFlatStart',...
    'pnameBGStart','-append');

%% Get file parameters
% Get voltage scan range for input file
frameTime_ms = getParameter([pnameInput,fnameInput],'Frame length (ms):');
frameStepSize = getParameter([pnameInput,fnameInput],'Step size (um):');
inputVoltX = getParameter([pnameInput,fnameInput],'Scan Voltage (X):');
inputVoltY = getParameter([pnameInput,fnameInput],'Scan Voltage (Y):');
inputLaserPercent = getParameter([pnameInput,fnameInput],'Laser Percent:');

if isnan(inputVoltX)
    inputVoltX = getParameter([pnameInput,fnameInput],...
        'Scanning Voltage (V):');
    inputVoltY = inputVoltX;
end
if ~isnan(frameTime_ms)
    tSeries = 1;
else
    tSeries = 0;
end

umPerVolt = 70;
umPerEWVolt = 5;
inputProperties = struct('frameTime',[],'stepSize',[],'umWidth',[],...
    'umHeight',[]);
inputProperties.frameTime = frameTime_ms;
inputProperties.stepSize = abs(umPerEWVolt*frameStepSize);
inputProperties.umWidth = umPerVolt*inputVoltX;
inputProperties.umHeight = umPerVolt*inputVoltX;

% Get voltage scan range of flat image
flatVoltX = getParameter([pnameFlat,fnameFlat],'Scan Voltage (X):');
flatVoltY = getParameter([pnameFlat,fnameFlat],'Scan Voltage (Y):');
flatLaserPercent = getParameter([pnameFlat,fnameFlat],'Laser Percent:');
if isnan(flatVoltX)
    flatVoltX = getParameter([pnameInput,fnameInput],...
        'Scanning Voltage (V):');
    flatVoltY = flatVoltX;
end

% Get voltage scan range of background image
if ~isnan(imFiberBGPath)
    bgVoltX = getParameter([pnameBG,fnameBG],'Scan Voltage (X):');
    bgVoltY = getParameter([pnameBG,fnameBG],'Scan Voltage (Y):');
    bgLaserPercent = getParameter([pnameBG,fnameBG],'Laser Percent:');
    if isnan(bgVoltX)
        bgVoltX = getParameter([pnameInput,fnameInput],...
            'Scanning Voltage (V):');
        bgVoltY = bgVoltX;
    end
end

%% Open files
infoFlat = imfinfo(imFiberFlatPath);
% Check if flat file is multi-page file
if numel(infoFlat)>1
    error('Flat file should only have one image');
end
% Read in flat file
imFlat = double(imread(imFiberFlatPath,'tif',1));
% (Optional) Resize flat file
imFlat = imresize(imFlat,[1000,1000],'bicubic');
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
imInputPadGauss = imgaussfilt(mean(imInputPad(:,:,:),3),1);
close(hWait);

%% Run findFiberCores function on input image
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

%% Run findFiberCores function on flat image
SE = strel('disk',3);
% Set threshold for flat image pixels
cont = 0;
while ~cont
    [centroidsFlat,imFlatBW] = getFiberCentroids(imFlat,...
        SE, threshFlat);
    imFlatBW = centerPadCrop(imFlatBW,1000,1000);
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

%% Load background file if it exists
if ~isnan(imFiberBGPath)
    % BG file pixel ratios
    infoBG = imfinfo(imFiberBGPath);
    xBG = infoBG(1).Height;
    yBG = infoBG(1).Height;
    BGVoltPerPixel = bgVoltX/xBG;
    zoomScaleBG = BGVoltPerPixel/flatVoltPerPixel;
    xBGSc = round(xBG*zoomScaleBG);
    yBGSc = round(yBG*zoomScaleBG);
    % Import BG file
    numImagesBG = numel(infoBG);
    imBGPad = double(zeros(yFlat,xFlat,numImagesBG));
    imBG = double(zeros(yBGSc,xBGSc,numImagesBG));
    hWait = waitbar(0, sprintf('Importing %d background images...',numImagesBG));
    for ii = 1:numImagesBG
        imBGTemp = imread(imFiberBGPath,'tif',ii);
        imBGTemp = imresize(imBGTemp,[xBGSc,yBGSc],'bicubic');
        imBG(:,:,ii) = imBGTemp;
        imBGPad(:,:,ii) = centerPadCrop(imBGTemp,xFlat,yFlat);
        waitbar(ii/numImagesBG,hWait);
    end
    close(hWait);
    imBGPadGauss = imgaussfilt(mean(imBGPad(:,:,:),3),1);
    
    % Run getFiberCentroids function on BG image
    SE = strel('disk',3);
    % Set threshold for BG image pixels
    cont = 0;
    while ~cont
       [BGCentroids,imBGBW] = getFiberCentroids(imBGPadGauss,...
            SE, threshBG);
        figChk = figure;
        pos = get(figChk,'Position');
        set(figChk,'Position',[0.5*pos(1),0.5*pos(2),pos(3)*2,pos(4)*1]);
        subplot(1,2,1); imshow(imBGPadGauss./max(imBGPadGauss(:)));
        subplot(1,2,2); imshow(imBGBW);
        text(10,10,sprintf('threshBG = %d',threshBG),'Color','r');
        threshUser = input('Enter new threshold value or C to continue [C]: ');
        if ischar(threshUser) || isempty(threshUser)
            cont = 1;
        else
            threshBG = threshUser;
        end
        close(figChk);
    end
    save('fcp_config.mat','threshBG','-append');
end

%% Registration
% Select starting position
startX = 0;
startY = 0;
% Run translation registration code
regDistance = round(75*(flatVoltX/inputVoltX+flatVoltY/inputVoltY)/2);
[xTrans, yTrans] = rigidAlignFiber(imInputBW,...
    imFlatBW,startX,startY,regDistance,20);
% + TEST + % Check quality of image registration
% Apply transformation to binary flat image
imFlatBWReg = imtranslate(imFlatBW,[xTrans, yTrans]);
% Check accuracy of alignment
overlap = mean(mean(imInputBW.*imFlatBWReg));
fprintf('\nOverlap number: %.3f\n',overlap);
% Show alignment overlap images
figChk = figure;
imshowpair(imInputBW, imFlatBWReg,'Scaling','joint');
getkeywait(4);
% pause;
close(figChk);

%% Generate core structure with registered flat core coordinates
% Set sampling radius
sampleRad = 4;
% Translate centroids according to registration
centroidsInput = [centroidsFlat(:,1) + xTrans , centroidsFlat(:,2) + yTrans];

% Crop centroid positions to size of input image
centroidsInput = centroidsInput((centroidsInput(:,1)<=cropIn.EndX &...
    centroidsInput(:,1)>cropIn.StartX) & (centroidsInput(:,2)<=cropIn.EndY...
    & centroidsInput(:,2)>cropIn.StartY),:);
centroidsInput(:,1) = ceil(centroidsInput(:,1))-cropIn.StartX;
centroidsInput(:,2) = ceil(centroidsInput(:,2))-cropIn.StartY;
numCores = size(centroidsInput,1);
% Acquire core values from each image
inputVal = getCentroidValues(imInput(:,:,:),centroidsInput,sampleRad);
inputVal(13039,:) = 0;

%% Corrector matrix creation
% Get values from flat image
imFlatReg = imtranslate(imFlat,[xTrans,yTrans]);
imFlatCrop = centerPadCrop(imFlatReg,xInputSc,yInputSc);
flatVal = getCentroidValues(imFlatCrop,centroidsInput,sampleRad);
% Create flat image with values
imFlatVal = makeFiberImage(centroidsInput,flatVal,...
    yInputSc,xInputSc,strel('disk',sampleRad));
% Retrieve correction values
[correctorVal] = getCorrectorValues(imFlatVal, centroidsInput, cropIn);
% Apply corrector values to each core
outputVal = inputVal.*correctorVal;

% Filter cores independently
outputValFilt = filterValSeries(outputVal,tSeries,correctorVal);


%% Process files
% Check for BG File
if ~isnan(imFiberBGPath)
    % Perform BG file registration
    % Select starting position
    startX = xTrans;
    startY = yTrans;
    % Run translation registration code
    [xTransBG, yTransBG] = rigidAlignFiber(imFlatBWReg,...
        imBGBW,startX,startY,round(regDistance/2),20);
    % Check quality of image registration
    % Apply transformation to binary flat image
    imBGBWReg = imtranslate(imBGBW,[xTransBG, yTransBG]);
    % Check accuracy of alignment
    overlap = mean(mean(imInputBW.*imFlatBWReg));
    fprintf('\nOverlap number: %.3f\n',overlap);
    % Show alignment overlap images
    figChk = figure;
    imshowpair(imBGBWReg, imFlatBWReg,'Scaling','joint');
    getkeywait(4);
    close(figChk);

    %% Scaling level for background subtraction (multiplicative)
    BGValScale = (inputLaserPercent/bgLaserPercent)^2*0.4;
    % Apply transformation to original background image
    imBGReg = imtranslate(imBGPad,[xTransBG+0, yTransBG+0]);
    imBGRegCrop = centerPadCrop(imBGReg,xInputSc,yInputSc);
    BGVal = getCentroidValues(imBGRegCrop,centroidsInput,sampleRad);
    %
    inputValSubBG = inputVal - mean(BGVal,2)*BGValScale;
    inputValSubBG(inputValSubBG<0)=0;
    % - TEST - %%
    
    imTest = makeFiberImage(centroidsInput,mean(inputValSubBG,2),...
        yInputSc,xInputSc,strel('disk',4));
    figure(4); imagesc(imTest);
    colormap(gca, jet(256)); colorbar; 
    %%
    % Corrector matrix with background subtracted flat
    BGValScaleFlat = (flatLaserPercent/bgLaserPercent)^2;
    flatValSubBG = flatVal-mean(BGVal,2)*BGValScaleFlat;
    % Create flat image with values
    imFlatValSubBG = makeFiberImage(centroidsInput,flatValSubBG,...
        yFlat,xFlat,strel('disk',sampleRad));
    % Retrieve correction values
    [correctorValSubBG] = getCorrectorValues(imFlatValSubBG, centroidsInput, cropIn);
    % Apply corrector values to each core
    outputValSubBG = inputValSubBG.*correctorValSubBG;
    
    % Filter cores independently
    outputValFiltSubBG = filterValSeries(outputValSubBG,tSeries,correctorValSubBG,'frameLength',5);
    
    % Equalize means of output and input images
    inputValSubBG = inputValSubBG./mean(inputValSubBG(:));
    outputValSubBG = outputValSubBG./mean(outputValSubBG(:));
    outputValFiltSubBG = outputValFiltSubBG./mean(outputValFiltSubBG(:));

%     % Generate output images using core values
%     seIm = strel('disk',4);
%     imInputValSubBG = makeFiberImage(centroidsInput,inputValSubBG,...
%         yInputSc,xInputSc,seIm);
%     imOutputValFiltSubBG = makeFiberImage(centroidsInput,outputValFiltSubBG,...
%         yInputSc,xInputSc,seIm);
%     % Rescale core images
%     imInputValSubBG = imInputValSubBG./max(imInputValSubBG(:));
%     imOutputValFiltSubBG = imOutputValFiltSubBG./max(imOutputValFiltSubBG(:));
    
    %% Grid and write data
    % Grid data with nearest interp
%     subfactor = 1;
%     method = 'nearest';
%     imOutputNearestSubBG = gridFiberCores(centroidsInput,outputValFiltSubBG,...
%         yInputSc,xInputSc,subfactor,method);
%     imOutputNearestSubBG = imOutputNearestSubBG./max(imOutputNearestSubBG(:));
%     
%     subfactor = 2;
%     method = 'natural';
%     
%     imOutputNatural2xSubBG = gridFiberCores(centroidsInput,outputValSubBG,...
%         yInputSc,xInputSc,subfactor,method);
%     imOutputNatural2xSubBG = imOutputNatural2xSubBG./max(imOutputNatural2xSubBG(:));
    subfactor = 2;
    method = 'nearest';
    imOutputNearest2xSubBG = gridFiberCores(centroidsInput,outputValSubBG,...
        yInputSc,xInputSc,subfactor,method);
    imOutputNearest2xSubBG = imOutputNearest2xSubBG./max(imOutputNearest2xSubBG(:));
%     
%     imOutputNearest1xSubBG = gridFiberCores(centroidsInput,outputValSubBG,...
%         yInputSc,xInputSc,1,'nearest');
%     imOutputNearest1xSubBG = imOutputNearest1xSubBG./max(imOutputNearest1xSubBG(:));
%     
    
%     imOutputCoreSet = gridFiberCores(centroidsInput,outputValFiltSubBG,...
%         yInputSc,xInputSc,1,'coreset');
    
    imOutputFiltNatural2xSubBG = gridFiberCores2(centroidsInput,outputValFiltSubBG,...
        yInputSc,xInputSc,subfactor,method);
    imOutputFiltNatural2xSubBG = imOutputFiltNatural2xSubBG./max(imOutputFiltNatural2xSubBG(:));
    
    chanNum = str2double(fnameInput(strfind(lower(fnameInput),'.tif')-1));
    if (chanNum~=1) && (chanNum~=2)
        chanNum=0;
    end
    fileNum = sscanf(fnameInput,'%d_',1);
    saveFileName = sprintf('%d_CH%d',fileNum,chanNum);

    saveFolderName = [pnameInput,'ProcessedSimple\',sprintf('%d_Output',fileNum)];
    status = mkdir(saveFolderName);
    bitDepth = 16;
%     writeFiberImages(saveFolderName,saveFileName,bitDepth,inputProperties,...
%             imInputValSubBG,imOutputValFiltSubBG,imOutputNearestSubBG,imOutputNatural2xSubBG);
%     writeFiberImages(saveFolderName,saveFileName,bitDepth,inputProperties,imOutputNatural2xSubBG,imOutputFiltNatural2xSubBG);
    writeFiberImages(saveFolderName,saveFileName,bitDepth,inputProperties,imOutputNearest2xSubBG);


else
 
    %% Equalize means of output and input images
    inputVal = inputVal./mean(inputVal(:));
    outputVal = outputVal./mean(outputVal(:));
    outputValFilt = outputValFilt./mean(outputValFilt(:));

    % %% Generate output images using core values
    % seIm = strel('disk',4);
    % imInputVal = makeFiberImage(centroidsInput,inputVal,...
    %     yInputSc,xInputSc,seIm);
    % imOutputValFilt = makeFiberImage(centroidsInput,outputValFilt,...
    %     yInputSc,xInputSc,seIm);
    % % Rescale core images
    % imInputVal = imInputVal./max(imInputVal(:));
    % imOutputValFilt = imOutputValFilt./max(imOutputValFilt(:));

    %% Grid and write data
    % Grid data with nearest interp
%     subfactor = 1;
%     method = 'nearest';
%     imOutputNearest = gridFiberCores(centroidsInput,outputValFilt,...
%         yInputSc,xInputSc,subfactor,method);
%     imOutputNearest = imOutputNearest./max(imOutputNearest(:));
    
    subfactor = 2;
    method = 'nearest';
    imOutputNearest2x = gridFiberCores(centroidsInput,outputVal,...
        yInputSc,xInputSc,subfactor,method);
    imOutputNearest2x = imOutputNearest2x./max(imOutputNearest2x(:));
    
    % Grid data with natural interp
%     subfactor = 2;
%     method = 'natural';
% 
%     imOutputNatural2x = gridFiberCores(centroidsInput,outputVal,...
%         yInputSc,xInputSc,subfactor,method);
%     imOutputNatural2x = imOutputNatural2x./max(imOutputNatural2x(:));
% 
%     imOutputNatural2xFilt = gridFiberCores(centroidsInput,outputValFilt,...
%         yInputSc,xInputSc,subfactor,method);
%     imOutputNatural2xFilt = imOutputNatural2xFilt./max(imOutputNatural2xFilt(:));

    chanNum = str2double(fnameInput(strfind(lower(fnameInput),'.tif')-1));
    if (chanNum~=1) && (chanNum~=2)
        chanNum=0;
    end
    fileNum = sscanf(fnameInput,'%d_',1);
    saveFileName = sprintf('%d_CH%d',fileNum,chanNum);

    saveFolderName = [pnameInput,'ProcessedSimple\',sprintf('%d_Output',fileNum)];
    status = mkdir(saveFolderName);
    bitDepth = 16;
%
    % writeFiberImages(saveFolderName,saveFileName,bitDepth,inputProperties,...
    %         imInputVal,imOutputValFilt,imOutputNearest,imOutputNatural2x);
    writeFiberImages(saveFolderName,saveFileName,bitDepth,inputProperties,imOutputNearest2x);
end