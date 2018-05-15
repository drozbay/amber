%% AMBER - Artefactual Multiphoton Bundle Effect Removal
% Author: Baris N. Ozbay, University of Colorado Denver
% Version: 1.0
% Date of latest release: 5/14/2018
%
% This code is an example of the use of the included functions to remove
% the heterogeneity artifact caused by performing multiphoton imaging
% through a non-regular coherent fiber-bundle.
%
% !! NOTE !! %
% ---------- %
% Notes are included in areas of this code that will need to be altered to
% function with the user's image files and system.
% ---------- %
%%
clear all;
close all;

% !! NOTE !! %
% ---------- %
% The writeFiberImages function uses the bioformats toolbox for matlab. The
% following path will need to be updated with the location of the
% bioformats_package.jar if the user wants to use this function.
% ---------- %
javaaddpath('bioformats_package.jar');

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
    startX = 0;
    startY = 0;
    save('fcp_config.mat','pnameInputStart','pnameFlatStart',...
        'threshInput','threshFlat','startX','startY');
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

% Save configuration
save('fcp_config.mat','pnameInputStart','pnameFlatStart','-append');

%% Get file parameters
% Get imaging dimensions for input file
% !! NOTE !! %
% ---------- %
% This code uses the getParameter function, which searches a text file that
% starts with the same name as the *.tif file containing the original
% images for specific parameters. The method of acquiring these parameters
% may need to be altered to extract metadata from specific microscopy
% imaging formats.
% ---------- %
% (Galvanometer scanning voltage is used as the input here. Volts may be
% replaced with any lateral unit, such as microns, but must be consistent)
frameTime_ms = getParameter([pnameInput,fnameInput],'Frame length (ms):');
frameStepSize = getParameter([pnameInput,fnameInput],'Step size (um):');
inputDimX = getParameter([pnameInput,fnameInput],'Scan Voltage (X):');
inputDimY = getParameter([pnameInput,fnameInput],'Scan Voltage (Y):');

% If there is no frame time parameter, assume stack is not a time series
if ~isnan(frameTime_ms)
    tSeries = 1;
else
    tSeries = 0;
end

% !! NOTE !! %
% ---------- %
% Set microns per dimensional units (In this example: microns per volt on
% galvanometers and microns per volt on electrowetting lens). Use a
% calibrated scaling accounting for the magnification of optics after the
% fiber bundle.
umPerDim = 50;
umPerZDim = 5;
inputProperties = struct('frameTime',[],'stepSize',[],'umWidth',[],...
    'umHeight',[]);
inputProperties.frameTime = frameTime_ms;
inputProperties.stepSize = abs(umPerZDim*frameStepSize);
inputProperties.umWidth = umPerDim*inputDimX;
inputProperties.umHeight = umPerDim*inputDimY;

% Get imaging dimensions of flat image
flatDimX = getParameter([pnameFlat,fnameFlat],'Scan Voltage (X):');
flatDimY = getParameter([pnameFlat,fnameFlat],'Scan Voltage (Y):');

%% Import images from files
infoFlat = imfinfo(imFiberFlatPath);
% Check if flat file is multi-page file
if numel(infoFlat)>1
    error('Flat file should only have one image');
end
% Read in flat file
imFlat = double(imread(imFiberFlatPath,'tif',1));
% Resize flat file to standard dimensions
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
inputDimPerPixel = inputDimX/xInput;
flatDimPerPixel = flatDimX/xFlat;
zoomScale = inputDimPerPixel/flatDimPerPixel;
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

%% Get centroid locations from input image
SE = strel('disk',3);
% Set threshold for input pixels based on user input
cont = 0;
disp('Set threshold level until at least 50% of fiber cores are identified');
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
disp('Set threshold level until almost all cores are identified with no noise');
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

%% Rigid registration of flat image onto input image
% Select starting position
startX = 0;
startY = 0;
% Get the range around the starting position to search depending on the FOV
% of the input image relative to the flat image. (May need to set this
% value manually if search extent is not large enough)
regDistance = round(75*(flatDimX/inputDimX+flatDimY/inputDimY)/2);
[xTrans, yTrans] = rigidAlignFiber(imInputBW,...
    imFlatBW,startX,startY,regDistance,20);
% + TEST + % Check quality of image registration
% An obvious convergence on a peak value of correlation should be visible.
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

%% Further processing
% Filter cores independently
outputValFilt = filterValSeries(outputVal,tSeries,correctorVal);

% Equalize means of output and input images
inputVal = inputVal./mean(inputVal(:));
outputVal = outputVal./mean(outputVal(:));
outputValFilt = outputValFilt./mean(outputValFilt(:));

%% Generate output image stacks
% (Optional) Create cleaned representations of original fiber images using
% structuring element to represent the cores.
seIm = strel('disk',4);
imInputVal = makeFiberImage(centroidsInput,inputVal,...
    yInputSc,xInputSc,seIm);
imOutputVal = makeFiberImage(centroidsInput,outputValFilt,...
    yInputSc,xInputSc,seIm);
% Normalize image values
imInputVal = imInputVal./max(imInputVal(:));
imOutputVal = imOutputVal./max(imOutputVal(:));

% Create images of fiber cores as uniform grid using interpolation
% !! NOTE !! %
% ---------- %
% This code uses the gridFiberCores function with two example settings to
% generate images based on nearest neighbor interpolation and natural
% nearest neighbor interpolation. See the gridFiberCores function for
% further details and options.
% ---------- %
subfactor = 1; % Set sub-sampling factor to 1
method = 'nearest'; % Nearest neighbor interpolation
imOutputNearest1x = gridFiberCores(centroidsInput,outputValFilt,...
    yInputSc,xInputSc,subfactor,method);
imOutputNearest1x = imOutputNearest1x./max(imOutputNearest1x(:));
    
subfactor = 2; % Set sub-sampling factor to 2
method = 'natural'; % Natural nearest-neighbor interpolation
imOutputNatural2x = gridFiberCores(centroidsInput,outputValFilt,...
    yInputSc,xInputSc,subfactor,method);
imOutputNatural2x = imOutputNatural2x./max(imOutputNatural2x(:));

%% Write output files
[~,saveFilePrefix,~] = fileparts(fnameInput);
saveFolderName = [pnameInput,'Processed'];
status = mkdir(saveFolderName);
bitDepth = 16;
writeFiberImages(saveFolderName,saveFilePrefix,bitDepth,inputProperties,...
    imInputVal,imOutputVal,imOutputNearest1x,imOutputNatural2x);
