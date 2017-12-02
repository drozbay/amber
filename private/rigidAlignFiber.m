function [xTrans, yTrans] = rigidAlignFiber(imInputBW,imFlatBW, startX, startY,roughDist,fineDist)

% B Ozbay 4/24/2017
% Inputs:
% imCoreBW - Image of sample cores with uniform intensity
% imFlatBW - Image of flat fiber image with uniform intensity
% startX - Initial position guess of X translation
% startY - Initial position guess of Y translation
% roughDist - Range in pixels to perform rough translation search
% fineDist - Range in pixels to perform fine translation registration
% Outputs:
% xTrans - x translation result
% yTrans - y translation result
% 

% Reduce size to speed function
% and convert values to double between 0 and 1
imRed = 0.5;
fixed = imresize(double(imInputBW),imRed);
fixed = fixed-min(min(fixed));
fixed = fixed/(max(max(fixed)));
fixed(fixed<0.7) = 0;
moving = imresize(double(imFlatBW),imRed);
moving = moving-min(min(moving));
moving = moving/(max(max(moving)));
moving(moving<0.7) = 0;

%% Set up rough estimate settings
% X translation settings
xMid = round(startX*imRed);
xDistance = roughDist;%floor(size(moving,2)/10);
xSpacing = 2;%floor(xDistance/15);

% Y translation settings
yMid = round(startY*imRed);
yDistance = roughDist;%floor(size(moving,1)/10);
ySpacing = 2;%floor(yDistance/15);

% Set up overlap matrix
xLength = floor(xDistance/xSpacing);
xx = (0:xLength-1)*xSpacing-floor(xDistance/2)+xMid;
yLength = floor(yDistance/ySpacing);
yy = (0:yLength-1)*ySpacing-floor(yDistance/2)+yMid;
overlapRough = zeros(yLength, xLength);

%% Perform rough estimation
hWait = waitbar(0, sprintf('Running rough estimation with xDistance = %d pixels...',xDistance/imRed));
% figure(10);
for ii = 1:xLength
    for jj = 1:yLength
        moved = imtranslate(moving,[xx(ii) yy(jj)]);
        multTemp = moved.*fixed;
        overlapRough(ii,jj) = mean(mean(multTemp));
%       imshow(multTemp,'InitialMagnification',100);
    end
    waitbar(ii/xLength,hWait);
end
close(hWait);
% Extract best translation
maxRough = max(max(overlapRough));
[ii, jj] = find(overlapRough==maxRough);
xTransRough = round(xx(ii)/imRed);
yTransRough = round(yy(jj)/imRed);
figRigidAlign = figure;
subplot(1,2,1); imagesc(xx/imRed,yy/imRed,overlapRough');
xlim([-1,1]*xDistance/imRed+startX); ylim([-1,1]*xDistance/imRed+startY);


%% Set up precision estimate settings

% Restore size to original for precise evaluation
fixed = double(imInputBW);
fixed = fixed-min(min(fixed));
fixed = fixed/(max(max(fixed)));
moving = double(imFlatBW);
moving = moving-min(min(moving));
moving = moving/(max(max(moving)));

% X translation settings
xMid = round(xTransRough);
xDistance = fineDist;
xSpacing = 1;

% Y translation settings
yMid = round(yTransRough);
yDistance = fineDist;
ySpacing = 1;

% Set up overlap matrix
xLength = floor(xDistance/xSpacing);
xx = (0:xLength-1)*xSpacing-floor(xDistance/2)+xMid;
yLength = floor(yDistance/ySpacing);
yy = (0:yLength-1)*ySpacing-floor(yDistance/2)+yMid;
overlapPrecise = zeros(yLength, xLength);

%% Perform precise estimation
hWait = waitbar(0, sprintf('Running precise estimation with xDistance %d...',xDistance));
for ii = 1:xLength
    for jj = 1:yLength
        moved = imtranslate(moving,[xx(ii) yy(jj)]);
        multTemp = moved.*fixed;
        overlapPrecise(ii,jj) = mean(mean(multTemp));
%         figure(9);
%         imshow(multTemp);
    end
    waitbar(ii/xLength,hWait);
end
close(hWait);
% Extract best translation
maxPrecise = max(max(overlapPrecise));
[ii, jj] = find(overlapPrecise==maxPrecise);
xTrans = xx(ii);
yTrans = yy(jj);
subplot(1,2,2); imagesc(xx,yy,overlapPrecise');
xlim([-1,1]*xDistance+xTransRough); ylim([-1,1]*yDistance+yTransRough);

