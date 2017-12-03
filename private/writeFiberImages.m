function writeFiberImages(folderName,baseFilename,bitDepth,varargin)
%WRITEFIBERIMAGES
% Writes images to files at specified bitdepth with suffix based on input
% file type (ie. CH1 or CH2)
% INPUT:
% folderName - String for output folder path
% baseFileName - String for base file name
% bitDepth - Integer that specifies either 8 or 16 bit output image
% varargin - List of images to write to files

if bitDepth~=8 && bitDepth~=16
    warning('bitDepth is not 8 or 16, defaulting to 8');
    bitDepth = 8;
end

if strcmp(folderName(end),'\')
    saveFolderName = [folderName];
else
    saveFolderName = [folderName,'\'];
end
mkdir(saveFolderName);

numFiles = nargin-3;

switch bitDepth
    case 8
        for ii = 1:numFiles
            imageToWrite = varargin{ii};
            variableName = inputname(ii+3);
            if max(imageToWrite(:))>1
                warning('Max value of %s is %d and will be clipped',variableName,max(imageToWrite(:)));
            end
            filenameToWrite = [baseFilename,'_',variableName];
            numImages = size(imageToWrite,3);
            imwrite(uint8(imageToWrite(:,:,1)*255),[saveFolderName,'\',filenameToWrite,'.tif']);
            for jj = 2:numImages
                imwrite(uint8(imageToWrite(:,:,jj)*255),[saveFolderName,'\',filenameToWrite,'.tif'],'writemode','append');
            end
        end
    case 16
        for ii = 1:numFiles
            imageToWrite = varargin{ii};
            variableName = inputname(ii+3);
            if max(imageToWrite(:))>1
                warning('Max value of %s is %d and will be clipped',variableName,max(imageToWrite(:)));
            end
            filenameToWrite = [baseFilename,'_',variableName];
            numImages = size(imageToWrite,3);
            imwrite(uint16(imageToWrite(:,:,1)*2^16-1),[saveFolderName,'\',filenameToWrite,'.tif']);
            for jj = 2:numImages
                imwrite(uint16(imageToWrite(:,:,jj)*2^16-1),[saveFolderName,'\',filenameToWrite,'.tif'],'writemode','append');
            end
        end
end