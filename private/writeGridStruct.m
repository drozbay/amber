function writeGridStruct(folderName,baseFilename,bitDepth,stOutGrid)
%WRITEGRIDSTRUCT

fileStartStr = 'imGridOut';
normalizeImages = 1;

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

ff = 0;
if isfield(stOutGrid,'yesFilt')
    filtType = 'Filt';
    if isfield(stOutGrid.yesFilt,'resHi')
        resType = 'HiRes';
        if isfield(stOutGrid.yesFilt.resHi,'natInterp')
            interpType = 'Nat';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resHi.natInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.yesFilt.resHi,'linInterp')
            interpType = 'Lin';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resHi.linInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.yesFilt.resHi,'neaInterp')
            interpType = 'Nea';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resHi.neaInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
    end
    if isfield(stOutGrid.yesFilt,'resMed')
        resType = 'MedRes';
        if isfield(stOutGrid.yesFilt.resMed,'natInterp')
            interpType = 'Nat';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resMed.natInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.yesFilt.resMed,'linInterp')
            interpType = 'Lin';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resMed.linInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.yesFilt.resMed,'neaInterp')
            interpType = 'Nea';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resMed.neaInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
    end
    if isfield(stOutGrid.yesFilt,'resNorm')
        resType = 'NormRes';
        if isfield(stOutGrid.yesFilt.resNorm,'natInterp')
            interpType = 'Nat';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resNorm.natInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.yesFilt.resNorm,'linInterp')
            interpType = 'Lin';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resNorm.linInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.yesFilt.resNorm,'neaInterp')
            interpType = 'Nea';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.yesFilt.resNorm.neaInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
    end
end
if isfield(stOutGrid,'nonFilt')
    filtType = '';
    if isfield(stOutGrid.nonFilt,'resHi')
        resType = 'HiRes';
        if isfield(stOutGrid.nonFilt.resHi,'natInterp')
            interpType = 'Nat';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resHi.natInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.nonFilt.resHi,'linInterp')
            interpType = 'Lin';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resHi.linInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.nonFilt.resHi,'neaInterp')
            interpType = 'Nea';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resHi.neaInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
    end
    if isfield(stOutGrid.nonFilt,'resMed')
        resType = 'MedRes';
        if isfield(stOutGrid.nonFilt.resMed,'natInterp')
            interpType = 'Nat';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resMed.natInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.nonFilt.resMed,'linInterp')
            interpType = 'Lin';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resMed.linInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.nonFilt.resMed,'neaInterp')
            interpType = 'Nea';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resMed.neaInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
    end
    if isfield(stOutGrid.nonFilt,'resNorm')
        resType = 'NormRes';
        if isfield(stOutGrid.nonFilt.resNorm,'natInterp')
            interpType = 'Nat';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resNorm.natInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.nonFilt.resNorm,'linInterp')
            interpType = 'Lin';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resNorm.linInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
        if isfield(stOutGrid.nonFilt.resNorm,'neaInterp')
            interpType = 'Nea';
            ff = ff + 1;
            imageCellArr{ff} = stOutGrid.nonFilt.resNorm.neaInterp;
            filenameCellArr{ff} = [fileStartStr,'_',filtType,resType,interpType];
        end
    end
end


if normalizeImages
    for ii = 1:ff
        imageCellArr{ii} = imageCellArr{ii}./max(imageCellArr{ii}(:));
    end
end
        
        
switch bitDepth
    case 8       
        for ii = 1:ff
            imageToWrite = imageCellArr{ii};
            filenameToWrite = [baseFilename,'_',filenameCellArr{ii}];
            if max(imageToWrite(:))>1
                warning('Max value of %s is %d and will be clipped',filenameCellArr{ii},max(imageToWrite(:)));
            end
            numImages = size(imageToWrite,3);
            imwrite(uint8(imageToWrite(:,:,1)*255),[saveFolderName,'\',filenameToWrite,'.tif']);
            for jj = 2:numImages
                imwrite(uint8(imageToWrite(:,:,jj)*255),[saveFolderName,'\',filenameToWrite,'.tif'],'writemode','append');
            end
        end
    case 16
        for ii = 1:ff
            imageToWrite = imageCellArr{ii};
            filenameToWrite = [baseFilename,'_',filenameCellArr{ii}];
            if max(imageToWrite(:))>1
                warning('Max value of %s is %d and will be clipped',filenameCellArr{ii},max(imageToWrite(:)));
            end
            numImages = size(imageToWrite,3);
            imwrite(uint16(imageToWrite(:,:,1)*2^16-1),[saveFolderName,'\',filenameToWrite,'.tif']);
            for jj = 2:numImages
                imwrite(uint16(imageToWrite(:,:,jj)*2^16-1),[saveFolderName,'\',filenameToWrite,'.tif'],'writemode','append');
            end
        end
end