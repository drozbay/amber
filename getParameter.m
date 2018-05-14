function [ paramOut ] = getParameter( filePath, paramLabel )
% getParameter - B. Ozbay (11/15/2017)
% [paramOut] = getParameter(filePath, paramLabel)
% Parses a parameter text file associated with image at filePath for
% a parameter identified by paramLabel and returns the number immediately
% following the label
%
% INPUTS:
% filePath - String path of image file associated with parameter file.
% Parameter file should end with '_parameter.txt'.
% paramLabel - String inside text file to find
%
% OUTPUTS:
% paramOut - Returns a number of type DOUBLE if paramLabel is found,
% otherwise returns NaN

[paramPath,paramName,~] = fileparts(filePath);
cont = 0;
while ~cont
    % Check if parameter file associated with filePath exists.
    % If not, decrement the length of the filePath string and check again.
    if exist([paramPath,'\',paramName,'parameter.txt'], 'file') == 2
        cont = 1;
        fileID = fopen([paramPath,'\',paramName,'parameter.txt']);
        paramText = fileread([paramPath,'\',paramName,'parameter.txt']);
        paramText(strfind(paramText,paramLabel)+length(paramLabel):end);
        if ~isempty(paramText(strfind(paramText,paramLabel)+length(paramLabel):end))
            C = textscan(paramText(strfind(paramText,paramLabel)+length(paramLabel):end),'%f');
            paramOut = C{1};
        else
            paramOut = NaN(1);
        end
        fclose(fileID);
    elseif length(paramName) < 2
        cont = 1;
        paramOut = NaN(1);
    else
        paramName = paramName(1:end-1);
    end
end

