function [ paramOut ] = getParameter( filePath, paramLabel )
%GETPARAMETER
[paramPath,paramName,~] = fileparts(filePath);
cont = 0;
while ~cont
    [paramPath,'\',paramName,'parameter.txt']
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

