function valFilt = filterValSeries( valIn, tSeries, varargin )
% filterValSeries - B. Ozbay (04/26/2017)
% [valFilt] = filterValSeries( valIn, tSeries, frameLength, sgOrder, gaussWindow )
% Performs a filter on each individual core value of a core values series.
% Specifying tSeries indicates whether the filter should use a temporal
% Savitsky-Golay filter to preserve some transient features of a time
% series. Otherwise a simple Gaussian filter is used.
%
% INPUTS:
% valIn - NxM matrix of core values for a series with N cores and M frames
% tSeries - Boolean value defining whether or not to interpret the input as
% a time series or not. If tSeries==1, a Savitsky-Golay filter will be
% used. Otherwise, a regular Gaussian window will be used.
% correctorVal (optional) - For time series only, series of correction
% values for each core used to weight amount of filtering
% Optional name-value pairs:
% frameLength - (def:9) For time series only, number of frames to use for
% Savitsky-Golay filter.
% sgOrder - (def:3) For time series only, order of Savitsky-Golay filter
% gaussWindow - (def:3) For non-time series, size of Gaussian filter window
%
% OUTPUTS:
% valFilt - NxM matrix of filtered core values 

% Parse Inputs
p = inputParser;

defaultCorrectorVal = 0;
defaultFrameLength = 9;
defaultSgOrder = 3;
defaultGaussWindow = 3;

addRequired(p, 'valIn', @isnumeric);
addRequired(p, 'tSeries', @isnumeric);
addOptional(p, 'correctorVal', defaultCorrectorVal, @isnumeric)
addParameter(p, 'frameLength', defaultFrameLength, @isnumeric);
addParameter(p, 'sgOrder', defaultSgOrder, @isnumeric);
addParameter(p, 'gaussWindow', defaultGaussWindow, @isnumeric);

parse(p,valIn,tSeries,varargin{:})

correctorVal = p.Results.correctorVal;
frameLength = p.Results.frameLength;
sgOrder = p.Results.sgOrder;
gaussWindow = p.Results.gaussWindow;

if tSeries && ~max(correctorVal(:))
    error('If filtering time series you need correctorVal input');
end

%% Perform filtering
valFilt = zeros(size(valIn));
numCores = size(valIn,1);
numImages = size(valIn,2);
% If not a time series, filter by simple 3D Gaussian window
if (numImages>1) && ~(tSeries)
    gWin = gausswin(gaussWindow);
    for ii = 1:numCores
        valFilt(ii,:) = nanconv(valIn(ii,:),...
            gWin,'edge','1d');
    end
% If time series, filter with Savitzky-Golay temporal filter
% based on fiber-core sensitivity weight
elseif (numImages>1) && (tSeries)
    correctorVal = correctorVal./min(correctorVal(:));
    correctorVal = correctorVal/max(correctorVal(:));
    for ii = 1:numCores
        filtWeight = correctorVal(ii,1);
        tempValFilt = sgolayfilt(valIn(ii,:),...
            sgOrder,frameLength);
        valFilt(ii,:) = filtWeight*tempValFilt +...
            (1-filtWeight)*valIn(ii,:);
    end
else
    valFilt = valIn;
    disp('valIn not an image series. Values are unchanged.');
end

end

