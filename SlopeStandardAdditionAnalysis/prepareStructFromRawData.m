function [ data ] = prepareStructFromRawData( rawData, concVect, tptotalVal, tpVal, twVect, technique )
% This funtion prepares raw data from 1 ms measurement for further
% processing by standardAddition()
% rawData - is data from measurement (nrofpoints x nrofmeasurements)
% concVect - is vector containing the concentration of analyte of each
%            column of rawData
% tptotal - is total number of samples per measurement point (usually
%           tp+tw)
% tpVal - value of tp to be used for callibration (usually low 1-5)
% twVect - vector of values of tw's to prepare the final calibration set
%          (eg. [ 5 10 15 20] - provides 4 sets for calibration)
% technique - voltammetric technique: 'sc' | 'np' | 'dp' | 'sqw'

    if ( numel(concVect) ~= size(rawData,2) ) 
        error('concVect has to describe every column of rawData');
    end

    pos=1;
    for col=1:size(rawData,2)
        for twnum=1:numel(twVect);
            data.Y(:,pos) = setTpTw(rawData(:,col),tptotalVal, tpVal, twVect(twnum), technique);
            data.CONC(1,pos) = concVect(col);
            data.SENS(1,pos) = twnum;
            pos=pos+1;
        end
    end
end

