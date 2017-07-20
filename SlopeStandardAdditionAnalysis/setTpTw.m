function [ res, partial1, partial2 ] = setTpTw( sig, realms, tp, tw, type)
% setTpTw is simple function to process RAW data form electrochemical
% analyzer. I can process data from SCV, DPV and NPV techniques.
% sig - signal with raw data (readout from A/D converter)
% realms - total probing time of one step (i.e. in SCV it is tp, in DPV and
% NPV it is 2*tp
% tw - new wait time (tp + tw =< realtp)
% tp - new tp time
% type - type of technique: 'sc', 'dp', 'dpasv', 'np', 'npasv', 'sqw'
    if ( realms < tw+tp )
        return;
    end

    for (i=1:realms:size(sig,1)-realms) 
        tres( ((i+realms-1)/realms) ) = sum(sig(i+tw:i+tw+tp-1))/tp;
    end;

    if (strcmp(type,'sc'))
        res=tres;
        res=partial1;
        res=partial2;
    elseif (strcmp(type,'dp') || strcmp(type,'np') || strcmp(type,'dpasv') || strcmp(type,'npasv') )
        for (i=1:2:length(tres)-1) 
            res(floor(i/2)+1) = tres(i+1)-tres(i);
            partial1(floor(i/2)+1) = tres(i);
            partial2(floor(i/2)+1) = tres(i+1);
        end;
    elseif strcmp(type,'sqw')
        for (i=2:2:length(tres)-1)
            res(floor(i/2)+1) = tres(i)-tres(i+1);
            partial1(floor(i/2)+1) = tres(i);
            partial2(floor(i/2)+1) = tres(i+1);
        end;        
    end
end