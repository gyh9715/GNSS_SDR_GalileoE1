function [invertSymbol ,activeChnList] = findInvertSymbol(trackResults)
%find the postive and negative of data with the aid of pliot channel data
%   Inputs:
%       trackResults    - output from the tracking function.
%
%   Outputs:
%       invertSymbol    - Invert information for each channel.
%                       - 0 : dont need to invert the data.
%                       - 1 : need to invert the data.

% Secondary codes search can be delayed to a later point in the tracking results
% to avoid noise due to tracking loop transients 
searchStartOffset = 0;

%--- Generate the secondary codes ----------------------------------------
secondaryCodes = [0 0 1 1   1 0 0 0   0 0 0 0   1 0 1 0     1 1 0 1   1 0 0 1    0];

%convert seconday codes from logic level to signal level
secondaryCodes (secondaryCodes >  0)  = -1;
secondaryCodes (secondaryCodes == 0)  =  1;

%--- Make a list of channels excluding not tracking channels --------------
activeChnList = find([trackResults.status] ~= '-');

for channelNr = activeChnList
    
    % Read output from tracking. It contains the secondary codes. The start
    % of record is skiped here to avoid tracking loop transients.
%     bits=trackResults(channelNr).Pilot_I_P(1 + searchStartOffset : end);
    bits=trackResults(channelNr).I_P(1 + searchStartOffset : end);
    
    % Now threshold the output and convert it to -1 and +1 
    bits(bits > 0)  =  1;
    bits(bits <= 0) = -1;
    
    % Correlate tracking output with the secondary codes
    tlmXcorrResult = xcorr(bits, secondaryCodes);
    
    %% Find all starting points off all preamble like patterns ================
    clear index
    clear index2

    xcorrLength = (length(tlmXcorrResult) +  1) /2;

    %--- Find at what index/ms the preambles start ------------------------
    index = find(...
        abs(tlmXcorrResult(xcorrLength : xcorrLength * 2 - 1)) > 20)' + ...
        searchStartOffset;
    if isempty(index)
        activeChnList = setdiff(activeChnList, channelNr);
        disp('Secondary codes not found.')
%         return
    else
    % define a flag the record the signal level of secondary code
    flag = sum(secondaryCodes == bits(index(20):index(20)+24));
    % notice that pilot channel and data channel has opposite signal level
    % in transmitting terminal.
    if flag == 0
        invertSymbol(channelNr)=0;
    else
        invertSymbol(channelNr)=1;
    end
    end   
end
end