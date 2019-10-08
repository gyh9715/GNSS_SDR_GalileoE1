function [decodeFECResult] = decodeFEC(decodeInterResult, settings, PageStart , activeChnList)
%decode interleaving from the data channel
%   Inputs:
%       decodeInterResult  - symbols after the processing of interleaving decoder
%       settings           - Receiver settings.
%       firstPage          - the array contains positions of the first
%                            preamble in each channel. The position is 4ms count 
%                            since start of tracking. Corresponding value will
%                            be set to 0 if no valid preambles were detected in
%                            the channel.
%       activeChnList      - list of channels containing valid preambles
%
%
%   Outputs:
%       decodeFECResult    - bits after the processing of decoding FEC with
%                            Viterbi decoder.

for channelNr = activeChnList
    
    % The position of first preamble in Galileo message.
    FirstPageStartNr = PageStart(channelNr);
    
    % The number of pages in data channel. (One page is 1 second, and the 
    % nominal page of Galileo I/NAV is 2 pages)
    PageNr = idivide(settings.msToProcess/4-FirstPageStartNr,settings.symbolRate ,'floor');
    
    for i =1:PageNr
        symbolRate=settings.symbolRate-10;
        symbols = decodeInterResult(channelNr).data(symbolRate * (i-1)+1:symbolRate * i);
        
        % Convert tracking results to symbol.(240 symbolds)
        symbols (symbols == 1 ) =  -1;
        symbols (symbols == 0 ) =   1;
        
        %Take into account the NOT gate in G2 polynomial (Galileo ICD Figure 13, FEC encoder),figure 13. Convolutional Coding Scheme
        symbols(2:2:end)=-1*symbols(2:2:end);
        %make a offset to the delay of viterbi decoder
        tbl=32;
        Symbols_vitdec=[symbols,ones(1,2*tbl)];
        %Generator polynomials
        trellis = poly2trellis(7,[171 133]);
        %implement Viterbi decoder
        decodeBits = vitdec(Symbols_vitdec,trellis,tbl,'cont','unquant');
        decodeBits = decodeBits(tbl+1:end);
        %Record result of decoding FEC
        decodeFECResult(channelNr).data(symbolRate/2 * (i-1)+1:symbolRate/2 * i)=decodeBits;
    end
end
end
