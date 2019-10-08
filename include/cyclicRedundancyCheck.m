function [CRCresult] = cyclicRedundancyCheck(decodeFECResult, activeChnList)
%decode interleaving from the data channel
%   Inputs:
%       decodeFECResult    - bits after the processing of decoding FEC with
%                            Viterbi decoder.
%       activeChnList      - list of channels containing valid preambles
%
%   Outputs:
%       CRCresult(channelNr).result(i)  - record if the CRC is successful
%                                         for channelNr in word i
%       CRCresult(channelNr).Word(i)  - record 128 bits data to compose 
%                                       every single word in Galileo message. 


for channelNr = activeChnList
    % Calculate the total nominal page number need to be processed.
    PageNr = length(decodeFECResult(channelNr).data)/240;
    % adjust the start point according to page type (odd or even).
    if decodeFECResult(channelNr).data(1) == 0
        startPoint =1;
    else
        startPoint =121;
    end
    % take CRC for process every nominal page
    for i =1:PageNr
        %extract data bits of one word in Galileo meassage.
        subFrameBits = decodeFECResult(channelNr).data(startPoint+(i-1)*240:startPoint-1+i*240);
        subFrameBitsCRC = [subFrameBits(1:120-6),subFrameBits(121:240-14)];
        % CRC Generator polynomials
        Polynomials=zeros(1,25);
        indexx=[0,1, 3, 4, 5, 6, 7,10, 11, 14, 17, 18, 23, 24]+1;
        Polynomials(1,indexx)=1;
        Polynomials= flip(Polynomials,2);
        % implement CRC.
        R=length(Polynomials)-1; 
        [q,r] = deconv(subFrameBitsCRC,Polynomials);
        r2=mod(r(end-R+1:end),2);
        %Record the status and the data bits for every word if CRC is passed.
        if sum(r2==0)==length(r2)
            CRCresult(channelNr).result(i) = 0;
            WordContent =[subFrameBits(3:120-6),subFrameBits(123:138)];
            CRCresult(channelNr).Word(1+128*(i-1):128*i)=WordContent;
        else
        %Record the status if CRC is not passed.
            CRCresult(channelNr).result(i) = 1;
        end
    end
end    
end
