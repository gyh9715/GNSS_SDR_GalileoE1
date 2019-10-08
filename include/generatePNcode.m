function [PNcode] = generatePNcode(PRN,dataType,codeType)
% generate pn code for pilot and data channel,
% generate BPSK-PN code for ASPeCT tracking method.
%   Inputs:
%       PRN             - PRN number of the sequence.
%       dataType        - recoginze data and pilot component
%                       data component PN code :   1
%                       pilot component PN code :  2
%       codeType        - generate PN code for ASPeCT tracking method
%                       E1 PN code with sub-carrier:    1
%                       E1 PN code without sub-carrier: 2
%   Outputs:
%       PNcode         - local ranging code according to the requirements
%                      from the input.


% generate 4092*2 bits E1 ranging code with sub-carrier
if codeType == 1
    if dataType == 2 %data component
        BOC11SubCarrier = [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1];
        BOC61SubCarrier = [-1,1,-1,1,-1,1,-1,1,-1,1,-1,1];
    else             %pilot component
        BOC11SubCarrier = [1,1,1,1,1,1,-1,-1,-1,-1,-1,-1];
        BOC61SubCarrier = [1,-1,1,-1,1,-1,1,-1,1,-1,1,-1];
    end
    CBOCSubCarrier = sqrt(10/11)*BOC11SubCarrier + sqrt(1/11)*BOC61SubCarrier;
    RangingCode = generateE1RangingCode(PRN,dataType);
    for i=1:4092
        PNcode((i-1)*12+1:i*12) = RangingCode(i)*CBOCSubCarrier;
    end

% generate 4092*2 bits E1 ranging code without sub-carrier
else
    SubCarrier=[1,1,1,1,1,1,1,1,1,1,1,1,];
    RangingCode = generateE1RangingCode(PRN,dataType);
    for i=1:4092
        PNcode((i-1)*2+1:i*2) = RangingCode(i)*SubCarrier;
    end
end