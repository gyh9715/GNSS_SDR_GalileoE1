function  Nco=loopFilter2(T,a3,b3,BL,Error,oldError_1,oldError_2,oldNco_1,oldNco_2)
% second-order filter for loop, third-order system for DLL
%   Inputs:
%       T               - coherent integration period .
%       a3,b3           - filter parameters
%       BL              - code error which is the output of the code
%                         discriminator.
%       oldError_1      - carrier error of the pervious period
%       oldError_2      - carrier error of 2 periods before
%       oldNco_1        - carrier NCO output of the pervious period
%       oldNco_2        - carrier NCO output of 2 periods before
%   Outputs:
%       Nco             - code NCO output of the prompt period

% Solve natural frequency
wn=4*(a3*b3-1)/(a3*b3*b3+a3*a3-b3)*BL;

a=wn^3;
b=a3*wn^2;
c=b3*wn;

aa=c+b*T/2+a*T*T/4;
bb=a*T*T/2-2*c;
cc=c-b*T/2+a*T*T/4;

Nco=aa*Error+bb*oldError_1+cc*oldError_2+2*oldNco_1-oldNco_2;