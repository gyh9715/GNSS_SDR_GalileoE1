function  Nco=loopFilter1(T,a2,BL,Error,oldError_1,oldNco_1)
% first-order filter for loop, second-order system for DLL
%   Inputs:
%       T               - coherent integration period .
%       a2              - 2*damping ratio
%       BL              - code error which is the output of the code
%                         discriminator.
%       oldError_1      - code error of the pervious period 
%       oldNco_1        - code NCO output of the pervious period 
%   Outputs:
%       Nco             - code NCO output of the prompt period

% Solve natural frequency
wn=4*a2/(1+a2*a2)*BL;

b0=a2*wn+T*wn*wn/2;
b1=-a2*wn+T*wn*wn/2;

Nco=b0*Error+b1*oldError_1+oldNco_1;