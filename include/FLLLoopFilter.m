function  Nco=FLLLoopFilter(T,a2,BL,Error,oldError_1,oldError_2,oldNco_1,oldNco_2)
%first-order filter for loop, second-order system for DLL

wn=4*a2/(1+a2*a2)*BL;

a=wn^2;
b=a2*wn;
c=0;

aa=c+b*T/2+a*T*T/4;
bb=a*T*T/2-2*c;
cc=c-b*T/2+a*T*T/4;

Nco=aa*Error+bb*oldError_1+cc*oldError_2+2*oldNco_1-oldNco_2;