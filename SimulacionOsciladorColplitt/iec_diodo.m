function I = iec_diodo(v)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Is = 2.52e-0009;
%Nr = 1.752; %
Is = 2.39e-14; %reverse saturatión current
Nr = 1.004; % Reverse current emission coefficient
Vt = 26e-03;
Gmin = 10e-12;
Br=7.946;% Ideal maximum reverse beta

if v >-5*Vt
    Iec = Is*(exp(v/(Nr*Vt))-1)+v*Gmin;
else 
    Iec = -Is+v*Gmin;
    %Iec = Is*(exp(vd/(Nr*Vt))-1)+vd*Gmin;
end 

I = Iec;
end