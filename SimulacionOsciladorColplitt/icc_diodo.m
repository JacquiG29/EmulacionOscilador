function I = icc_diodo(v)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
%Is = 2.52e-0009;
%Nf = 1.752; %
Vin=1;
Is = 2.39e-14; %reverse saturatión current
Nf = 1.008; % Forward current emission coefficient
Vt = 26e-03;
Gmin = 10e-12;
Bf=294.3;% Ideal maximum forward beta

if v >-5*Vt
    Icc = Is*(exp(v/(Nf*Vt))-1)+v*Gmin;
else
    Icc = -Is+v*Gmin;
    %Icc = Is*(exp(vd/(Nf*Vt))-1)+vd*Gmin;
end
I=Icc;
end