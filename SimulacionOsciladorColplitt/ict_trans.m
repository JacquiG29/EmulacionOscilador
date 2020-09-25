function Ict = ict_trans(Iec,Icc,v1)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
Vaf = 63.2;
Ict = (Icc - Iec)*(1/(1+v1/Vaf));

end