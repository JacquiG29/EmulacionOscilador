function L = inductancia(frecuencia,Ceq)
n = frecuencia

switch n
    case 261
        L=407.06e-3;
    case 440
        L=143.92e-3;
    case 1000
        L=27.860e-3;
    case 5000
        L=1.11e-3;
    case 10000
        L=278.63e-6;
    case 15000
        L=0.00012384;
    otherwise
        L = (1/(2*pi*frecuencia))^2*(1/Ceq);
end
end

