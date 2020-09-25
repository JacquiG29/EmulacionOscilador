clear
close all
%% Import Audio
factor = 12;
fs=44100
%% Parametros simulación
t0 = 0; % tiempo inicial
tf1 = 1; % tiempo final
dt = 1/(factor*fs); % período de muestreo
K = (tf1-t0)/dt; % número de iteraciones
K = round(K);
t = t0:dt:tf1;

%% Parámetros del circuito
Vin=1;
frecuencia=261;%Opciones de frecuencia: 261,440,1000,5000,10000,15000
system=1;%1 nivel 1, 2 nivel 2, 3 nivel 3 simulink
size_system=3; %3 para espacio de estados de 3x3 y 5 para espacios de estado de 5x5
%Espacio de estados de 5x5 solo aplica para system=3

R1 = 100;
R2 = 1*10^(3);
C1 = 10*10^(-6);
C2 = 1*10^(-6);
Ceq=(C1*C2)/(C1+C2);
L = inductancia(frecuencia,Ceq)
fcirc=1/(2*pi*sqrt(L*Ceq))

Rb=1;
Re=0.4683;
Rc=0.85;
Cbc=3.728e-12;
Cbe=1.358e-11 ;

Br=7.946;
Bf=294.3;
%% Matrices del sistema
if system==1
    %obtenido a mano
    A = [ 0       0                                1/C1;...
        0       (1/(C2*R2))*((R1/(R2+R1))-1)     -1/C2;...
        -1/L     1/L                              0];
    
    B=[0;...
        1/(C2*(R2+R1));...
        0];
    
    Q=[-1/C1                 -1/C1     0;...
        R1/(C2*(R2+R1))       0        -R1/(C2*(R2+R1));...
        0                     0        0];
    C=[0,R1/(R2+R1),0];
    
    D=R2/(R2+R1);
    
    R=[(R2*R1)/(R2+R1),0,-(R2*R1)/(R2+R1)];
    
    
    E=[ 1       -R1/(R1+R2)       0;...
        1       0                 0];
    
    F=[-R2/(R1+R2);...
        0];
    
    G=[-(R2*R1)/(R1+R2)        0            (R2*R1)/(R1+R2);...
        0                     0            0];
   
elseif system==2
    %obtenido a mano
    A = [ 0       0                                1/C1;...
        0       (1/(C2*R2))*((R1/(R2+R1))-1)     -1/C2;...
        -1/L     1/L                              0];
    
    B=[0;...
        1/(C2*(R2+R1));...
        0];
    
    Q=[-1/C1                 -1/C1     0;...
        R1/(C2*(R2+R1))       0        -R1/(C2*(R2+R1));...
        0                     0        0];
    
    C=[0,R1/(R2+R1),0];
    
    D=R2/(R2+R1);
    
    R=[(R2*R1)/(R2+R1),0,-(R2*R1)/(R2+R1)];
    
    E=[ 1       -R1/(R1+R2)       0;...
        1       0                 0];
    
    F=[-R2/(R1+R2);...
        0];
    
    G=[-(Rc+Rb+((R2*R1)/(R1+R2)))        -Rb            (Rc+((R2*R1)/(R1+R2)));...
        -Rb                              -(Rb+Re)        -Re];
elseif system==4
     %Obtenidas por simulink con capacitores
    load('niv3.mat');
    if frecuencia==261
        A = f_261.A;
        B = f_261.B;
        C = f_261.C;
        D = f_261.D;
    elseif frecuencia==440
        A = f_440.A;
        B = f_440.B;
        C = f_440.C;
        D = f_440.D;
    elseif frecuencia==1000
        A = f_1000.A;
        B = f_1000.B;
        C = f_1000.C;
        D = f_1000.D;
    elseif frecuencia==5000
        A = f_5000.A;
        B = f_5000.B;
        C = f_5000.C;
        D = f_5000.D;
    elseif frecuencia==10000
        A = f_10000.A;
        B = f_10000.B;
        C = f_10000.C;
        D = f_10000.D;
    elseif frecuencia==15000
        A = f_15000.A;
        B = f_15000.B;
        C = f_15000.C;
        D = f_15000.D;
    end
    
    E=C(2:3,:);
    F=D(2:3,4);
    G=D(2:3,1:3);
end

%% Funciones
fv = @(vd) [iec_diodo(vd(1))/Br;icc_diodo(vd(2))/Bf;ict_trans(iec_diodo(vd(1)),icc_diodo(vd(2)),vd(1))];

if size_system==3
    lim = @(z,miu,v) [A * [z(1);z(2);z(3)] + B * miu(1) + Q * [z(4);z(5);z(6)]; ...
        ((eye(size(G,2)) - (num_jacobian(fv,v,G)*G))^(-1) *...
        (num_jacobian(fv,v,G)) *...
        (E*A * [z(1);z(2);z(3)] + E*B * miu(1) + E*Q * [z(4);z(5);z(6)] + F * miu(2)))];
    
    x0 = [0,0,0]'; % condiciones iniciales
    
elseif size_system==5
    lim = @(z,miu,v) [A * [z(1);z(2);z(3);z(4);z(5)] + B * miu(1) + Q * [z(6);z(7);z(8)]; ...
        %((eye(size(G,2)) - (num_jacobian(fv,v,G)*G))^(-1) *...
        ((num_jacobian(fv,v,G)) *...
        (E*A *[z(1);z(2);z(3);z(4);z(5)]+ E*B * miu(1) + E*Q *  [z(6);z(7);z(8)]+ F * miu(2)))];
    
    x0 = [0,0,0,0,0]'; % condiciones iniciales
end

%% Inicialización y condiciones iniciales
I0 = [0,0,0]';
z0 = [x0;I0];
v0 = [0,0]';
u0 = 0;
ud = 0;
miu0 = [u0;ud];
y0 = 0;
v = v0;
x = x0;
I = I0;
z = z0; % vector de estado
miu = miu0; % vector de entradas
u = u0;
y = y0; % vector de salidas
% Arrays para almacenar las trayectorias de las variables de estado,
% entradas y salidas del sistema
FV = zeros(numel(fv([0;0])),K+1);
Z = zeros(numel(z),K+1);
V = zeros(numel(v),K+1);
U = zeros(numel(u),K+1);
Y = zeros(numel(y),K+1);
% Inicialización de arrays
Z(:,1) = z0;
V(:,1) = v0;
U(:,1) = u0;
Y(:,1) = y0;

%% Solución recursiva del sistema LTI
for k = 1:K
    
    %u=s(k);
    u=Vin;
    
    miu = [u; 0];
    v = E*x+F*u+G*I;%Calculo el valor actual de los voltajes.
    % Se actualiza el estado del sistema y su salida
    k1 = lim(z,miu,v);
    k2 = lim(z+(dt/2)*k1, miu, v);
    k3 = lim(z+(dt/2)*k2, miu, v);
    k4 = lim(z+dt*k3, miu, v);
    z = z + (dt/6)*(k1+2*k2+2*k3+k4);
   
    if size_system==3
        x = [z(1);z(2);z(3)];
        I = [z(4);z(5);z(6)];%Extraigo la aproximación de la corriente.
    elseif size_system==5
        x = [z(1);z(2);z(3);z(4);z(5)];
        I = [z(6);z(7);z(8)];%Extraigo la aproximación de la corriente.
    end
    
    v = E*x+F*u+G*I; %Calculo el voltaje de la siguiente iteración utilizando
    %la aproximación de la corriente de la siguiente
    %iteración.
    I = fv(v);       %Calculo la corriente real de la siguiente iteración
    
    y = C*x+D*u+R*I; %Extraigo la salida.
    z = [x;I];       %Actualizo la corriente
    
    % Se almacenan los valores actuales en los arrays
    FV(:,k+1) = fv(v);
    Z(:,k+1) = z;
    U(:,k+1) = u;
    V(:,k+1) = v;
    Y(:,k+1) = y;
end
%% grafica
t = t0:dt:tf1;  %definir tiempo
figure;
plot(t,Y,'-r',t,U,'-g');
legend('Vout','Vin')
ylabel('Voltaje (V)')
xlabel('Tiempo (s)')
title('Entra vs. salida')

figure;
plot(t,Z(4,:),'-r',t,Z(5,:),'-g',t,Z(6,:),'-b');
legend('Id1','Id2','Ict')
ylabel('Corriente (A)')
xlabel('Tiempo (s)')
title('Corriente de los diodos')

figure;
subplot(2,1,1);
plot(t,Z(1,:),'-r',t,Z(2,:),'-g');
legend('X1','X2')
ylabel('Voltaje (V)')
xlabel('Tiempo (s)')
title('Variables de estado')

subplot(2,1,2);
plot(t,Z(3,:));
legend('X3')
xlabel('tiempo(s)')
ylabel('corriente(A)')


%% Guardar Audio
Yout=Y;

switch(factor)
    case 2
        scale = 24;
    case 4
        scale = 12;
    case 8
        scale = 6;
    case 12
        scale = 4;
    case 16
        scale = 3;
end

% *********************************************************************
% A Partir de esta sección se deben de incluir los archivos de audio de
% LTspice ubicados en la carpeta ¨Archivos de audio de LTSPice para comparación¨
% si se desea obtener una comparación mediante el análisis
% de fourier o exportar los datos para realizar las comparaciones con el
% archivo llamado comp_niv
% *********************************************************************
%% Analisis de Fourier
if frecuencia==261
    audiowrite('mat_frec261.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec261.wav');
    [data2, fs2] = audioread('lt_frec261.wav');
    data2=downsample(data2,scale);
    
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==1
        save('data_261_niv1','data1')
    elseif system==2
        save('data_261_niv2','data1')
    end
elseif frecuencia==440
    audiowrite('mat_frec440.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec440.wav');
    [data2, fs2] = audioread('lt_frec440.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==1
        save('data_440_niv1','data1')
    elseif system==2
        save('data_440_niv2','data1')
    end
elseif frecuencia==1000
    audiowrite('mat_frec1000.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec1000.wav');
    [data2, fs2] = audioread('lt_frec1000.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==1
        save('data_1000_niv1','data1')
    elseif system==2
        save('data_1000_niv2','data1')
    end
elseif frecuencia==5000
    audiowrite('mat_frec5000.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec5000.wav');
    [data2, fs2] = audioread('lt_frec5000.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==1
        save('data_5000_niv1','data1')
    elseif system==2
        save('data_5000_niv2','data1')
    end
elseif frecuencia==10000
    audiowrite('mat_frec10000.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec10000.wav');
    [data2, fs2] = audioread('lt_frec10000.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==1
        save('data_10000_niv1','data1')
    elseif system==2
        save('data_10000_niv2','data1')
    end
elseif frecuencia==15000
    audiowrite('mat_frec15000.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec15000.wav');
    [data2, fs2] = audioread('lt_frec15000.wav');
    data2=downsample(data2,scale);
    
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==1
        save('data_15000_niv1','data1')
    elseif system==2
        save('data_15000_niv2','data1')
        
    end
end
%%
%MATLAB
m1 = length(data1);
n1 = pow2(nextpow2(m1));% number of samples
data_fft1 = fft(data1,n1);
f1 = (0:n1-1)*(fs1/n1); % frequency range
power1 = abs(data_fft1).^2/n1;

%LTSpice
m2 = length(data2);
n2 = pow2(nextpow2(m2));% number of samples
data_fft2 = fft(data2,n2);
f2 = (0:n2-1)*((fs2/scale)/n2);     % frequency range
power2 = abs(data_fft2).^2/n2;

figure;
plot(f1(5:floor(n1/2)),power1(5:floor(n1/2)))
xlabel('Frequency')
ylabel('Power')
hold on
plot(f2(5:floor(n2/2)),power2(5:floor(n2/2)))
legend('Matlab','LTSpice')
hold off

%A bajas frecuencias coincides las frecuencias fundamentales
%A alta frecuencia no coinciden pero no se presenta ruido en la simulacion
