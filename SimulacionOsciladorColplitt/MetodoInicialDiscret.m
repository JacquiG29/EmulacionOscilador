clear
close all
%% Import Audio
factor = 12;
fs=44100;
%% Parametros simulación
t0 = 0; % tiempo inicial
tf1 = 1; % tiempo final
dt = 1/(factor*fs); % período de muestreo
K = (tf1-t0)/dt; % número de iteraciones
K = round(K);
t = t0:dt:tf1;

%% Parámetros del circuito
Vin=1;
system=3;%1 nivel 1, 2 nivel 2, 3 nivel 3 simulink
size_system=5; %3 para espacio de estados de 3x3 y 5 para espacios de estado de 5x5 
%Espacio de estados de 5x5 solo aplica para system=3
frecuencia=261;%Opciones de frecuencia: 261,440,1000,5000,10000,15000

R1 = 100;
R2 = 1*10^(3);
C1 = 10*10^(-6);
C2 = 1*10^(-6);
Ceq=(C1*C2)/(C1+C2);
L = inductancia(frecuencia,Ceq);

Ceq=(C1*C2)/(C1+C2)
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
    
    B=[Q,B];
    
    C=[0,R1/(R2+R1),0];
    
    D=R2/(R2+R1);
    
    R=[(R2*R1)/(R2+R1),0,-(R2*R1)/(R2+R1)];
    
    D=[R,D];
    
    E=[ 1       -R1/(R1+R2)       0;...
        1       0                 0];
    
    F=[-R2/(R1+R2);...
        0];
    
    G=[-(Rc+Rb+((R2*R1)/(R1+R2)))        -Rb            (Rc+((R2*R1)/(R1+R2)));...
        -Rb                              -(Rb+Re)        -Re];
elseif system==3
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
%Si cambio signos de Re logro obtener salida pero presenta ruido
%si dejo los signos diverge esto solo pasa para frecuencias bajas

%% Discretización del sistema
use_euler = false; % verdadero si se quiere discretizar con euler, falso si
% se quiere ZOH
if(use_euler)%Euler
    Ad = eye(length(A))+A*dt;
    Bd = B*dt;
    Cd = C;
    Dd = D;
else %ZOH
    Ad = expm(A*dt);
    fun=@(s)expm(s*A);%definir funcion a integrar
    Bd = (integral(fun,0,dt,'ArrayValued',true))*B;%opcion elegida porque se trabaja con matrices
    Cd = C;
    Dd = D;
end

%% Funciones
fv = @(vd) [iec_diodo(vd(1))/Br;icc_diodo(vd(2))/Bf;ict_trans(iec_diodo(vd(1)),icc_diodo(vd(2)),vd(1))];

if size_system==3
    x0 = [0,0,0]'; % condiciones iniciales
elseif size_system==5
    x0 = [0,0,0,0,0]'; % condiciones iniciales
end

%% Inicialización y condiciones iniciales
v0 = [0,0]';
u0 = [0;0;0;Vin];
I0 = [0,0,0]';
y0 =  C*x0;
v = v0;
x = x0;
I = I0;
u = u0;
y = y0; % vector de salidas
% Arrays para almacenar las trayectorias de las variables de estado,
% entradas y salidas del sistema
FV = zeros(numel(fv([0;0])),K+1);
X = zeros(numel(x),K+1);
U = zeros(numel(u),K+1);
Y = zeros(numel(y),K+1);
% Inicialización de arrays
X(:,1) = x0;
U(:,1) = u0;
Y(:,1) = y0;

%% Solución recursiva del sistema LTI
for k = 1:K
    v = E*x+F*u(4)+G*u(1:3);
    u(1:3)=fv(v);
    u(4)=Vin;
    
    % Se actualiza el estado del sistema y su salida
    x = Ad*X(:,k)+Bd*U(:,k);
    y = Cd*X(:,k)+Dd*U(:,k);
    
    % Se almacenan los valores actuales en los arrays
    FV(:,k+1) = fv(v);
    X(:,k+1) = x;
    U(:,k+1) = u;
    Y(:,k+1) = y;
end
%% grafica
t = t0:dt:tf1;  %definir tiempo
figure(1)

subplot(2,1,1);
plot(t,X(1,:),t,X(2,:));
title('Variables de estado')
legend('V_c1','V_c2')
xlabel('tiempo(s)')
ylabel('voltaje(V)')

subplot(2,1,2);
plot(t,X(3,:));
legend('i_L')
xlabel('tiempo(s)')
ylabel('corriente(mA)')

figure(3)

subplot(3,1,1);
plot(t,Y(1,:));
title('Voltaje de salida')

subplot(3,1,2);
plot(t,U(1,:),t,U(2,:),t,U(3,:));
legend('Id1','Id2','Ict')
title('Entradas de corriente')

subplot(3,1,3);
plot(t,U(4,:))
title('Vin')

%% Guardar Audio
Yout=Y(1,:);

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
    audiowrite('mat_frec261_zoh.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec261_zoh.wav');
    [data2, fs2] = audioread('lt_frec261.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==3
        save('data_261_niv3','data1')
    end
elseif frecuencia==440
    audiowrite('mat_frec440_zoh.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec440_zoh.wav');
    [data2, fs2] = audioread('lt_frec440.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==3
        save('data_440_niv3','data1')
    end
elseif frecuencia==1000
    audiowrite('mat_frec1000_zoh.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec1000_zoh.wav');
    [data2, fs2] = audioread('lt_frec1000.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==3
        save('data_1000_niv3','data1')
    end
elseif frecuencia==5000
    audiowrite('mat_frec5000_zoh.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec5000_zoh.wav');
    [data2, fs2] = audioread('lt_frec5000.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==3
        save('data_5000_niv3','data1')
    end
elseif frecuencia==10000
    audiowrite('mat_frec10000_zoh.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec10000_zoh.wav');
    [data2, fs2] = audioread('lt_frec10000.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,D);
    if system==3
        save('data_10000_niv3','data1')
    end
elseif frecuencia==15000
    audiowrite('mat_frec15000_zoh.wav',Yout,factor*fs)
    [data1, fs1] = audioread('mat_frec15000_zoh.wav');
    [data2, fs2] = audioread('lt_frec15000.wav');
    data2=downsample(data2,scale);
    D = finddelay(data1,data2)
    data1=delayseq(data1,5);
    if system==3
        save('data_15000_niv3','data1')
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

figure(4)
plot(f1(5:floor(n1/64)),power1(5:floor(n1/64)))
xlabel('Frequency')
ylabel('Power')
hold on
plot(f2(5:floor(n2/64)),power2(5:floor(n2/64)))
legend('Matlab','Original')
hold off

%A bajas frecuencias coincides las frecuencias fundamentales
%A alta frecuencia no coinciden pero no se presenta ruido en la simulacion