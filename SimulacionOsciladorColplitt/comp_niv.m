% Comparación entre niveles de complejidad
close all
frecuencia=261;
factor = 12;
fs=44100;

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

%%
if frecuencia==261
    load('data_261_niv1.mat');
    data_niv1=data1;
    load('data_261_niv2.mat');
    data_niv2=data1;
    load('data_261_niv3.mat');
    data_niv3=data1;
    [data2, fs2] = audioread('lt_frec261.wav');
    data2=downsample(data2,scale);
    data2=delayseq(data2,0,factor*fs);%261
    %data_niv3=delayseq(data_niv3,2.2e-3,factor*fs);%f2
    %data_niv3=delayseq(data_niv3,0.5e-3,factor*fs);%f4
    %data_niv3=delayseq(data_niv3,0.15e-3,factor*fs);%f8
    timeLimits =  [6.515200216e-1 6.833616739e-1];
elseif frecuencia==440
    load('data_440_niv1.mat');
    data_niv1=data1;
    load('data_440_niv2.mat');
    data_niv2=data1;
    load('data_440_niv3.mat');
    data_niv3=data1;
    [data2, fs2] = audioread('lt_frec440.wav');
    data2=downsample(data2,scale);
    data2=delayseq(data2,1.04e-3,factor*fs);%440
    timeLimits =  [6.529997542e-1 6.709209595e-1];%440
elseif frecuencia==1000
    load('data_1000_niv1.mat');
    data_niv1=data1;
    load('data_1000_niv2.mat');
    data_niv2=data1;
    load('data_1000_niv3.mat');
    data_niv3=data1;
    [data2, fs2] = audioread('lt_frec1000.wav');
    data2=downsample(data2,scale);
    data2=delayseq(data2,0.22e-3,factor*fs);%1000
    data_niv3=delayseq(data_niv3,-0.12e-3,factor*fs);
    timeLimits =  [8.665127134e-1 8.730365383e-1];%1000
elseif frecuencia== 5000
    load('data_5000_niv1.mat');
    data_niv1=data1;
    load('data_5000_niv2.mat');
    data_niv2=data1;
    load('data_5000_niv3.mat');
    data_niv3=data1;
    [data2, fs2] = audioread('lt_frec5000.wav');
    data2=downsample(data2,scale);
    data2=delayseq(data2,-0.1e-3,factor*fs);
    data_niv3=delayseq(data_niv3,0.009e-3,factor*fs);
    timeLimits =  [7.384226037e-1 7.399745317e-1];
elseif frecuencia==10000
    load('data_10000_niv1.mat');
    data_niv1=data1;
    load('data_10000_niv2.mat');
    data_niv2=data1;
    load('data_10000_niv3.mat');
    data_niv3=data1;
    [data2, fs2] = audioread('lt_frec10000.wav');
    data2=downsample(data2,scale);

    timeLimits =  [1.242803445 1.243597870];
    data2=delayseq(data2,-1.04e-05,factor*fs);%5000 y 10000
    data_niv3=delayseq(data_niv3,-2.005e-3,factor*fs);%10000
elseif frecuencia==15000
    load('data_15000_niv1.mat');
    data_niv1=data1;
    load('data_15000_niv2.mat');
    data_niv2=data1;
    load('data_15000_niv3.mat');
    data_niv3=data1;
    [data2, fs2] = audioread('lt_frec15000.wav');
    data2=downsample(data2,scale);

    timeLimits =  [1.732906175 1.733438813];
    data2=delayseq(data2,(-88235+2.15)/(12*fs),factor*fs);%5000 y 10000
    data_niv3=delayseq(data_niv3,0.007e-3,factor*fs);%10000
    %data_niv2=delayseq(data_niv2,0.94);
end

frequencyLimits = [0 22050]; % Hz
%%
data2_ROI = data2(:);
sampleRate = factor*fs; % Hz
startTime = 0; % seconds
minIdx = ceil(max((timeLimits(1)-startTime)*sampleRate,0))+1;
maxIdx = floor(min((timeLimits(2)-startTime)*sampleRate,length(data2_ROI)-1))+1;
data2_ROI = data2_ROI(minIdx:maxIdx);
[Pdata2_ROI, Fdata2_ROI] = pspectrum(data2,sampleRate, 'FrequencyLimits',frequencyLimits);

%%
data_niv1_ROI = data_niv1(:);
sampleRate = factor*fs; % Hz
startTime = 0; % seconds
minIdx = ceil(max((timeLimits(1)-startTime)*sampleRate,0))+1;
maxIdx = floor(min((timeLimits(2)-startTime)*sampleRate,length(data_niv1_ROI)-1))+1;
data_niv1_ROI = data_niv1_ROI(minIdx:maxIdx);
[Pdata_niv1_ROI, Fdata_niv1_ROI] = pspectrum(data_niv1,sampleRate,'FrequencyLimits',frequencyLimits);

%%
data_niv2_ROI = data_niv2(:);
sampleRate = factor*fs; % Hz
startTime = 0; % seconds
minIdx = ceil(max((timeLimits(1)-startTime)*sampleRate,0))+1;
maxIdx = floor(min((timeLimits(2)-startTime)*sampleRate,length(data_niv2_ROI)-1))+1;
data_niv2_ROI = data_niv2_ROI(minIdx:maxIdx);
[Pdata_niv2_ROI, Fdata_niv2_ROI] = pspectrum(data_niv2,sampleRate,'FrequencyLimits',frequencyLimits);
%%
data_niv3_ROI = data_niv3(:);
sampleRate = factor*fs; % Hz
startTime = 0; % seconds
minIdx = ceil(max((timeLimits(1)-startTime)*sampleRate,0))+1;
maxIdx = floor(min((timeLimits(2)-startTime)*sampleRate,length(data_niv3_ROI)-1))+1;
data_niv3_ROI = data_niv3_ROI(minIdx:maxIdx);
[Pdata_niv3_ROI, Fdata_niv3_ROI] = pspectrum(data_niv3,sampleRate,'FrequencyLimits',frequencyLimits);
%% Gráfica en amplitud
plot(data2_ROI,'--','LineWidth',1)
hold on
plot(data_niv1_ROI,'-o','MarkerIndices',1:505:length(data_niv1_ROI),'LineWidth',1)
plot(data_niv2_ROI,'-d','MarkerIndices',1:505:length(data_niv2_ROI),'LineWidth',1)
plot(data_niv3_ROI,'-s','MarkerIndices',1:505:length(data_niv3_ROI),'LineWidth',1)
hold off

xlabel('Tiempo(ms)')
ylabel('Voltaje (V)')
legend({'V_{out} LTspice','V_{out} Nivel 1','V_{out} Nivel 2','V_{out} Nivel 3'})
saveas(gcf,'niv_f12_261_amp.png')
%% Gráfica en frecuencia
figure;
pspectrum(data2_ROI,sampleRate, ...
    'FrequencyLimits',frequencyLimits);
hold on
pspectrum(data_niv1_ROI,sampleRate, ...
    'FrequencyLimits',frequencyLimits);
pspectrum(data_niv2_ROI,sampleRate, ...
    'FrequencyLimits',frequencyLimits);
pspectrum(data_niv3_ROI,sampleRate, ...
    'FrequencyLimits',frequencyLimits);
hold off
xlabel('Frecuencia(kHz)')
ylabel('PowerSpectrum (dB)') 
title('Análisis en frecuencia')

legend({'V_{out} LTspice','V_{out} Nivel 1','V_{out} Nivel 2','V_{out} Nivel 3'})
saveas(gcf,'niv_f12_261_frec.png')
%% RMSE
RMSE_amp1=sqrt(mean((data2_ROI-data_niv1_ROI).^2));
RMSE_frec1=sqrt(mean((Pdata2_ROI-Pdata_niv1_ROI).^2));


RMSE_amp2=sqrt(mean((data2_ROI-data_niv2_ROI).^2));
RMSE_frec2=sqrt(mean((Pdata2_ROI-Pdata_niv2_ROI).^2));

RMSE_amp3=sqrt(mean((data2_ROI-data_niv3_ROI).^2));
RMSE_frec3=sqrt(mean((Pdata2_ROI-Pdata_niv3_ROI).^2));

RMSE_aplitud=[RMSE_amp1,RMSE_amp2,RMSE_amp3]
RMSE_frecuencia=[RMSE_frec1,RMSE_frec2,RMSE_frec3]
%%
m1 = length(data2);
n1 = pow2(nextpow2(m1));% number of samples
data_fft1 = fft(data2,n1);
f1 = (0:n1-1)*((fs2/scale)/n1);     % frequency range
power1 = abs(data_fft1).^2/n1;

m2 = length(data_niv1);
n2 = pow2(nextpow2(m2));% number of samples
data_fft2 = fft(data_niv1,n2);
f2 = (0:n2-1)*(sampleRate/n2); % frequency range
power2 = abs(data_fft2).^2/n2;

m3 = length(data_niv2);
n3 = pow2(nextpow2(m3));% number of samples
data_fft3 = fft(data_niv2,n3);
f3 = (0:n3-1)*(sampleRate/n3); % frequency range
power3 = abs(data_fft3).^2/n3;

m4 = length(data_niv3);
n4 = pow2(nextpow2(m4));% number of samples
data_fft4 = fft(data_niv3,n4);
f4 = (0:n4-1)*(sampleRate/n4); % frequency range
power4 = abs(data_fft4).^2/n4;

figure;
plot(f1(5:floor(n1/32)),power1(5:floor(n1/32)))
xlabel('Frequency')
ylabel('Power')
hold on
plot(f2(5:floor(n2/32)),power2(5:floor(n2/32)))
plot(f3(5:floor(n3/32)),power3(5:floor(n3/32)))
plot(f4(5:floor(n4/32)),power4(5:floor(n4/32)))
legend({'V_{out} LTspice','V_{out} Nivel 1','V_{out} Nivel 2','V_{out} Nivel 3'})
hold off
