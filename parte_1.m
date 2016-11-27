%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 1)SEÑAL CUADRADA
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Generación de señal cuadrada

frecuencia_signal = 500;    % Frecuencia de 500 Hz
v_pp = 1.5;                 % Voltaje peak to peak
time = 0.3;                % Tiempo a graficar
n_puntos = 100000;            % Número de puntos

[t ,signal, sampling_period] = square_signal(frecuencia_signal, v_pp, time, n_puntos);

figure, plot(t, signal)
xlim([0 0.02])

T = sampling_period;
Fs = 1/T;
L = n_puntos;
%--------------------------------------------------------------------------
%Display fft por código(mitad del espectro positiva , sin réplica)

Y = fft(signal);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure, plot(f,P1)
title('Espectro de señal moduladora')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 3000]);
%--------------------------------------------------------------------------
%Display fft por función (mitad del espectro positiva , sin réplica)

fhs(signal,Fs,0.02)
title(('Espectro de señal moduladora'))
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 2)MODULACIÓN DE LA SEÑAL
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Generación de señal modulada

amplitud = 3/2;         %Vpp = 3 voltios
mod_frec = 10*1000;     %Frecuencia modulación 10 kHz
moduladora = signal;

deltaf1=3000;
deltaf2=100;
deltaf=deltaf1;
kf = k_f (deltaf, v_pp/2);

modulated = F_modulator(amplitud, mod_frec, kf, moduladora, T, t);
%--------------------------------------------------------------------------
%Display de señal modulada

%En el tiempo
figure, plot(t,modulated), title('Señal modulada')
xlim([0 0.01])
ylim([-2 2])

%En frecuencia
fhs(modulated,Fs,0.3)
title('Espectro de la señal modulada')


%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 3)ADICIÓN DE AWGN
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Generación de señal ruidosa

%SNR = 15;        %Niveles de ruido
SNR = 30;
modulada_con_ruido = ruido_awgn(modulated, amplitud, SNR); %Fn propia

mod_con_ruido = awgn(modulated, SNR);                      %Fn predefinida
%--------------------------------------------------------------------------
%Display de señal ruidosa (matlab y propia)

%En el tiempo
figure, plot(t,mod_con_ruido), title('Señal modulada con ruido gaussiano')
xlim([0 0.02])

figure, plot(t,modulada_con_ruido),  title('Señal modulada con ruido gaussiano')
xlim([0 0.02])

%En frecuencia
fhs(mod_con_ruido,Fs,0.3)
title('Espectro de la señal modulada con ruido')

fhs(modulada_con_ruido,Fs,0.3)
title('Espectro de la señal modulada con ruido')
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
% 4)DEMODULACIÓN DE LA SEÑAL
%--------------------------------------------------------------------------
%--------------------------------------------------------------------------
%Filtrado pasa bajos 

input=modulada_con_ruido;
d = fdesign.lowpass('N,Fc',100,20000,Fs); %Parámetros: orden de filtro y frecuencia de corte
Hd = design(d);
output = filter(Hd,input);
fvtool(Hd)                                         %Respuesta del filtro
plot(psd(spectrum.periodogram,output,'Fs',Fs))    %Espectro post-filtro
figure,plot(t,input)
xlim([0 0.02])
figure,plot(t,output)
xlim([0 0.02])
%-------------------------------------------------------------------------
%Filtrado pasa banda

d2 = fdesign.bandpass('N,F3dB1,F3dB2',2,7000,14000,Fs); %Parámetros: orden de filtro y frecuencias de corte
Hd2 = design(d2);
output2 = filter(Hd2,input);
fvtool(Hd2)                                         %Respuesta del filtro
plot(psd(spectrum.periodogram,output2,'Fs',Fs))    %Espectro post-filtro
figure,plot(t,input)
xlim([0 0.02])
figure,plot(t,output2)
xlim([0 0.02])
%-------------------------------------------------------------------------
%Segundo intento: filtros youtube

%generate the noisy signal which will be filtered /señal a filtrar
x=input;
% This script is available at https://dadorran.wordpress.com search for
% filtering matlab demo
% plot(x)
% title('Noisy signal')
% xlabel('Samples');
% ylabel('Amplitude')
 
%plot magnitude spectrum of the signal
X_mags = abs(fft(x));
X_mags = X_mags/(1.1*max(X_mags));
figure
plot(X_mags)
xlabel('DFT Bins')
ylabel('Magnitude')
 
% plot first half of DFT (normalised frequency)
num_bins = length(X_mags);
plot([0:1/(num_bins/2 -1):1], X_mags(1:num_bins/2))
xlabel('Normalised frequency (\pi rads/sample)')
ylabel('Magnitude')
 
% %design a second order filter using a butterworth design technique 
% [b a] = butter(2, 0.3, 'low')
%  
% %plot the frequency response (normalised frequency)
% H = freqz(b,a, floor(num_bins/2));
% hold on
% plot([0:1/(num_bins/2 -1):1], abs(H),'r');
%  
% %filter the signal using the b and a coefficients obtained from
% %the butter filter design function
% x_filtered = filter(b,a,x);
%  
% %plot the filtered signal
% figure(2)
% plot(x_filtered,'r')
% title('Filtered Signal - Using Second Order Butterworth')
% xlabel('Samples');
% ylabel('Amplitude')
%  
% %Redesign the filter using a higher order filter
% [b2 a2] = butter(20, 0.3, 'low')
%  
% %Plot the magnitude spectrum and compare with lower order filter
% H2 = freqz(b2,a2, floor(num_bins/2));
% figure(10)
% hold on
% plot([0:1/(num_bins/2 -1):1], abs(H2),'g');
%  
% %filter the noisy signal and plot the result
% x_filtered2 = filter(b2,a2,x);
% figure
% plot(x_filtered2,'g')
% title('Filtered Signal - Using 20th Order Butterworth')
% xlabel('Samples');
% ylabel('Amplitude')
 
%Use a band reject filter in place of a low pass filter
[b_stop a_stop] = butter(1, [0.05 0.07], 'bandpass'); %Parámetros: Orden y frecs de corte
 
%plot the magnitude spectrum
H_stopband = freqz(b_stop,a_stop, floor(num_bins/2));
hold on
plot([0:1/(num_bins/2 -1):1], abs(H_stopband),'c');
 
%plot filtered signal
x_filtered_stop = filter(b_stop,a_stop,x);
figure;
plot(x_filtered_stop,'c')
title('Filtered Signal - Using Bandpass')
xlabel('Samples');
ylabel('Amplitude')
 
% %Use matlabs built-in buttord function to get the optimum order to meet a specification
% [N Wn] = buttord(0.1, 0.5, 5, 40)
%  
% %use the N and Wn values obtained above to design the filter in the usual way
% [b3 a3] = butter(N, Wn, 'low');
%  
% %plot the magnitude spectrum
% H3 = freqz(b3,a3, floor(num_bins/2));
% figure(10);
% hold on
% plot([0:1/(num_bins/2 -1):1], abs(H2),'k');
% figure(10)
%  
% %filter the signal and plot the ouput of the filter
% x_filtered3 = filter(b3,a3,x);
% figure(5);
% plot(x_filtered3,'k')
% title(['Filtered Signal - Using ' num2str(N) ' th Order Butterworth'])
% xlabel('Samples');
% ylabel('Amplitude')
% % script available from https://dadorran.wordpress.com
% % comparison with other filter design techniques (chebyshev and elliptical)
% [b_butter a_butter] = butter(4, 0.2, 'low');
% H_butter = freqz(b_butter, a_butter);
%  
% [b_cheby a_cheby] = cheby1(4, 0.5, 0.2, 'low');
% H_cheby = freqz(b_cheby, a_cheby);
%  
% [b_ellip a_ellip] = ellip(4, 0.5, 40, 0.2, 'low');
% H_ellip = freqz(b_ellip, a_ellip);
%  
% %plot each filter to compare 
% figure(11)
% norm_freq_axis = [0:1/(512 -1):1];
% plot(norm_freq_axis, abs(H_butter))
% hold on
% plot(norm_freq_axis, abs(H_cheby),'r')
% plot(norm_freq_axis, abs(H_ellip),'g')
% legend('Butterworth', 'Chebyshev', 'Elliptical')
% xlabel('Normalised Frequency');
% ylabel('Magnitude')
%  
% %plot in dB for verification that spec is met
% figure(12);
% plot(norm_freq_axis, 20*log10(abs(H_butter)))
% hold on
% plot(norm_freq_axis, 20*log10(abs(H_cheby)),'r')
% plot(norm_freq_axis, 20*log10(abs(H_ellip)),'g')
% legend('Butterworth', 'Chebyshev', 'Elliptical')
% xlabel('Normalised Frequency ');
% ylabel('Magnitude (dB)')

%--------------------------------------------------------------------------
%Demodulación con función de Matlab (descomentar la que se quiere filtrar)

% demod_in=modulada_con_ruido;
demod_in=output;
% demod_in=output2;
% demod_in=x_filtered_stop


z = fmdemod(demod_in,mod_frec,Fs,deltaf);
figure
plot(t,z), title('Señal demodulada mediante fmdemod')
xlim([0 0.02])
%--------------------------------------------------------------------------
%Demodulación por bloques

%Derivada de la señal
deri=diff(demod_in);
t2=t(1:length(t)-1);
figure
plot(t2,deri)
title('Derivada de la señal modulada')
xlabel('Tiempo')
ylabel('Amplitud')
xlim([0 0.02])

%plot magnitude spectrum of the signal
x2=abs(deri);
X2_mags = abs(fft(x2));
X2_mags = X2_mags/(1.1*max(X2_mags));
figure
plot(X2_mags)
xlabel('DFT Bins')
ylabel('Magnitude')
 
%plot first half of DFT (normalised frequency)
num_bins2 = length(X2_mags);
plot([0:1/(num_bins2/2 -1):1], X_mags(1:num_bins2/2))
xlabel('Normalised frequency (\pi rads/sample)')
ylabel('Magnitude')

%Redesign the filter using a higher order filter
[b2 a2] = butter(1, 0.01, 'low') %Orden y frec de corte
 
%Plot the magnitude spectrum and compare with lower order filter
H2 = freqz(b2,a2, floor(num_bins2/2));
hold on
plot([0:1/(num_bins2/2 -1):1], abs(H2),'g');
 
%filter the noisy signal and plot the result
x_filtered2 = filter(b2,a2,x2);
figure
plot(x_filtered2,'g')
title('Filtered Signal - Using 20th Order Butterworth')
xlabel('Samples');
ylabel('Amplitude')

d3 = fdesign.lowpass('N,Fc',200,500,Fs); % Orden y frec de corte
designmethods(d3)
Hd3 = design(d3);
output3 = filter(Hd3,x2);
fvtool(Hd3)                                         %Respuesta del filtro
plot(psd(spectrum.periodogram,output3,'Fs',Fs))    %Espectro post-filtro
figure,plot(t2,x2)
xlim([0 0.02])
figure,plot(t2,output3)
xlim([0 0.02])






 





