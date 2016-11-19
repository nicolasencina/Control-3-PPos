frecuencia_signal = 500;    % Frecuencia de 500 Hz
v_pp = 1.5;                 % Voltaje peak to peak
time = 0.3;                % Tiempo a graficar
n_puntos = 100000;            % Número de puntos

%sampling_period = time/(n_puntos - 1);
%sampling_freq = 1/sampling_period

[t ,signal, sampling_period] = square_signal(frecuencia_signal, v_pp, time, n_puntos);

% figure, plot(t, signal)

T = sampling_period;
Fs = 1/T;
L = n_puntos;

Y = fft(signal);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
% figure, plot(f,P1)
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f (Hz)')
% ylabel('|P1(f)|')
% xlim([0 3000]);

% integral=cumsum(signal)*T;
% figure
% subplot(1,2,1), plot(t,signal)
% subplot(1,2,2), plot(t,integral)

amplitud = 3/2;         %Vpp = 3 voltios
mod_frec = 10*1000;     %Frecuencia modulación 10 kHz
moduladora = signal;
kf = 628.3;
kf = 837.75;
modulated = F_modulator(amplitud, mod_frec, kf, moduladora, T, t);

figure, plot(t,modulated)
xlim([0 0.001])

SNR = 50;
%SNR = 30;
modulada_con_ruido = ruido_awgn(modulated, amplitud, SNR);

mod_con_ruido = awgn(modulated, SNR); 

figure, plot(t,mod_con_ruido)
xlim([0 0.001])

figure, plot(t,modulada_con_ruido)
xlim([0 0.001])

%Demodulación con función de Matlab
deltaf1=75;
deltaf2=100;
deltaf=deltaf1;        %Selección de desviación
z = fmdemod(modulated,mod_frec,Fs,deltaf)
figure
plot(t,z)




