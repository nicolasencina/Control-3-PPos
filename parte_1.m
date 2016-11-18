frecuencia_signal = 500;    % Frecuencia de 500 Hz
v_pp = 1.5;                 % Voltaje peak to peak
time = 0.005;                % Tiempo a graficar
n_puntos = 1000;            % Número de puntos

%sampling_period = time/(n_puntos - 1);
%sampling_freq = 1/sampling_period;

[t ,signal, sampling_period] = square_signal(frecuencia_signal, v_pp, time, n_puntos);
figure, plot(t, signal)

T = sampling_period;
Fs = 1/T
L = n_puntos;

Y = fft(signal);

P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;
figure, plot(f,P1)
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')

integral=cumsum(signal)*T;
figure
subplot(1,2,1), plot(t,signal)
subplot(1,2,2), plot(t,integral)





