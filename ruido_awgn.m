function  [signal_with_noise] = ruido_awgn(signal, A_signal, SNR)

N = length(signal);
signal_power_DB = 20*log10(A_signal);
noise_power_DB = signal_power_DB - SNR;

ruido = wgn(1,N, noise_power_DB, 1);
signal_with_noise = signal + ruido;

end