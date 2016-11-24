function y = fhs(signal,Fs,alfa)
    L=length(signal);
    Y = fft(signal);

    P2 = abs(Y/L);
    P1 = P2(1:L/2+1);
    P1(2:end-1) = 2*P1(2:end-1);
    
    f = Fs*(0:(L/2))/L;
    figure, plot(f,P1)
    title('Single-Sided Amplitude Spectrum of X(t),2')
    xlabel('f (Hz)')
    ylabel('|P1(f)|')
    flim=floor(Fs/2);
    lim=flim*alfa;
    xlim([0 lim]);
end