[tiempo ,signal] = square_signal(500, 1.5, 0.02 ,1000);

plot(tiempo, signal)
title('Señal cuadrada, 1.5 [Vpp], 500 [Hz] ')
xlabel('Tiempo')
ylabel('Voltaje [V]')