function [modulated] = F_modulator(amplitud, mod_frec, kf, moduladora, T, time)

integral=cumsum(moduladora)*T;

frec_instantanea = 2*pi*mod_frec*time + kf*integral;
plot(frec_instantanea)
modulated = amplitud*cos(frec_instantanea);
end

