function [time, signal, delta_T] = square_signal(frec, vpp, time, points)
    period = 1/frec;
    delta_T = time/(points-1);
    time = linspace(0, time, points);
    
    time_counter = 0;
    voltaje_actual = vpp/2;
    
    for i=1:points
        signal(i) = voltaje_actual;
        if time_counter > (period/2)
            time_counter = 0;
            if voltaje_actual == vpp/2
                voltaje_actual = -vpp/2;
            else
                voltaje_actual = vpp/2;
            end
        end
        time_counter = time_counter + delta_T;
    end
end 