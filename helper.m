
classdef    helper 
methods     ( Static = true )
    function obj = objective(signal, observation)
        obj = sum((signal.*observation))^2;
    end
    
    function signal = sine(tau, delta, t, Ts, a, T)
        signal = cat(2, zeros(1, floor(tau/delta)), a*sin(2*pi*t/Ts));
    end
    
    function signal = trapezoid(tau, delta, t, Ts, a, T)
        signal = cat(2, 0:delta:a, a*ones(1,floor((Ts-2*a)/delta)), a:-delta:0);
        padding_start = zeros(1, floor(tau/delta));
        signal = cat(2, padding_start, signal);    
    end
    
    function crlb = CRLB_sine(amplitude, time_period, A, sigma)
        crlb = zeros(2);
        crlb(1,1) = (sigma^2*time_period)/(2*(A^2)*pi^2);
        crlb(2,2) = (2*sigma^2)/(time_period);
    end
    
    function crlb = CRLB_trapezoid(amplitude, time_period, A, sigma)
        crlb = zeros(2);
        crlb(1,1) = (sigma^2)/(2*(amplitude)*A^2);
        crlb(2,2) = (sigma^2)/(time_period-2*amplitude + 2*amplitude^3/3);
    end
    
end
end
