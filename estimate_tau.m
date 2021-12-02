function [tau, crlb, y, signal_raw] = estimate_tau(A, tau, sigma, T, varargin)
    Ts = 5;
    T_total = T;
    delta = 0.1; % sampling resolution
    t = 0:delta:Ts; % sampling time 
    a = 1; % amplitude of the generated signal
    SIGNAL_TYPE = 'trapezoid';

    % Read the optional parameters
    if (rem(length(varargin),2)==1)
        error('Optional parameters should always go by pairs');
    else
        for i=1:2:(length(varargin)-1)
            switch upper(varargin{i})
                case 'DELTA'
                    delta = varargin{i+1};
                case 'SIGNAL_TYPE'
                    varargin{i+1};
                    SIGNAL_TYPE = varargin{i+1};
                otherwise
                    % Hmmm, something wrong with the parameter string
                    error(['Unrecognized option: ''' varargin{i} '''']);
            end;
        end;
    end;


    if isequal(SIGNAL_TYPE, 'trapezoid')
        signal = trapezoid(tau, delta, t, Ts, a, T);
        signal_raw = trapezoid(0, delta, t, Ts, a, T);
        crlb_mat = CRLB_trapezoid(a, Ts, A, sigma);
        crlb = crlb_mat(1,1);
    elseif isequal(SIGNAL_TYPE, 'sine')
        signal = sine(tau, delta, t, Ts, a, T);
        signal_raw = sine(0, delta, t, Ts, a, T);
        crlb_mat = CRLB_sine(a, Ts, A, sigma);
        crlb = crlb_mat(1,1);    
    end

    % observation
    y = A*signal + randn(1, length(signal))*sigma;
    
    gridsize = 0.1*sqrt(crlb);
    tau_grid = 0:gridsize:10;
    tau_grid(2:end) = tau_grid(2:end) + (rand(1,length(tau_grid)-1) - 0.5)*gridsize;
    obj_list = zeros(1, length(tau_grid));

    for i=1:length(tau_grid)
        tau_candidate = tau_grid(i);
        signal_candidate = trapezoid(tau_candidate, delta, t, Ts, a, T);

        % match the size of the singal_candidate with the observation by padding it with zeros
        max_length = max(length(y), length(signal_candidate));
        signal_candidate = cat(2, signal_candidate, zeros(1, max_length-length(signal_candidate)));
        observation = cat(2, y, zeros(1, max_length-length(y)));

        obj = objective(signal_candidate, observation);
        obj_list(i) = obj;
    end
    [val, max_ind] = max(obj_list);
    tau = tau_grid(max_ind);
end

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
    crlb(1,1) = (2*sigma^2)/((amplitude^2)*(A^2)*time_period);
    crlb(2,2) = (2*sigma^2)/((amplitude^2)*time_period);
end

function crlb = CRLB_trapezoid(amplitude, time_period, A, sigma)
    crlb = zeros(2);
    crlb(1,1) = (sigma^2)/(2*(amplitude)*A^2);
    crlb(2,2) = (sigma^2)/(amplitude*(time_period-amplitude));
end
