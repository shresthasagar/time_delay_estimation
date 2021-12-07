
clear;
clc;

tau = 5;
sigma = 1;
T = 10;
A = 5;
A_list = [ 0.5, 1, 2, 4, 6, 8, 10];
% A_list = [0.1, 0.2, 0.5, 1]

num_trials = 1000;
delta = 0.01;
Ts = 4;

t = 0:delta:Ts; % sampling time 
a = 1; % amplitude of the generated signal
SIGNAL_TYPE = 'trapezoid';

if isequal(SIGNAL_TYPE, 'trapezoid')
    s_t = helper.trapezoid(0,delta, t, Ts, a, T);
elseif isequal(SIGNAL_TYPE, 'sine')
    s_t = helper.sine(0,delta, t, Ts, a, T);
end

s_flip = flip(s_t);

for a_ind=1:length(A_list)
    A = A_list(a_ind)
    error = zeros(1,num_trials);
    for i=1:num_trials
        tau = 5 + randn(1)*0.01; 
        if isequal(SIGNAL_TYPE, 'trapezoid')
            s_delay = A*helper.trapezoid(tau, delta, t, Ts, a, T);
        elseif isequal(SIGNAL_TYPE, 'sine')
            s_delay = A*helper.sine(tau, delta, t, Ts, a, T); 
        end
        s_delay = s_delay + randn(size(s_delay))*(sigma/sqrt(delta));
        
        conv_out = conv(s_flip, s_delay);
        [m,b] = max(conv_out);
        tau_hat = b*delta - Ts - delta;
        error(i) = tau_hat - tau;
    end
    if isequal(SIGNAL_TYPE, 'trapezoid')
        crlb_mat = helper.CRLB_trapezoid(a, Ts, A, sigma);
        crlb(a_ind) = crlb_mat(1,1)
    elseif isequal(SIGNAL_TYPE, 'sine')
        crlb_mat = helper.CRLB_sine(a, Ts, A, sigma);
        crlb(a_ind) = crlb_mat(1,1);
    end
    bias_list(a_ind) = mean(error);
    mse(a_ind) = mean((error).^2);
end
mse 
crlb
%%
figure(1)
t = tiledlayout(1, 2)
t.TileSpacing = 'compact';
t.Padding = 'compact';

snr = (A_list.^2/sigma^2);
snr = 10.*log10(snr);
nexttile
semilogy(snr, mse, linewidth=2);
hold on;
semilogy(snr, crlb, linewidth=2);
legend('MSE', 'CRLB')
% title("MSE vs. SNR")
xlabel('SNR (dB)', fontsize=16);
ylabel('MSE/CRLB', fontsize=16);
ylim([0, max(mse)+1]);
ax=gca;
ax.FontSize = 16;


nexttile
semilogx(snr, bias_list, linewidth=2);
% title("MSE vs. SNR");
xlabel('SNR (dB)', fontsize=16);
ylabel('Bias', fontsize=16);
% ylim([min(bias_list)-1, max(bias_list)+1]);

ax=gca;
ax.FontSize = 16;


set(gcf, 'PaperPosition', [0 0 10 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [10 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'result/result_trapezoid', 'pdf') %Save figure
