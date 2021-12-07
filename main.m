tau = 5;
sigma = 1;
T = 10
A = 1
A_list = [0.1, 0.2, 0.5, 1, 2, 5, 10]
% A_list = [0.1, 0.2, 0.5, 1]

num_trials = 200

mse_list = zeros(1, length(A_list));
crlb_list = zeros(1, length(A_list));
bias_list = zeros(1, length(A_list));
for ind=1:length(A_list)
    ind
    A = A_list(ind);
    error = zeros(1, num_trials);
    for i=1:num_trials
        [tau_hat, crlb] = estimate_tau(A, tau, sigma, T, 'signal_type', 'trapezoid');
        error(i) = tau_hat-tau;
    end
    bias_list(ind) = mean(error);
    mse_list(ind) = mean(error.^2);
    crlb_list(ind) = crlb;
end
snr = (A_list/sigma);


%%
figure(1)
t = tiledlayout(1, 2)
t.TileSpacing = 'compact';
t.Padding = 'compact';

snr = 10.*log10(snr);
nexttile
semilogy(snr, mse_list, linewidth=2);
hold on;
semilogy(snr, crlb_list, linewidth=2);
legend('MSE', 'CRLB')
% title("MSE vs. SNR")
xlabel('SNR', fontsize=16);
ylabel('MSE', fontsize=16);
ylim([min(mse_list)-1, max(mse_list)+1]);
ax=gca;
ax.FontSize = 16;

% nexttile
% semilogx(snr, crlb_list, linewidth=2)
% % title("MSE vs. SNR")
% xlabel('SNR', fontsize=16)
% ylabel('CRLB', fontsize=16)
% ylim([min(crlb_list)-1, max(crlb_list)+1])

% ax=gca;
% ax.FontSize = 16


nexttile
semilogx(snr, bias_list, linewidth=2);
% title("MSE vs. SNR");
xlabel('SNR', fontsize=16);
ylabel('Bias', fontsize=16);
ylim([min(bias_list)-1, max(bias_list)+1]);

ax=gca;
ax.FontSize = 16;


set(gcf, 'PaperPosition', [0 0 15 5]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [15 5]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'result/result_trapezoid', 'pdf') %Save figure


%%
% figure(2)
% t = tiledlayout(1, 2)
% t.TileSpacing = 'compact';
% t.Padding = 'compact';

% [~, ~, y, sig] = estimate_tau(5, tau, sigma, T, 'signal_type', 'sine');

% nexttile
% sig = cat(2,sig, zeros(1, length(y)-length(sig)));
% x = linspace(0, T, length(sig));
% plot(x, sig, linewidth=2);
% % title("Obser")
% xlabel('t', fontsize=16);
% ylabel('s(t)', fontsize=16);
% ylim([min(sig)-1, max(sig)+1]);
% ax=gca;
% ax.FontSize = 16;


% nexttile
% x = linspace(0, T, length(y));
% plot(x, y, linewidth=2);
% title("y(t) with A=5, tau=5 seconds");
% xlabel('t', fontsize=16);
% ylabel('observation', fontsize=16);
% ylim([min(y)-1, max(y)+1]);
% ax=gca;
% ax.FontSize = 16;

% set(gcf, 'PaperPosition', [0 0 10 6]); %Position plot at left hand corner with width 5 and height 5.
% set(gcf, 'PaperSize', [10 6]); %Set the paper to have width 5 and height 5.
% saveas(gcf, 'result/sinusoid', 'pdf'); %Save figure
