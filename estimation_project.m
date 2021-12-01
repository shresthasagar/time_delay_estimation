tau = 5;
sigma = 1;
T = 10
A = 1
A_list = [0.1, 0.2, 0.5, 1, 2, 5, 10, 20, 50, 100]
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
t = tiledlayout(3, 1)
t.TileSpacing = 'compact';
t.Padding = 'compact';

nexttile
semilogx(snr, mse_list, linewidth=2)
% title("MSE vs. SNR")
xlabel('SNR', fontsize=16)
ylabel('MSE', fontsize=16)
ylim([min(mse_list)-1, max(mse_list)+1])
ax=gca;
ax.FontSize = 16

nexttile
semilogx(snr, crlb_list, linewidth=2)
% title("MSE vs. SNR")
xlabel('SNR', fontsize=16)
ylabel('CRLB', fontsize=16)
ylim([min(crlb_list)-1, max(crlb_list)+1])

ax=gca;
ax.FontSize = 16


nexttile
semilogx(snr, bias_list, linewidth=2)
% title("MSE vs. SNR")
xlabel('SNR', fontsize=16)
ylabel('Bias', fontsize=16)
ylim([min(bias_list)-1, max(bias_list)+1])

ax=gca;
ax.FontSize = 16


set(gcf, 'PaperPosition', [0 0 8 15]); %Position plot at left hand corner with width 5 and height 5.
set(gcf, 'PaperSize', [8 15]); %Set the paper to have width 5 and height 5.
saveas(gcf, 'project_result', 'pdf') %Save figure