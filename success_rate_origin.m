%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')
addpath('~/code/matlab/common/prob_tools')

addpath('baselines/FGNSR/matlab')
addpath('baselines/PrecondSPA/')
addpath('fw_core')
addpath('utils')

M=50; N=200; K=40;
hparams = [5:1:13];
num_trials = 20;
alpha = ones(K, 1);

success_count = [];

Tracking = struct;
Tracking.lambda = -3;
Tracking.M = M;
Tracking.N = N;
Tracking.K = K;
Tracking.num_trials = num_trials;
Tracking.hparams = hparams;

global_tic = tic;
output_dir = sprintf('./results/SNR/origin/');
if ~isdir(output_dir)
    mkdir(output_dir)
    fprintf('Created output directory\n');
end

fprintf('Running X2\n')
num_exps = numel(hparams);
for exp_ind=1:num_exps
    fprintf('Exp %d/%d -  Hparam: %f\n', exp_ind, num_exps, ...
        hparams(exp_ind))
    Tracking.fw(exp_ind).success_count = [];
    Tracking.fw(exp_ind).mse = [];
    Tracking.fw(exp_ind).mse = [];
    Tracking.fw(exp_ind).Werr = [];
    Tracking.fw(exp_ind).RAE = [];
    Tracking.fw(exp_ind).MRSA = [];

    Tracking.fg(exp_ind).success_count = [];
    Tracking.fg(exp_ind).mse = [];
    Tracking.fg(exp_ind).Werr = [];
    Tracking.fg(exp_ind).RAE = [];
    Tracking.fg(exp_ind).MRSA = [];

    Tracking.spa(exp_ind).success_count = [];
    Tracking.spa(exp_ind).mse = [];
    Tracking.spa(exp_ind).Werr = [];
    Tracking.spa(exp_ind).RAE = [];
    Tracking.spa(exp_ind).MRSA = [];

    Tracking.optimal(exp_ind).Werr = [];
    Tracking.optimal(exp_ind).mse = [];

    for trial_ind=1:num_trials
        if mod(trial_ind, 10) == 0
            fprintf('Trial %d\n', trial_ind);
            fprintf('Total elapsed time: %.2f s\n', duration);
        end
            
        param = hparams(exp_ind);
        [X, pure_pixel_set, W, H, V] = generate_data(M, N, K, 'SNR', param, 'alpha', alpha);

        W_hat = X(:, pure_pixel_set);
        Tracking.optimal(exp_ind).Werr = mse_measure(W_hat, W)*size(W, 2)/norm(W, 'fro')^2;
        H_hat = computeApproxError(X, pure_pixel_set);
        Tracking.optimal(exp_ind).mse(trial_ind) = norm(X- W_hat*H_hat, 'fro')/norm(X, 'fro');

        % ==============================
        % FRANK-WOLFE
        fw_tic = tic;
        options = struct;
        options.maxIters = 100;
        options.verbose = 0;
        options.lambda = Tracking.lambda;
        options.debug = 0;
        options.N = K;
        options.backend = 'mex';
        options.dualGapThreshold=0;
        options.epsilon = 1e-5;
        [C_hat, fw_tracking] = fw(X, options);
        [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);

        Tracking.fw(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.fw(exp_ind).duration(trial_ind) = toc(fw_tic);
        Tracking.fw(exp_ind).lambda(trial_ind) = fw_tracking.lambda;

        W_hat = X(:, Lambda_hat);
        H_hat = computeApproxError(X, Lambda_hat);
        Tracking.fw(exp_ind).Werr(trial_ind) = mse_measure(W_hat, W)*size(W, 2)/norm(W, 'fro')^2;
        Tracking.fw(exp_ind).mse(trial_ind) = norm(X- W_hat*H_hat, 'fro')/norm(X, 'fro');
        Tracking.fw(exp_ind).RAE(trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        Tracking.fw(exp_ind).MRSA(trial_ind) = compute_MRSA(W_hat, W);


        % ==============================
        % FAST GRADIENT
        fg_tic = tic;
        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        Tracking.fg(exp_ind).success_count(trial_ind) = ...
            all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        Tracking.fg(exp_ind).duration(trial_ind) = toc(fg_tic);

        Lambda_hat = K_fg;
        W_hat = X(:, Lambda_hat);
        H_hat = computeApproxError(X, Lambda_hat);
        Tracking.fg(exp_ind).Werr(trial_ind) = mse_measure(W_hat, W)*size(W, 2)/norm(W, 'fro')^2;
        Tracking.fg(exp_ind).mse(trial_ind) = norm(X- W_hat*H_hat, 'fro')/norm(X, 'fro');
        Tracking.fg(exp_ind).RAE(trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        Tracking.fg(exp_ind).MRSA(trial_ind) = compute_MRSA(W_hat, W);

        % ==============================
        % SPA
        i_tic = tic;
        Lambda_hat = SPAselect(X, K);
        Tracking.spa(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.spa(exp_ind).duration(trial_ind) = toc(i_tic);

        duration = toc(global_tic);
        Tracking.duration = duration;

        W_hat = X(:, Lambda_hat);
        H_hat = computeApproxError(X, Lambda_hat);
        Tracking.spa(exp_ind).Werr(trial_ind) = mse_measure(W_hat, W)*size(W, 2)/norm(W, 'fro')^2;
        Tracking.spa(exp_ind).mse(trial_ind) = norm(X- W_hat*H_hat, 'fro')/norm(X, 'fro');
        Tracking.spa(exp_ind).RAE(trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
        Tracking.spa(exp_ind).MRSA(trial_ind) = compute_MRSA(W_hat, W);
        % save('./results/Werr.mat')
    end
end

fw = arrayfun(@(x) mean(x.RAE), Tracking.fw);
fg = arrayfun(@(x) mean(x.RAE), Tracking.fg);
spa = arrayfun(@(x) mean(x.RAE), Tracking.spa);

figure();
linewidth = 1.4;
markersize = 9;
marker_pools = {'s', 'o', '^', 'x', 's', 'd', '^', 'v'};
linestyle_pools = {'-', '--', '-.'};

lines{1} = fw;
lines{2} = fg;
lines{3} = spa;
names = {'\texttt{MERIT}', '\texttt{FastGradient}', '\texttt{SPA}'};
for i=1:numel(lines)
    plot(hparams, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
        'DisplayName', names{i}, ...
        'LineWidth', linewidth,  ...
        'MarkerSize', markersize ...
        );
    hold on
end
xlabel('SNR');
ylabel('\texttt{Relative approximation error}', 'Interpreter',  'latex');
legend('Location', 'northwest', 'Interpreter', 'latex');
axis tight;
set(gca, 'FontSize', 15);
saveas(gcf, './results/old_setting/RAE.eps', 'epsc')


figure();
linewidth = 1.4;
markersize = 9;
marker_pools = {'s', 'o', '^', 'x', 's', 'd', '^', 'v'};
linestyle_pools = {'-', '--', '-.'};

fw = arrayfun(@(x) mean(x.MRSA), Tracking.fw);
fg = arrayfun(@(x) mean(x.MRSA), Tracking.fg);
spa = arrayfun(@(x) mean(x.MRSA), Tracking.spa);
lines{1} = fw;
lines{2} = fg;
lines{3} = spa;
names = {'\texttt{MERIT}', '\texttt{FastGradient}', '\texttt{SPA}'};
for i=1:numel(lines)
    plot(hparams, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
        'DisplayName', names{i}, ...
        'LineWidth', linewidth,  ...
        'MarkerSize', markersize ...
        );
    hold on
end
xlabel('SNR');
% yticks([0.1 0.30 1.0 5.0 10.0])
ylabel('MRSA');
legend('Location', 'northeast', 'Interpreter', 'latex');
axis tight;
set(gca, 'FontSize', 15);
saveas(gcf, './results/old_setting/MRSA.eps', 'epsc')




fprintf('Running X2\n')
fw = arrayfun(@(x) mean(x.success_count), Tracking.fw);
fg = arrayfun(@(x) mean(x.success_count), Tracking.fg);
spa = arrayfun(@(x) mean(x.success_count), Tracking.spa);

figure();
linewidth = 1.4;
markersize = 9;
marker_pools = {'s', 'o', '^', 'x', 's', 'd', '^', 'v'};
linestyle_pools = {'-', '--', '-.'};

lines{1} = fw;
lines{2} = fg;
lines{3} = spa;
names = {'\texttt{MERIT}', '\texttt{FastGradient}', '\texttt{SPA}'};
for i=1:numel(lines)
    plot(hparams, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
        'DisplayName', names{i}, ...
        'LineWidth', linewidth,  ...
        'MarkerSize', markersize ...
        );
    hold on
end
xlabel('SNR');
ylabel('success rate');
legend('Location', 'northwest', 'Interpreter', 'latex');
axis tight;
set(gca, 'FontSize', 15);
saveas(gcf, './results/old_setting/succ_rate.eps', 'epsc')

% % % ------------- MSE
% fw = arrayfun(@(x) mean(x.mse), Tracking.fw);
% fg = arrayfun(@(x) mean(x.mse), Tracking.fg);
% spa = arrayfun(@(x) mean(x.mse), Tracking.spa);
% figure();
% linewidth = 1.4;
% markersize = 9;
% marker_pools = {'s', 'o', '^', 'x', 's', 'd', '^', 'v'};
% linestyle_pools = {'-', '--', '-.'};
%
% lines{1} = fw;
% lines{2} = fg;
% lines{3} = spa;
% names = {'\texttt{MERIT}', '\texttt{FastGradient}', '\texttt{SPA}'};
% for i=1:numel(lines)
%     plot(hparams, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
%         'DisplayName', names{i}, ...
%         'LineWidth', linewidth,  ...
%         'MarkerSize', markersize ...
%         );
%     hold on
% end
% xlabel('SNR');
% % yticks([0.1 0.30 1.0 5.0 10.0])
% ylabel('MSE');
% legend('Location', 'northeast', 'Interpreter', 'latex');
% axis tight;
% set(gca, 'FontSize', 15);
% saveas(gcf, 'mse.eps', 'epsc')
%
%
%
% % W_Error
% fw = arrayfun(@(x) mean(x.Werr), Tracking.fw);
% fg = arrayfun(@(x) mean(x.Werr), Tracking.fg);
% spa = arrayfun(@(x) mean(x.Werr), Tracking.spa);
% optimal = arrayfun(@(x) mean(x.Werr), Tracking.optimal);
% figure();
% linewidth = 1.4;
% markersize = 9;
% marker_pools = {'s', 'o', '^', 'x', 's', 'd', '^', 'v'};
% linestyle_pools = {'-', '--', '-.', '-'};
%
% lines{1} = fw;
% lines{2} = fg;
% lines{3} = spa;
% % lines{4} = optimal;
% names = {'\texttt{MERIT}', '\texttt{FastGradient}', '\texttt{SPA}', 'optimal'};
% for i=1:numel(lines)
%     plot(hparams, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
%         'DisplayName', names{i}, ...
%         'LineWidth', linewidth,  ...
%         'MarkerSize', markersize ...
%         );
%     hold on
% end
% xlabel('SNR');
% % yticks([0.1 0.30 1.0 5.0 10.0])
% ylabel('MSE$_{W}$', 'interpreter', 'latex');
% legend('Location', 'northeast', 'Interpreter', 'latex');
% axis tight;
% set(gca, 'FontSize', 15);
% saveas(gcf, 'Werr.eps', 'epsc')
