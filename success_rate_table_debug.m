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

% M=80; N=200; 
% list_K = [40, 50, 60, 70];

M=50; N=200; 
list_K = [10, 20, 30, 40];

list_SNR = [5];
num_trials = 1;

Tracking = struct;
Tracking.lambda = -3;
Tracking.M = M;
Tracking.N = N;
Tracking.num_trials = num_trials;

global_tic = tic;
% output_dir = './results/R1_table/';
% if ~isdir(output_dir)
%     mkdir(output_dir)
%     fprintf('Created output directory\n');
% end

fprintf('Running X2\n')
for i=1:numel(list_K)
    for j=1:numel(list_SNR)
        for trial_ind=1:num_trials
            if mod(trial_ind, 10) == 0
                fprintf('K=%d - SNR=%d - Trial %d\n', list_K(i), list_SNR(j), trial_ind);
                fprintf('Total elapsed time: %.2f s\n', duration);
            end
            K = list_K(i);
            SNR = list_SNR(j);
            alpha = 1*ones(K, 1);

            W = rand(M, K);
            H = zeros(K, N);
            H(:, 1:K) = eye(K);
            H(:, K+1:end) = dirichlet_rnd(alpha, N-K);
            Y = W*H;
            SNR = 10^(SNR/10);
            noise = randn(size(Y)); 
            sigma2 = sum(vecnorm(Y, 2, 1).^2) / M / N / SNR;
            noise = sqrt(sigma2)*noise;
            X = Y + noise;
            indices = randperm(N);
            X = X(:, indices);
            H = H(:, indices);
            r_pure_pixel_set = [];
            pure_pixel_set = 1:K;
            for ii=1:numel(pure_pixel_set)
                r_pure_pixel_set(end+1) = find(indices == pure_pixel_set(ii));
            end
            pure_pixel_set = r_pure_pixel_set;


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

            Tracking.fw(i, j, trial_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
            W_hat = X(:, Lambda_hat);
            H_hat = computeApproxError(X, Lambda_hat);
            Tracking.fw_W(i, j, trial_ind) = mse_measure(W_hat, W)*size(W, 2);
            Tracking.fw_RAE(i, j, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
            Tracking.fw_MRSA(i, j, trial_ind) = compute_MRSA(W_hat, W);


            options.lambda = 0;
            [C_hat, fw_tracking] = fw(X, options);
            [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);
            Tracking.fw0(i, j, trial_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
            W_hat = X(:, Lambda_hat);
            H_hat = computeApproxError(X, Lambda_hat);
            Tracking.fw0_W(i, j, trial_ind) = mse_measure(W_hat, W)*size(W, 2);
            Tracking.fw0_RAE(i, j, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
            Tracking.fw0_MRSA(i, j, trial_ind) = compute_MRSA(W_hat, W);

            % ==============================
            % FAST GRADIENT
            fg_tic = tic;
            [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
                'debug', 0, 'verbose', false);
            Tracking.fg(i, j, trial_ind) = all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
            W_hat = X(:, K_fg);
            Tracking.fg_W(i, j, trial_ind) = mse_measure(W_hat, W)*size(W, 2);
            H_hat = computeApproxError(X, K_fg);
            Tracking.fg_W(i, j, trial_ind) = mse_measure(W_hat, W)*size(W, 2);
            Tracking.fg_RAE(i, j, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
            Tracking.fg_MRSA(i, j, trial_ind) = compute_MRSA(W_hat, W);


            % ==============================
            % SPA
            i_tic = tic;
            Lambda_hat = SPAselect(X, K);
            Tracking.spa(i, j, trial_ind) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
            W_hat = X(:, Lambda_hat);
            Tracking.spa_W(i, j, trial_ind) = mse_measure(W_hat, W)*size(W, 2);
            H_hat = computeApproxError(X, Lambda_hat);
            Tracking.spa_W(i, j, trial_ind) = mse_measure(W_hat, W)*size(W, 2);
            Tracking.spa_RAE(i, j, trial_ind) = 1-norm(X - W_hat*H_hat, 'fro')/norm(X, 'fro');
            Tracking.spa_MRSA(i, j, trial_ind) = compute_MRSA(W_hat, W);

            duration = toc(global_tic);
            Tracking.duration = duration;

            % save(sprintf('%s/success_rate_%.2f.mat', output_dir, alpha(1)), 'Tracking')
        end
    end
end
% save(sprintf ('./results/R1_table/m4070.mat'))
% load(sprintf ('./results/R1_table/m4070.mat'))
% save(sprintf ('./results/R1_table/m1040.mat'))
load(sprintf ('./results/R1_table/m1040.mat'))

mean(Tracking.fw, 3)
mean(Tracking.fw0, 3)
mean(Tracking.fg, 3)
mean(Tracking.spa, 3)

mean(Tracking.fw_RAE, 3)
mean(Tracking.fw0_RAE, 3)
mean(Tracking.fg_RAE, 3)
mean(Tracking.spa_RAE, 3)

mean(Tracking.fw_MRSA, 3)
mean(Tracking.fw0_MRSA, 3)
mean(Tracking.fg_MRSA, 3)
mean(Tracking.spa_MRSA, 3)
% save('./results/R1_table/m405060.mat')
%
% mean(Tracking.fw, 3)'
% mean(Tracking.fg, 3)'
% mean(Tracking.spa, 3)'
% set(groot, 'defaultAxesTickLabelInterpreter','latex');
% set(groot, 'defaultLegendInterpreter','latex');
% set(groot, 'DefaultLineMarkerSize', 9)
% set(groot, 'DefaultLineLineWidth', 1.4)
% set(groot, 'DefaultAxesFontSize', 14);
%
% figure();
% plot(list_SNR, mean(Tracking.fw(1, :, :), 3), '-s', 'DisplayName',  '\texttt{MERIT}')
% hold on
% plot(list_SNR, mean(Tracking.fg(1, :, :), 3), '-o', 'DisplayName',  '\texttt{FastGradient}')
% plot(list_SNR, mean(Tracking.spa(1, :, :), 3), '-^', 'DisplayName',  '\texttt{SPA}')
% legend('Location',  'SouthEast')
% saveas(gcf, './results/R1_table/succ_K_40.eps', 'epsc')
%
% figure();
% plot(list_SNR, mean(Tracking.fw(2, :, :), 3), '-s', 'DisplayName',  '\texttt{MERIT}')
% hold on
% plot(list_SNR, mean(Tracking.fg(2, :, :), 3), '-o', 'DisplayName',  '\texttt{FastGradient}')
% plot(list_SNR, mean(Tracking.spa(2, :, :), 3), '-^', 'DisplayName',  '\texttt{SPA}')
% legend('Location',  'SouthEast')
% saveas(gcf, './results/R1_table/succ_K_50.eps', 'epsc')
%
% figure();
% plot(list_SNR, mean(Tracking.fw(3, :, :), 3), '-s', 'DisplayName',  '\texttt{MERIT}')
% hold on
% plot(list_SNR, mean(Tracking.fg(3, :, :), 3), '-o', 'DisplayName',  '\texttt{FastGradient}')
% plot(list_SNR, mean(Tracking.spa(3, :, :), 3), '-^', 'DisplayName',  '\texttt{SPA}')
% legend('Location',  'SouthEast')
% saveas(gcf, './results/R1_table/succ_K_60.eps', 'epsc')

