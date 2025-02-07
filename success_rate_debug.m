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

M=100; N=400; K=70;
hparams = [5:1:10];
num_trials = 5;
% num_trials = 5;
alpha = 1*ones(K, 1);

success_count = [];

Tracking = struct;
Tracking.lambda = -3;
Tracking.M = M;
Tracking.N = N;
Tracking.K = K;
Tracking.num_trials = num_trials;
Tracking.hparams = hparams;

global_tic = tic;
output_dir = sprintf( './results/SNR/%d-%d/', M, K);
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
    Tracking.gp(exp_ind).success_count = [];
    Tracking.fg(exp_ind).success_count = [];
    Tracking.spa(exp_ind).success_count = [];

    for trial_ind=1:num_trials
        if mod(trial_ind, 10) == 0
            fprintf('Trial %d\n', trial_ind);
            fprintf('Total elapsed time: %.2f s\n', duration);
        end
            
        param = hparams(exp_ind);
        [X, pure_pixel_set, W, H, V] = generate_data(M, N, K, 'SNR', param, 'alpha', alpha);


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


        % ==============================
        % FAST GRADIENT
        fg_tic = tic;
        [C_hat_fg, K_fg, ~, ~, ~, ~, t] = fgnsr(X, K, 'maxiter', 150, ...
            'debug', 0, 'verbose', false);
        Tracking.fg(exp_ind).success_count(trial_ind) = ...
            all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        Tracking.fg(exp_ind).duration(trial_ind) = toc(fg_tic);


        % ==============================
        % SPA
        i_tic = tic;
        Lambda_hat = SPAselect(X, K);
        Tracking.spa(exp_ind).success_count(trial_ind) = ...
            all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
        Tracking.spa(exp_ind).duration(trial_ind) = toc(i_tic);

        duration = toc(global_tic);
        Tracking.duration = duration;


    end
end

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
% yticks([0.1 0.30 1.0 5.0 10.0])
ylabel('success rate');
legend('Location', 'northwest', 'Interpreter', 'latex');
axis tight;
set(gca, 'FontSize', 15);


