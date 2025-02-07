%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD/')
addpath('~/code/matlab/common/prox_ops')

addpath('baselines/FGNSR/matlab')
addpath('baselines/PrecondSPA/')
addpath('fw_core')
addpath('utils')

%% UNIFORM A AND DIRICHLET S AND GAUSSIAN NOISE
% M = 40; N = 300; K = 20;
M=50; N=200; K=40;
hparams = [0.1:0.5:5];
num_trials = 50;

success_count = [];

Tracking = struct;
Tracking.M = M;
Tracking.N = N;
Tracking.K = K;
Tracking.num_trials = num_trials;
Tracking.hparams = hparams;

global_tic = tic;
output_dir = './';
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
        end
            
        param = hparams(exp_ind);
        [X, pure_pixel_set, W, H, V] = generate_data(M, N, K, 'SNR', 8);


        % ==============================
        % FRANK-WOLFE
        fw_tic = tic;
        options = struct;
        options.maxIters = 100;
        options.verbose = 0;
        options.lambda = hparams(exp_ind);
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
            'debug', 0, 'verbose', false, 'mu', hparams (exp_ind));
        Tracking.fg(exp_ind).success_count(trial_ind) = ...
            all(sort(K_fg(:)) == sort(pure_pixel_set(:)));
        Tracking.fg(exp_ind).duration(trial_ind) = toc(fg_tic);

        save(sprintf('%s/lambda.mat', output_dir), 'Tracking')
    end
    toc(global_tic)
end

load(sprintf('%s/lambda.mat', output_dir))
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
names = {'\texttt{MERIT}', '\texttt{FastGradient}'};
for i=1:numel(lines)
    plot(hparams, lines{i}, sprintf('%s%s', marker_pools{i}, linestyle_pools{i}), ...
        'DisplayName', names{i}, ...
        'LineWidth', linewidth,  ...
        'MarkerSize', markersize ...
        );
    hold on
end
xlabel('$\lambda$ in \texttt{MERIT} or $\mu$ in \texttt{FastGradient}', 'Interpreter', 'latex');
% yticks([0.1 0.30 1.0 5.0 10.0])
ylabel('success rate');
legend('Location', 'northwest', 'Interpreter', 'latex');
axis tight;
set(gca, 'FontSize', 15);
saveas(gcf, 'lambda.eps', 'epsc')
exportgraphics(gcf, 'lambda.png', 'resolution', 300);


