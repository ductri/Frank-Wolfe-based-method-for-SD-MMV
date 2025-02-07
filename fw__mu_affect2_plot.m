%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')

addpath('./synthetic_exp')
addpath('./topic_modeling/baselines/PrecondSPA')
addpath('./utils')

addpath('./fw_core/')

output_dir = './results/mu_affect/';
lines = {};
names = {};
counter = 1;

T = load(sprintf('%s/mu_sensitivity_new_1e-1.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\mu=1\mathrm{e}{-1}$';
% T = load(sprintf('%s/mu_sensitivity_pre_1e-1.mat', output_dir ));
% lines{counter} = [mean(T.success_rate, 3) lines{counter}];
% T = load(sprintf('%s/mu_sensitivity_post_1e-1.mat', output_dir ));
% lines{counter} = [ lines{counter} mean(T.success_rate, 3)];
counter = counter+1;


T = load(sprintf('%s/mu_sensitivity_new_1e-3.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\mu=1\mathrm{e}{-3}$';
% T = load(sprintf('%s/mu_sensitivity_pre_1e-3.mat', output_dir ));
% lines{counter} = [mean(T.success_rate, 3) lines{counter}];
% T = load(sprintf('%s/mu_sensitivity_post_1e-3.mat', output_dir ));
% lines{counter} = [ lines{counter} mean(T.success_rate, 3)];
counter = counter+1;

T = load(sprintf('%s/mu_sensitivity_new_1e-5.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\mu=1\mathrm{e}{-5}$';


% T = load(sprintf('results/lambda_affect/lambda_sensitivity_SPA.mat', output_dir ));
% lines{end+1} = mean(T.success_rate, 3);
% names{end+1} = '$\texttt{SPA}$';
%
noise = T.noise;
% T = load(sprintf('%s/mu_sensitivity_pre_1e-5.mat', output_dir ));
% lines{counter} = [mean(T.success_rate, 3) lines{counter}];
% noise = [T.noise noise];
% T = load(sprintf('%s/mu_sensitivity_post_1e-5.mat', output_dir ));
% lines{counter} = [ lines{counter} mean(T.success_rate, 3)];
% noise = [noise T.noise];

% T = load(sprintf('%s/mu_sensitivity_1e-7.mat', output_dir ));
% lines{end+1} = mean(T.success_rate, 3);
% names{end+1} = '$\mu=1e-7$';



figure('DefaultAxesFontSize', 18);
linewidth = 1.4;
% lines = mat2cell(success_rate, ones(1, size(success_rate, 1)),...
%     [size(success_rate, 2)]);

markers = {'o',  '+',  '*', '.', 'x', '_', 's', 'd'};
for i=1:numel(lines)
    plot(noise, lines{i}, sprintf('%s-', markers{i}), ...
        'DisplayName', names{i}, ...
        'LineWidth', linewidth,  ...
        'MarkerSize', 7 ...
        );
    hold on
end
xlabel('SNR');
ylabel('success rate');
legend('Location', 'southeast', 'Interpreter', 'latex');

% set(gca, 'FontSize', 15);

path_to_file = sprintf('%s/mu_sensitivity.eps', output_dir);
saveas(gcf, path_to_file, 'epsc')
fprintf(sprintf('Exported to %s\n', path_to_file));
