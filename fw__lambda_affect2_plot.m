%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')

addpath('./synthetic_exp')
addpath('./utils')

addpath('./fw_core/')

output_dir = './results/lambda_affect/';
lines = {};
names = {};
counter = 1;

T = load(sprintf('%s/lambda_sensitivity_new_0.0.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\lambda=0$';
counter = counter+1;

T = load(sprintf('%s/lambda_sensitivity_new_0.01.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\lambda=0.01$';
counter = counter+1;

T = load(sprintf('%s/lambda_sensitivity_new_0.05.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\lambda=0.05$';
counter = counter+1;

T = load(sprintf('%s/lambda_sensitivity_new_0.1.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\lambda=0.1$';
counter = counter+1;

T = load(sprintf('%s/lambda_sensitivity_new_1.0.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\lambda=1.0$';
counter = counter+1;

% T = load(sprintf('%s/lambda_sensitivity_new_2.0.mat', output_dir ));
% lines{end+1} = mean(T.success_rate, 3);
% names{end+1} = '$\lambda=2$';
% counter = counter+1;

T = load(sprintf('%s/lambda_sensitivity_new_5.0.mat', output_dir ));
lines{end+1} = mean(T.success_rate, 3);
names{end+1} = '$\lambda=5$';
counter = counter+1;

% T = load(sprintf('%s/lambda_sensitivity_new_10.0.mat', output_dir ));
% lines{end+1} = mean(T.success_rate, 3);
% names{end+1} = '$\lambda=10$';
% counter = counter+1;


% T = load(sprintf('%s/lambda_sensitivity_SPA.mat', output_dir ));
% lines{end+1} = mean(T.success_rate, 3);
% names{end+1} = '$\texttt{SPA}$';
% counter = counter+1;

noise = T.noise;

markers = {'o',  '+',  '*', '.', 'x', '_', 's', 'd'};
figure('DefaultAxesFontSize', 18);
linewidth = 1.4;
% lines = mat2cell(success_rate, ones(1, size(success_rate, 1)),...
%     [size(success_rate, 2)]);

for i=1:numel(lines)
    plot(noise, lines{i}, sprintf('%s-', markers{i}), ...
        'DisplayName', names{i}, ...
        'LineWidth', linewidth,  ...
        'MarkerSize', 10 ...
        );
    hold on
end
xlabel('SNR');
ylabel('success rate');
legend('Location', 'southeast', 'Interpreter', 'latex');

% set(gca, 'FontSize', 15);

path_to_file = sprintf('%s/lambda_sensitivity.eps', output_dir);
saveas(gcf, path_to_file, 'epsc')
fprintf(sprintf('Exported to %s\n', path_to_file));
