%% Clear all things
clc; clear; close all; path(pathdef);
addpath('~/code/matlab/common')
addpath('~/code/matlab/common/visualization/lines')
addpath('~/code/matlab/common/PGD')
addpath('~/code/matlab/common/prox_ops')
addpath('~/code/matlab/common/prob_tools')

addpath('./synthetic_exp')
addpath('./utils')

addpath('./fw_core/')

output_dir = './results/lambda_affect/';
codeName = 'new_all';
if true
    fprintf(sprintf('%s/lambda_sensitivity_%s.mat \n', output_dir, codeName));

    M = 50; N = 200; K = 40;
    noise = [2:0.5:25];
    
    numTrial = 50;
    success_rate = [];
    hparams = [0, 0.01, 0.05, 0.1, 1.0 , 5];

    totalRuns = numel(noise)*numel(hparams)*numTrial;
    counter = 1;

    t0 = tic;
    for i=1:numel(hparams)
        for noiseIndex=1:numel(noise)
            for trialIndex=1:numTrial
                % [X, pure_pixel_set, S] = generate_data(1, 3, M, N, L, 'SNR', noise(noiseIndex));
                W = rand(M, K);
                H = zeros(K, N);
                H(:, 1:K) = eye(K);
                H(:, K+1:end) = dirichlet_rnd(ones(K, 1), N-K);
                Y = W*H;
                SNR = 10^(noise(noiseIndex)/10);
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


                if mod(counter, 100) == 0
                    fprintf('Current run %d/%d  - Time: %f \n', counter, totalRuns, toc(t0));
                    
                end
                counter = counter + 1;

                options = struct;
                options.maxIters = 100;
                options.N = K;
                options.checkpointIterval = 1;
                options.firstIterMethod = 'NONE';
                options.verbose = 0;
                options.innerVerboseInterval = 6000;
                options.numThreads = 4;
                options.epsilon = 1e-5;

                options.backend = 'mex';

                options.lambda = hparams(i);
                [C_hat, Tracking] = fw(X, options);
                [~, Lambda_hat] = maxk(vecnorm(C_hat, Inf, 2), K, 1);
                success_rate(i, noiseIndex, trialIndex) = all(sort(Lambda_hat(:)) == sort(pure_pixel_set(:)));
            end
            save(sprintf('%s/lambda_sensitivity_%s.mat', output_dir, codeName));
        end
    end

end
