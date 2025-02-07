function [X, pure_pixel_set, W, H, noise] = generate_data(M, N, K, varargin) 
% -------------------------------------------------------------------------
% 
% Summary of this function goes here
% Detailed explanation goes here


% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% -------------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('epsilon', nan);
    p.addOptional('SNR', 20);
    p.addOptional('alpha', []);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    W = rand(M, K);
    % W = W./sum(W, 1);
    % W = W.*dirichlet_rnd(ones(1, M), K);
    % W = 100*randn(M, K);
    % W(W<0) = 0;

    H = zeros(K, N);
    H(:, 1:K) = eye(K);
    if numel(options.alpha) == 0
        options.alpha = ones(1, K);
    end
    H(:, K+1:end) = dirichlet_rnd(options.alpha, N-K);
    Y = W*H;
    SNR = options.SNR; SNR = 10^(SNR/10);
    noise = randn(size(Y)); 
    sigma2 = sum(vecnorm(Y, 2, 1).^2) / M / N / SNR;
    noise = sqrt(sigma2)*noise;
    X = Y + noise;

    % indices = randperm(N);
    % X = X(:, indices);
    % H = H(:, indices);
    % r_pure_pixel_set = [];
    % pure_pixel_set = 1:K;
    % for i=1:numel(pure_pixel_set)
    %     r_pure_pixel_set(end+1) = find(indices == pure_pixel_set(i));
    % end
    % pure_pixel_set = r_pure_pixel_set;
    pure_pixel_set = 1:K;
end
