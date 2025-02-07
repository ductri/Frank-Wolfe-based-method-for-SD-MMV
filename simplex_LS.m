function [out] = simplex_LS(A, b, varargin) 
% ----------------------------------------------------------------------
% 
% min_{x}    ||Ax - b||_F^2 
% subject to 1^T x = 1, x >= 0

% Author: Tri Nguyen (nguyetr9@oregonstate.edu)

% ----------------------------------------------------------------------
    
    % Check the options structure.
    p = inputParser;
    p.addOptional('verbose', 0);
    p.addOptional('debug', 0);
    p.addOptional('max_iters', 0);
    p.KeepUnmatched = true;
    p.parse(varargin{:});
    options = p.Results;

    ops = struct;
    ops.verbose = false;
    ops.debug = false;
    ops.max_iters = 300;
    ops.tol = 1e-10;
    p_fn = @(x, l) proj_simplex(x);
    step_size = 1/svds(A, 1)^2;

    init_point = rand(size(A, 2), 1);
    g_fn = @(x) A'*(A*x - b);
    ops.f_fn = @(x) norm(A*x - b)^2;
    ops.ground_truth = ones(size(init_point));
    [x, tracking] = pgd_fista(g_fn, p_fn, step_size, init_point, ops);
    out = norm(A*x - b);
end

