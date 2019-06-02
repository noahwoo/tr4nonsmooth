%
% trust region for non-smooth objective
% Input: 
% 
% Output:
% 

function [x obj iter eps info] = SmoothTR(obj_func, x0, opts)

% main loop for approximation smooth function
iter = 0;
info = 0;
stop = 0;

global gs_info;
% init trust region radius
delta_0 = 1;
while stop == 0
    % solve the optimization with mu = mu_k
    [x_k, obj_k, grad, delta_k, info] = tr_solve(obj_func, x0, delta_0, opts);
    
    % print the message 
    if opts.verbose >= 1
        fprintf('#Iteration-%d: objective = %f, norm(grad) = %f, radius = %f, mu = %f\n',...
            iter, obj_k, norm(grad), delta_k, opts.mu);
    end
    
    % set arguments for next smooth function
    % opts.mu_red = get_mu_red();
    opts.mu = opts.mu * opts.mu_red;
    x0 = x_k;
    % track status
    if info == 3 || info == 4
        gs_info.mu_fail = gs_info.mu_fail + 1;
    end
    delta_0 = delta_k;
    % break if tolerance met or maximal iteration number reached
    iter = iter + 1;
    if norm(grad) < opts.tol
        stop = 1;
        info = 1;
    elseif iter > opts.iter_max
        stop = 1;
        info = 3;
    end
    
    if stop == 1
        x = x_k; 
        obj = obj_k;
        eps = norm(grad);
        break;
    end
end
gs_info.mu_it = iter;

function coef = get_mu_red()
% get the way to do smoothing shrink
coef = 0.5;