% experiment for the  prostate cancer data
function prostate
clear all;
path(path, 'gqtpar');
verbose = 0;

% 
% global variables for solving
% 
% coefficient to balance the trade-off between theta and phi
global sigma;
% q in |t|^q
global kiu;
% each row of A corresponds to an instance of data
global A;
% y-value of each instance
global b; 
% hessian matrix
global HM;

% load data
dat = load('dataset/prostate.data.tt');
% data scaling
dat(:, 9) = dat(:, 9) - mean(dat(:, 9));
for j=1:8
    dat(:, j) = (dat(:, j) - mean(dat(:, j)))/std(dat(:, j));
end
% data spliting
T = dat(1:67,:);
F = dat(68:97,:);

A = T(:,1:8); 
b = T(:,9);
D = eye(8);

B = F(:,1:8); 
c = F(:,9);

% 
% options
%
opts = [];
opts.mu = 1e-2;
opts.mu_red = 0.01;
opts.iter_max = 100;
opts.tol = 1e-13;

% 0: par updated by iteration, 1: par estimated by gradient and radius
opts.par_way = 0;
opts.verbose = 1;

method = 0; % 0: linear loss, 1: log-linear loss
% objective function
obj_func = [];
if method == 0
    obj_func.obj  = @objective_lin;
    obj_func.grad = @gradient_lin;
    obj_func.hessian = @hessian_lin;
elseif method == 1
    obj_func.obj  = @objective_loglin;
    obj_func.grad = @gradient_loglin;
    obj_func.hessian = @hessian_loglin;
end

% stats about the alg.
global gs_info;
gs_info = [];
gs_info.mu_it = 0;
gs_info.mu_fail = 0;
gs_info.tr_it = 0;
gs_info.tr_fail = 0;
gs_info.gqt_it = 0;
gs_info.obj_eval = 0;
gs_info.grad_eval = 0;
gs_info.hessian_eval = 0;

% initial point 
x0 = zeros(8,1);
% x0 = [0.655; 0.208; 1.24e-7; 6.89e-7; 0.232; 7.42e-13; 3.67e-12; 1.71e-8];

% model arguments
kiu = 0.2;
sigmaa = [7.5];

% search for the maximal sigma
for k=1:length(sigmaa)
    sigma = sigmaa(k);
    % f = obj_func.obj(x0, 1e-5)
    % g = obj_func.grad(x0, 1e-5)
    [x obj iter eps info] = SmoothTR(obj_func, x0, opts);
    if verbose == 1
        fprintf('objective=%.8f, #iterate=%d, tolerance=%.8f\n', ...
                    obj, iter, eps);
    end
    for j=1:8
        fprintf('x(%d) = %.16f\n', j, x(j));
    end
    
    
    % %
    % adjust the optimal solution
    L = (sigma * kiu * (1-kiu) ./(2*sum(A.*A, 1))).^(1/(2-kiu));
    L = L';
    Lg = ((sigma * kiu) / (2 * norm(A,2) * sqrt(objective_lin(x,0))))^(1/(1-kiu));
    Lg = Lg * ones(length(x), 1);
    
    I = find(L < Lg);
    L(I) = Lg(I);
    L 
    I = find(abs(x) <= L);
    if length(I) > 0
        x(I) = 0;
    end
    % end-adjust
    % %
    
    % %
    % check the optimality
    % first order
    scaled_grad = 2*diag(x)*A'*(A*x-b) + sigma * kiu * abs(x).^kiu;
    nscaled_grad = norm(scaled_grad)
    % second order
    I = find(x > 0);
    H = 2*diag(x)*A'*A*diag(x) + sigma*kiu*(1-kiu)*diag(abs(x).^sigma);
    Hh = H(I, I);
    E = eig(Hh);
    I = find(E <= 0);
    if length(I) == 0
        fprintf('second order condition satisfied\n');
    end
    % end-check
    % %
        
    % x = [0.655; 0.208; 1.24e-7; 6.89e-7; 0.232; 7.42e-13; 3.67e-12; 1.71e-8];
    % x = [0.545; 0.237; 0; 0.098; 0.165; 0; 0; 0.059];
    v = B*x;
    MSE = mean((v-c).^2);
    fprintf('MSE=%f, lambda=%f, q=%f\n', MSE, sigma, kiu);
    fprintf('Solving stats:\n');
    fprintf('  #smooth: %d, #smooth_fail: %d\n', gs_info.mu_it, gs_info.mu_fail);
    fprintf('  #tr: %d, #tr_fail: %d\n', gs_info.tr_it, gs_info.tr_fail);
    fprintf('  #gqtpar: %d\n', gs_info.gqt_it);
    fprintf('  #obj_eval: %d, #grad_eval: %d, #hessian_eval: %d\n', ...
        gs_info.obj_eval, gs_info.grad_eval, gs_info.hessian_eval);
end

%
% objective, gradient and hessian for linear and log-linear regression
% 

% linear regression theta functions
% objective
function f = objective_lin( x, mu )

global sigma;
f = obj_theta_lin(x) + sigma * obj_phi_mu(x, mu);

% gradient
function g = gradient_lin( x, mu )

global sigma; 
g = grad_theta_lin(x) + sigma * grad_phi_mu(x, mu);

% hessian
function hessian_lin( x, mu )

global sigma;
global HM;

n = length(x);
HM = zeros(n, n);
hessian_theta_lin(x); 
hessian_phi_mu(x, mu, sigma);

% linear objective
function f = obj_theta_lin(x) 

global A;
global b;

[m n] = size(A);
f = sum((A*x - b).^2); 

% linear gradient
function g = grad_theta_lin(x)

global A;
global b;

g = 2*(A')*(A*x-b);

% linear hessian
function hessian_theta_lin(x) 

global A;
global HM;

[m n] = size(A);
HM = HM + 2*(A')*A;

% log-linear regression theta function
% objective
function f = objective_loglin( x, mu )

global sigma;
f = obj_theta_loglin(x) + sigma * obj_phi_mu(x, mu);

% gradient
function g = gradient_loglin( x, mu )

global sigma; 
g = grad_theta_loglin(x) + sigma * grad_phi_mu(x, mu);

% hessian
function hessian_loglin( x, mu )

global sigma;
global HM;

n = length(x);
HM = zeros(n, n);
hessian_theta_lin(x); 
hessian_phi_mu(x, mu, sigma);

% loglin objective
function f = obj_theta_loglin(x) 

global A;
global b;

[m n] = size(A);
f = log(sum((A*x - b).^2)+1); 

% loglin gradient
function g = grad_theta_loglin(x)

global A;
global b;
Ab = (A*x-b);
g = 2*(A')*Ab / (sum(Ab.^2)+1);

% loglin hessian
function hessian_theta_loglin(x) 

global A;
global HM;

Ab = (A*x - b); 

[m n] = size(A);
HM = HM + 2*(A')*A - 2* (A'*Ab)*(Ab'*A) / (sum(Ab.^2)+1)^2;

% phi functions
function f = obj_phi_mu( x, mu ) 

global kiu; 

f = sum((x.^2+mu^2).^(kiu/2));

% gradient
function g = grad_phi_mu( x, mu ) 

global kiu; 

g = (x.^2+mu^2).^(kiu/2-1)*kiu.*x;

% hessian
function hessian_phi_mu( x, mu, sigma )

global kiu; 
global HM; 

dd = kiu*(x.^2+mu^2).^(0.5*kiu-1) + ...
    (x.^2+mu^2).^(0.5*kiu-2)*kiu*(kiu-2).*x.^2;
HM = HM + sigma * diag(dd);