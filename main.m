%
% entry point of trust region method for non-smooth objective function 
% 
function main

path(path, 'gqtpar');

IM = imread('SheppLoganPhantom.png');
imshow(IM);

%% call the SmoothTR method

% setup the arguments
opts = [];
opts.mu = 1;
opts.mu_red = 0.5;
opts.iter_max = 50;
opts.tol = 1e-8;

% 0: par updated by iteration, 
% 1: par estimated by gradient and radius
opts.par_way = 0;
opts.verbose = 0;

%% setup the global arguments 
% each column of D corresponds to a combination 
% of free variables in sparse constraints
global D; 
% coefficient to balance the trade-off between theta and phi
global sigma;
% q in |t|^q
global kiu; 
% each row of A corresponds to a instance of data
global A; 
% y-value of each instance
global b; 
% theta functional
global theta;

n = 2;

% theta functions: 0: Linear Regression, 1: Logistic Regression
theta = 0;
%% A, b for linear regression
% A = [0.3 1; 1 0.1; 2 0.5];
% b = [2; 1; 3];
% D = [ 1 1; 1 0 ];
% kiu = 1;

%% sample estimation
% (g1, g2) for estimation
g1 = 0.1; g2 = 0.8;
A = [1 0; 0 1];
b = [0.1; 0.8];
D = [1; -1];
kiu = 0.5;

%% A, b for logistic regression
% A = [-0.2, -0.5; 2 1; 0.4 5];
% b = [-1; 1; 1];
% D = [ 1 1; 1 0 ];
% kiu = 1;

% plot the graph
% x = -1:0.1:1;
% y = -1:0.1:1; 
% len = length(x); 
% Z = zeros(len, len);
% for p=1:len
%    for q=1:len
%        v = [x(p); y(q)];
%        Z(p,q) = obj_theta_lr(v);
%    end
% end
% mesh(x,y,Z);

% balance coefficient
sigma = 0.2;

%% objective function
obj_func = [];
if theta == 0
  obj_func.obj  = @objective_lin;
  obj_func.grad = @gradient_lin; 
  obj_func.hessian = @hessian_lin;
elseif theta == 1
  obj_func.obj  = @objective_lr;
  obj_func.grad = @gradient_lr; 
  obj_func.hessian = @hessian_lr; 
end
    
%% stats about the alg.
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

%% solve the non-smooth optimization problem
x0 = zeros(n,1);
[x obj iter eps info] = SmoothTR(obj_func, x0, opts);
if info == 1 || info == 2
    fprintf('Optimal point found at %d iteration with: \n  tolerance %f\n  objective %f.\n',...
        iter, eps, obj);
    x
else
    fprintf('Failed to get optimal point\n');
    x
end

fprintf('Solving stats:\n');
fprintf('  #smooth: %d, #smooth_fail: %d\n', gs_info.mu_it, gs_info.mu_fail);
fprintf('  #tr: %d, #tr_fail: %d\n', gs_info.tr_it, gs_info.tr_fail);
fprintf('  #gqtpar: %d\n', gs_info.gqt_it);
fprintf('  #obj_eval: %d, #grad_eval: %d, #hessian_eval: %d\n', ...
    gs_info.obj_eval, gs_info.grad_eval, gs_info.hessian_eval);

% objective, gradient and hessian
% linear regression
function f = objective_lin( x, mu )
global sigma;
f = obj_theta_lin(x) + sigma * obj_phi_mu(x, mu);

function g = gradient_lin( x, mu )
global sigma; 
g = grad_theta_lin(x) + sigma * grad_phi_mu(x, mu);

function H = hessian_lin( x, mu )
global sigma;
H = hessian_theta_lin(x) + sigma * hessian_phi_mu(x, mu);

% linear regression theta functions
function f = obj_theta_lin(x) 
global A;
global b;
[m n] = size(A);
f = 0;

for k=1:m
    f = f + (A(k,:)*x - b(k))^2;
end

function g = grad_theta_lin(x)
global A;
global b;
[m n] = size(A);
g = zeros(n, 1);
for k=1:m
    g = g + 2*(A(k,:)*x - b(k))*A(k,:)';
end

function H = hessian_theta_lin(x) 
global A;
[m n] = size(A);
H = 2*(A')*A;

% logistic regression
function f = objective_lr( x, mu )
global sigma;
f = obj_theta_lr(x) + sigma * obj_phi_mu(x, mu);

function g = gradient_lr( x, mu )
global sigma; 
g = grad_theta_lr(x) + sigma * grad_phi_mu(x, mu);

function H = hessian_lr( x, mu )
global sigma;
H = hessian_theta_lr(x) + sigma * hessian_phi_mu(x, mu);

% logistic regression theta
function f = obj_theta_lr(x) 
global A;
global b;
[m n] = size(A);
f = 0;
for k=1:m
    f = f + log(1+exp(-b(k)*A(k,:)*x));
end

function g = grad_theta_lr(x)
global A;
global b;
[m n] = size(A);
g = zeros(n, 1);
for k=1:m
    eba = exp(-b(k)*A(k,:)*x);
    g = g - b(k)*(eba/(1+eba))*A(k,:)';
end

function H = hessian_theta_lr(x) 
global A;
global b;
[m n] = size(A);
dd = zeros(m,1);
for k=1:m
    eba = exp(-b(k)*A(k,:)*x);
    t = eba/(1+eba);
    dd(k) = t*(1-t);
end
H = (A')*diag(dd)*A;

% phi functions
function f = obj_phi_mu( x, mu ) 
global D;
global kiu; 

[n m] = size(D);
f = 0; 
for k = 1:m
    f = f+((D(:,k)'*x)^2 + mu)^(0.5*kiu);
end

function g = grad_phi_mu( x, mu ) 
global D;
global kiu; 

[n m] = size(D);
g = zeros(n, 1);
for k = 1:m
    dx = D(:,k)'*x;
    g = g + ((dx).^2 + mu)^(0.5*kiu-1) * kiu * dx * D(:,k);
end

function H = hessian_phi_mu( x, mu )
global D;
global kiu; 

[n m] = size(D);
dd = zeros(m,1);
for k = 1:m
    dx = D(:,k)'*x;
    d = kiu*(dx^2+mu)^(0.5*kiu-1) + ... 
        (dx^2+mu)^(0.5*kiu-2)*kiu*(kiu-2)*dx^2;
    dd(k) = d; 
end
H = D*diag(dd)*(D');