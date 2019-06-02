% wrapper for SmoothTR
function [x obj iter eps info] = ... 
    SmoothTRWrap(H, f, mD, mtheta, mphi, msigma, mkiu, malpha, mlambda) 
% Input: 
%  H:      matrix in linear combine
%  f:      observed image vector
%  mD:     difference matrix
%  mtheta: linear or logistic method
%  mphi: loss function for the regularization term
%  msigma: arugment to balance the theta and phi
%  mkiu, malpha, mlambda: opts to phi
% Output: 

% call the SmoothTR method
% setup the arguments
opts = [];
opts.mu = 1e-2;
opts.mu_red = 0.01;
opts.iter_max = 100;
opts.tol = 1e-14;

% 0: par updated by iteration, 
% 1: par estimated by gradient and radius
opts.par_way = 0;
opts.verbose = 0;

% setup the global arguments 
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
% hessian matrix
global HM; 
% alpha for phi
global alpha;
% lambda for phi
global lambda;

% arguments
theta = mtheta;
sigma = msigma; % balance coefficient
kiu = mkiu; 
alpha = malpha; 
phi = mphi;
lambda = mlambda;

% data and coefficient
A = H; 
b = f;
D = mD; 
[m n] = size(A); 

% objective function
obj_func = [];

obj_func.obj = @objective; 
obj_func.grad = @gradient; 
obj_func.hessian = @hessian; 

if theta == 1
    obj_func.obj_theta = @obj_theta_lin;
    obj_func.grad_theta = @grad_theta_lin;
    obj_func.hessian_theta = @hessian_theta_lin;
elseif theta == 2
    obj_func.obj_theta = @obj_theta_lr;
    obj_func.grad_theta = @grad_theta_lr;
    obj_func.hessian_theta = @hessian_theta_lr;
end

if phi == 3
    obj_func.obj_phi = @obj_phi_mu;
    obj_func.grad_phi = @grad_phi_mu;
    obj_func.hessian_phi = @hessian_phi_mu;
elseif phi == 2
    obj_func.obj_phi = @obj_phi_mu_log;
    obj_func.grad_phi = @grad_phi_mu_log;
    obj_func.hessian_phi = @hessian_phi_mu_log;
elseif phi == 1
    obj_func.obj_phi = @obj_phi_mu_lr;
    obj_func.grad_phi = @grad_phi_mu_lr;
    obj_func.hessian_phi = @hessian_phi_mu_lr;    
elseif phi == 4
    obj_func.obj_phi = @obj_phi_mu_hard;
    obj_func.grad_phi = @grad_phi_mu_hard;
    obj_func.hessian_phi = @hessian_phi_mu_hard;    
elseif phi == 5
    obj_func.obj_phi = @obj_phi_mu_scad;
    obj_func.grad_phi = @grad_phi_mu_scad;
    obj_func.hessian_phi = @hessian_phi_mu_scad;    
elseif phi == 6
    obj_func.obj_phi = @obj_phi_mu_mcp;
    obj_func.grad_phi = @grad_phi_mu_mcp;
    obj_func.hessian_phi = @hessian_phi_mu_mcp;    
end


% if theta == 0
%   obj_func.obj  = @objective_lin;
%   obj_func.grad = @gradient_lin; 
%   obj_func.hessian = @hessian_lin;
% elseif theta == 1
%   obj_func.obj  = @objective_lr;
%   obj_func.grad = @gradient_lr; 
%   obj_func.hessian = @hessian_lr;
% elseif theta == 2
%   obj_func.obj  = @objective_lin_log;
%   obj_func.grad = @gradient_lin_log; 
%   obj_func.hessian = @hessian_lin_log;  
% elseif theta == 3
%   obj_func.obj  = @objective_lr_log;
%   obj_func.grad = @gradient_lr_log; 
%   obj_func.hessian = @hessian_lr_log;    
% elseif theta == 4
%   obj_func.obj  = @objective_lin_lr;
%   obj_func.grad = @gradient_lin_lr; 
%   obj_func.hessian = @hessian_lin_lr;      
% elseif theta == 5
%   obj_func.obj  = @objective_lr_lr;
%   obj_func.grad = @gradient_lr_lr; 
%   obj_func.hessian = @hessian_lr_lr;
% end

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

% solve the non-smooth optimization problem
x0 = zeros(n,1);
% x0 = [1; 5; 0]; 
% x0 = [0.545; 0.237; 0; 0.098; 0.165; 0; 0; 0.059];
% f = objective_lin(x0, 0) 

[x obj iter eps info] = SmoothTR(obj_func, x0, opts);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% S
function y = S_mu(t, mu)
y = sqrt(t.*t+4*mu^2);
function g = S_mu_g(t, mu)
g = t./sqrt(t.*t+4*mu^2);
function h = S_mu_h(t, mu) 
h = 4*mu^2 ./ (t.^2+4*mu^2).^(1.5);

% S_hat
function y = S_h_mu(t, mu)
y = 0.5*(t + S_mu(t, mu));
function g = S_h_mu_g(t, mu) 
g = 0.5*(1+S_mu_g(t, mu));
function h = S_h_mu_h(t, mu) 
h = 0.5*S_mu_h(t, mu);

% S_check
function y = S_c_mu(t, mu)
y = 1 - S_h_mu(1-t, mu);
function g = S_c_mu_g(t, mu) 
g = S_h_mu_g(1-t, mu);

% linear regression theta functions
function f = obj_theta_lin(x) 
global A;
global b;
[m n] = size(A);
f = sum((A*x - b).^2); 

function g = grad_theta_lin(x)
global A;
global b;
g = 2*(A')*(A*x-b);

function hessian_theta_lin(x) 
global A;
global HM;
[m n] = size(A);
HM = HM + 2*(A')*A;

% logistic regression theta
function f = obj_theta_lr(x) 
global A;
global b;
[m n] = size(A);
Ax = A*x;
f = -sum(b.*(Ax) - log(exp(Ax) + 1));

function g = grad_theta_lr(x)
global A;
global b;

[m n] = size(A);
ea = exp(A*x);
g = - A'*(b - ea./(1+ea));

function H = hessian_theta_lr(x) 
global A;
global b;
global HM; 
[m n] = size(A);
dd = zeros(m,1);
ea = exp(A*x);
dd = ea./((1+ea).*(1+ea));
HM = HM + (A')*diag(dd)*A;

%% phi functions t^q
% 
function f = obj_phi_mu( x, mu ) 
global D;
global kiu; 

[m n] = size(D);
f = sum(((D*x).^2+4*mu^2).^(kiu/2));

function g = grad_phi_mu( x, mu ) 
global D;
global kiu; 

[m n] = size(D);
dx = D*x;
dx = (dx.^2+4*mu^2).^(kiu/2-1)*kiu.*dx;
g = (D')*dx;

function hessian_phi_mu( x, mu, sigma )
global D;
global kiu; 
global HM; 

dx = D*x;
dd = kiu*(dx.^2+4*mu^2).^(0.5*kiu-1) + ...
    (dx.^2+4*mu^2).^(0.5*kiu-2)*kiu*(kiu-2).*dx.^2;
HM = HM + sigma * (D')*diag(dd)*D;

%% phi function, log(alpha * t^q + 1)
% 
function f = obj_phi_mu_log(x, mu) 
global D;
global kiu; 
global alpha; 

f = sum(log(1+alpha*((D*x).^2+4*mu^2).^(kiu/2)));

function g = grad_phi_mu_log(x, mu)
global D;
global kiu; 
global alpha;

dx = D*x;
dx = alpha * kiu * dx .* (dx.*dx+4*mu^2).^(kiu/2-1) ./ ... 
    (1+alpha*(dx.*dx+4*mu^2).^(kiu/2));
g  = (D')*dx;

function hessian_phi_mu_log( x, mu, sigma )
global D;
global kiu;
global alpha;
global HM;

dx = D*x;
t  = (dx.*dx+4*mu^2).^0.5;
t1 = dx.*(dx.*dx+4*mu^2).^(-0.5);
t2 = (dx.*dx+4*mu^2).^(-0.5) - dx.*dx.*(dx.*dx+4*mu^2).^(-1.5);
at = (alpha*t.^kiu+1);

dd = kiu*alpha*((kiu-1)*t.^(kiu-2).*t1.^2 + t.^(kiu-1).*t2).*at;
dd = dd - (kiu*alpha*t.^(kiu-1).*t1).^2;
dd = dd./at.^2;

HM = HM + sigma * (D')*diag(dd)*D;

%% phi function, alpha * t^q / alpha * t^q + 1
% 
function f = obj_phi_mu_lr(x, mu) 
global D;
global kiu; 
global alpha; 

aq = alpha*((D*x).^2+4*mu^2).^(kiu/2);
f = sum(aq./(1+aq));

function g = grad_phi_mu_lr(x, mu)
global D;
global kiu; 
global alpha;

dx = D*x;
dx = alpha * kiu * dx .* (dx.*dx+4*mu^2).^(kiu/2-1) ./ ... 
    (1+alpha*(dx.*dx+4*mu^2).^kiu);
g  = (D')*dx;

function hessian_phi_mu_lr( x, mu, sigma )
global D;
global kiu;
global alpha;
global HM;

dx = D*x;
t  = (dx.*dx+4*mu^2).^0.5;
t1 = dx.*(dx.*dx+4*mu^2).^(-0.5);
t2 = (dx.*dx+4*mu^2).^(-0.5) - dx.*dx.*(dx.*dx+4*mu^2).^(-1.5);
at = (alpha*t.^kiu+1);

dd = alpha*kiu*((kiu-1)*t.^(kiu-2).*t1.^2 + t.^(kiu-1).*t2).*at;
dd = dd - 2*(alpha*kiu*t.^(kiu-1).*t1).^2;
dd = dd./at.^3;

HM = HM + sigma * (D')*diag(dd)*D;

%% phi function, lambda - S^h_mu^2(lambda - S_mu(t))/lambda 
% 
function f = obj_phi_mu_hard(x, mu)
global D;
global lambda;

t = D*x;
t1 = S_h_mu(lambda - S_mu(t, mu), mu);
f = sum( lambda - t1.*t1/lambda );

function g = grad_phi_mu_hard(x, mu) 
global D;
global lambda; 

t = D*x;
t1 = lambda - S_mu(t, mu); 
g = (2/lambda) * S_h_mu(t1, mu) .* S_h_mu_g(t1, mu) .* S_mu_g(t, mu);
g = (D')*g; 

function hessian_phi_mu_hard(x, mu, sigma) 
global D;
global lambda;
global HM;

t = D*x;
t1 = lambda-S_mu(t, mu);
t2 = S_mu_g(t, mu);
t3 = S_h_mu(t1, mu);
t4 = S_h_mu_g(t1, mu);
t5 = S_h_mu_h(t1, mu);

dd = -(2/lambda)*(t4 .* t2).^2;
dd = dd -(2/lambda) * t3 .* t5 .* t2.^2;
dd = dd +(2/lambda) * t3 .* t4 .* S_mu_h(t, mu);

HM = HM + sigma * (D')*diag(dd)*D;

%% phi function, scad
% 
function y = scad_obj_mu(t, mu) 
global lambda;
global alpha;

y = S_c_mu(S_h_mu(alpha - t./lambda, mu)./(alpha-1), mu);

function f = obj_phi_mu_scad(z, mu) 
global D;

t = D*z;
f = 0;
for k=1:length(t)
    tt = t(k);
    f = f + quad(@(x)scad_obj_mu(x, mu), 0, S_mu(tt, mu));
end

function g = grad_phi_mu_scad(z, mu)
global D;

t = D*z;
g = scad_obj_mu(S_mu(t, mu), mu) .* S_mu_g(t, mu);
g = (D')*g;

function hessian_phi_mu_scad(z, mu, sigma) 
global D;
global lambda; 
global HM;
global alpha 

t = D*z;
t1 = alpha - S_mu(t, mu)/lambda; 
t2 = S_h_mu(t1, mu); 
t3 = t2/(alpha-1); 

dd = -1/(lambda*(alpha-1)) * S_c_mu_g(t3, mu).*S_h_mu_g(t1, mu).*S_mu_g(t,mu).^2;
dd = dd + S_c_mu(t3, mu) .* S_mu_h(t, mu); 

HM = HM + sigma * (D')*diag(dd)*D;

%% phi function: Minimax concave plus
% 
function y = mcp_obj_mu(t, mu)
global lambda;
global alpha;

y = S_h_mu(1 - t./(alpha*lambda), mu);

function f = obj_phi_mu_mcp(z, mu) 
global D;
global alpha;
global lambda;

t = D*z;
smuT = S_mu(t, mu);
aT   = 1-smuT./(alpha*lambda);
aT2  = aT.*aT;
f = sum( 0.5*( 1 - aT2) + ... 
    0.5*(sqrt(1+4*mu^2) - aT.*sqrt(aT2+4*mu^2)) + ...
    2*mu^2*(log(1+sqrt(1+4*mu^2))-log(aT+sqrt(aT2+4*mu^2))));
f = 0.5*alpha*lambda*f;

% g = 0;
% for k=1:length(t)
%     tt = t(k);
%     g = g + quad(@(x)mcp_obj_mu(x, mu), 0, S_mu(tt, mu));
% end

function g = grad_phi_mu_mcp(z, mu)
global D; 
global alpha;
global lambda;

t = D*z; 
g = S_h_mu(1 - S_mu(t, mu)./(alpha*lambda), mu) .* S_mu_g(t, mu); 
g = (D')*g;

function hessian_phi_mu_mcp(z, mu, sigma)
global D;
global HM;
global alpha; 
global lambda;

t  = D*z;
t1 = 1 - S_mu(t, mu)./(alpha*lambda); 
dd = -1/(alpha*lambda) .* S_c_mu_g(t1, mu) .* S_mu_g(t, mu).^2 + ...
    S_c_mu(t1, mu) .* S_mu_h(t, mu);
HM = HM + sigma * (D')*diag(dd)*D;

%%%%%%%%%%%%%%%%%%%%%%%%%%
% objective, gradient and hessian
% linear regression
% function f = objective_lin( x, mu )
% global sigma;
% % mu = 0;
% % x = [3.2429; 5.4641];
% f = obj_theta_lin(x) + sigma * obj_phi_mu(x, mu);
% 
% function g = gradient_lin( x, mu )
% global sigma; 
% g = grad_theta_lin(x) + sigma * grad_phi_mu(x, mu);
% 
% function hessian_lin( x, mu )
% global sigma;
% global HM;
% n = length(x);
% HM = zeros(n, n);
% hessian_theta_lin(x); 
% hessian_phi_mu(x, mu, sigma);

%%%
% function f = objective_lin_log( x, mu )
% global sigma;
% % mu = 0;
% % x = [3.2429; 5.4641];
% f = obj_theta_lin(x) + sigma * obj_phi_mu_log(x, mu);
% 
% function g = gradient_lin_log( x, mu )
% global sigma; 
% g = grad_theta_lin(x) + sigma * grad_phi_mu_log(x, mu);
% 
% function hessian_lin_log( x, mu )
% global sigma;
% global HM;
% n = length(x);
% HM = zeros(n, n);
% hessian_theta_lin(x); 
% hessian_phi_mu_log(x, mu, sigma);

%%%
% function f = objective_lin_lr( x, mu )
% global sigma;
% % mu = 0;
% % x = [3.2429; 5.4641];
% f = obj_theta_lin(x) + sigma * obj_phi_mu_lr(x, mu);
% 
% function g = gradient_lin_lr( x, mu )
% global sigma; 
% g = grad_theta_lin(x) + sigma * grad_phi_mu_lr(x, mu);
% 
% function hessian_lin_lr( x, mu )
% global sigma;
% global HM;
% n = length(x);
% HM = zeros(n, n);
% hessian_theta_lin(x); 
% hessian_phi_mu_lr(x, mu, sigma);

%%%
% objective, gradient and hessian
% logistic regression
% function f = objective_lr( x, mu )
% global sigma;
% f = obj_theta_lr(x) + sigma * obj_phi_mu(x, mu);
% 
% function g = gradient_lr( x, mu )
% global sigma; 
% g = grad_theta_lr(x) + sigma * grad_phi_mu(x, mu);
% 
% function H = hessian_lr( x, mu )
% global sigma;
% global HM;
% 
% n = length(x);
% HM = zeros(n, n);
% 
% hessian_theta_lr(x); 
% hessian_phi_mu(x, mu, sigma);

%%%
% function f = objective_lr_log( x, mu )
% global sigma;
% f = obj_theta_lr(x) + sigma * obj_phi_mu_log(x, mu);
% 
% function g = gradient_lr_log( x, mu )
% global sigma; 
% g = grad_theta_lr(x) + sigma * grad_phi_mu_log(x, mu);
% 
% function hessian_lr_log( x, mu )
% global sigma;
% global HM;
% n = length(x);
% HM = zeros(n, n);
% hessian_theta_lr(x); 
% hessian_phi_mu_log(x, mu, sigma);

%%%
% function f = objective_lr_lr( x, mu )
% global sigma;
% f = obj_theta_lr(x) + sigma * obj_phi_mu_lr(x, mu);
% 
% function g = gradient_lr_lr( x, mu )
% global sigma; 
% g = grad_theta_lr(x) + sigma * grad_phi_mu_lr(x, mu);
% 
% function H = hessian_lr_lr( x, mu )
% global sigma;
% global HM;
% 
% n = length(x);
% HM = zeros(n, n);
% 
% hessian_theta_lr(x); 
% hessian_phi_mu_lr(x, mu, sigma);

% global objective
function f = objective(x, mu, theta_f, phi_f)
global sigma; 

f = theta_f(x) + sigma*phi_f(x, mu);

% global gradient
function g = gradient(x, mu, theta_g, phi_g)
global sigma;

g = theta_g(x) + sigma*phi_g(x, mu);

% global hessian
function hessian(x, mu, theta_h, phi_h)
global sigma;
global HM;

n = length(x);
HM = zeros(n, n);
theta_h(x); 
phi_h(x, mu, sigma);