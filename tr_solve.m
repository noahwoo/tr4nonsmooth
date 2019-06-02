% 
% solve the smooth objective with trust-region method
% 
function [x obj grad delta_k info] = tr_solve(obj_func, x0, delta_0, opts) 

% set the intial x-value and trust region radius
j = 0;
x_j = x0;
delta_j = delta_0;

% arguments for gqtpar
rtol = 0.1;
atol = 0;
maxit = 12;
par = 0;

% arguments for trust region
redf_eps = 1e-32;
bnd_eps  = 1e-10;
eta1 = 0.1;
eta2 = 0.9;
gama1 = 0.5;
gama2 = 1.5;
delta_max = 1e12;
delta_min = 1e-4; % OR setting this as 1e-1
step_ratio = 0.8; 
grad_red = 1;

global gs_info;

global HM;
% solve by trust region method
obj_j = obj_func.obj(x_j, opts.mu, obj_func.obj_theta, obj_func.obj_phi);
grad_j = obj_func.grad(x_j, opts.mu, obj_func.grad_theta, obj_func.grad_phi);
obj_func.hessian(x_j, opts.mu, obj_func.hessian_theta, obj_func.hessian_phi);

p_grad = grad_j;

% track the status
gs_info.obj_eval = gs_info.obj_eval + 1;
gs_info.grad_eval = gs_info.grad_eval + 1;
gs_info.hessian_eval = gs_info.hessian_eval + 1;

stop = 0;
info = 0;

% global g_proc_info; 
% global g_nit;
% global g_fid;
% % 
% % % coefficient to balance the trade-off between theta and phi
% global sigma;
% % q in |t|^q
% global kiu; 
% % each row of A corresponds to a instance of data
% global A; 
% % y-value of each instance
% global b; 
max_oit = 100;

x = x_j;
obj = obj_j;
grad = grad_j;
delta_k = 0;

while stop == 0
  % construct the quadratic form 
  bgrad = grad_j; % gradient 
  % H = hessian_j; % hession 
  if opts.par_way == 1
      par = norm(bgrad)/delta_j;
  end
  % call quadratic solver 
  [p, redf, par, nit, z, info] = ... 
    gqtparg(HM, bgrad, delta_j, rtol, atol, maxit, par);
  if opts.verbose >= 3 && (info == 3 || info == 4)
      % add the stat info
      fprintf('WARNING: quadratic solver terminated at INFO=%d\n', info);
  end
  % track the status
  gs_info.tr_it = gs_info.tr_it + 1;
  if info == 3 || info == 4
    gs_info.tr_fail = gs_info.tr_fail + 1;
  end
  gs_info.gqt_it = gs_info.gqt_it + nit;
  
  % check objective descent
  if redf <= 0 % redf_eps % TODO: check the relationship with rtol and atol
      x = x_j;
      obj = obj_j;
      grad = grad_j;
      delta_k = delta_j;
      
      if opts.verbose >= 3
          fprintf('  Subproblem terminated at iteration %d, objective= %f, redf = %f, norm(p) = %f, norm(grad) = %f info = %d\n',... 
              j, obj, redf, norm(p), norm(grad_j), info);
          % fprintf('  x(1)=%.10f, x(2)=%.10f\n', x(1), x(2));
      end
      info = 1;
      break;
  end
  obj_new = obj_func.obj(x_j+p, opts.mu, obj_func.obj_theta, obj_func.obj_phi);
  rho = (obj_j - obj_new) / redf;
  gs_info.obj_eval = gs_info.obj_eval + 1;
  
  % move to the next point
  if rho > eta1 
      x_j = x_j + p;
      p_grad = grad_j;
      % update gradient and objective 
      grad_j = obj_func.grad(x_j, opts.mu, obj_func.grad_theta, obj_func.grad_phi);
      obj_func.hessian(x_j, opts.mu, obj_func.hessian_theta, obj_func.hessian_phi);
      % track the status
      gs_info.grad_eval = gs_info.grad_eval + 1;
      gs_info.hessian_eval = gs_info.hessian_eval + 1;
      obj_j = obj_new;
  end
  
  % % record solution info
%   fprintf(g_fid, 'iteration=%d, x       = (%.8f, %.8f, %.8f)\n',...
%       g_nit, x_j(1), x_j(2), x_j(3));
%   grad_x = 2*diag(x_j) * A'*(A*x_j-b) + sigma*kiu*abs(x_j).^kiu; 
%   
%   fprintf(g_fid, 'iteration=%d, obj     = %.8f\n', g_nit, obj_func.obj(x_j,0));
%   fprintf(g_fid, 'iteration=%d, grad    = (%.16f, %.16f, %.16f)\n', ...
%       g_nit, grad_x(1), grad_x(2), grad_x(3));
%   
%   HM_x = 2*diag(x_j) * A' * A * diag(x_j) + sigma * kiu * (kiu-1) * diag(abs(x_j).^kiu);
%   
%   fprintf(g_fid, '                        %.8f %.8f %.8f\n', ...
%              HM_x(1,1), HM_x(1,2), HM_x(1,3));
%   fprintf(g_fid, 'iteration=%d, hessian = %.8f %.8f %.8f\n', ... 
%       g_nit, HM_x(2,1), HM_x(2,2), HM_x(2,3));
%   fprintf(g_fid, '                        %.8f %.8f %.8f\n', ...
%              HM_x(3,1), HM_x(3,2), HM_x(3,3));
%   E = eig(HM_x(1:2, 1:2)); 
%   fprintf(g_fid, 'iteration=%d, eigen   = (%.8f %.8f)\n', ...
%       g_nit, E(1), E(2));
%   fprintf(g_fid, 'iteration=%d, mu      = %.16f \n\n', g_nit, opts.mu);
  
%   g_nit = g_nit + 1;
  
  % obj_func.hessian(x_j, 0)
  % update trust region radius
  if rho > eta2 && norm(p) >= step_ratio * delta_j
      delta_j = max(delta_min, min([gama2 * delta_j, delta_max]));
  elseif rho < eta1 % shrink
      delta_j = gama1 * delta_j;
  end
  
  % verbose
  if opts.verbose >= 4
      fprintf('  #Iteration-%d: objective=%f, radius=%f, lambda=%f\n', ... 
          j, obj_j, delta_j, par);
  end
  
  % check gradient
  if (norm(p_grad) < grad_red*opts.mu && delta_j >= delta_min)
      stop = 1;
      info = 2;
  elseif delta_j < bnd_eps
      stop = 1;
      info = 3;
  end
  
  % check maximal iterations
  if j > max_oit
      stop = 1;
      info = 3;
  end
  
  if stop == 1
      x = x_j;
      obj = obj_j;
      grad = grad_j;
      delta_k = delta_j;
      % verbose
      if opts.verbose >= 3
          fprintf('  Subproblem terminated at iteration %d, objective= %f, info = %d\n',... 
              j, obj, info);
      end
      break;
  end
  % iter for the next
  j = j + 1;
end