% experiment 
function experiment

clear all;

path(path, 'gqtpar');

verbose = 1;
% exps = [1; 2; 3; 4; 5];
exps = [ 3 ];
lex  = length(exps); 

for k=1:lex
    if exps(k) == 1
        % 1. sythetic circle image data
        n = 32;
        figure;
        I = gen_img(n);
        R = 0.2*randn(n); 
        I = I + R;
        imshow(I);

        I = reshape(I, n*n, 1);
        D = gen_diff_mat(n);
        H = eye(n*n);

        mtheta = 1; msigma = 0.2; mkiu = 0.5;
        [x obj iter eps info] = ... 
            SmoothTRWrap(H, I, D, mtheta, msigma, mkiu); 

        if verbose > 0
            fprintf('objective=%.8f, #iterate=%d, tolerance=%.8f\n', ...
            obj, iter, eps);
        end

        figure;
        RI = reshape(x, n, n);
        imshow(RI);
    elseif exps(k) == 2 
        % 2. Shepp-Logan image skipped
    elseif exps(k) == 3 
        % 3. prostate cancer
        dat = load('dataset/prostate.data.tt');
        [m n] = size(dat);
        beta = [0.740; 0.367; 0; 0; 0; 0; 0; 0];
        % data scaling
        dat(:, 9) = dat(:, 9) - mean(dat(:, 9));
        for j=1:8
            dat(:, j) = (dat(:, j) - mean(dat(:, j)))/std(dat(:, j));
        end
        
        % data spliting
        T = dat(1:67,:); 
        F = dat(68:97,:);
        A = T(:,1:9-1); b = T(:,9);
        [m n] = size(A);
        D = eye(n);
        
        mtheta = 1; mphi = 3; 
        mkiu = 0.5; malpha = 3.7;
        msigmaa = [8.0]; % [5:0.5:20.0];% [8:1:8]; %5:0.5:20.0];
        
        maxMSE = 100;
        maxSigma = 0;
        maxAlpha = 0;
        maxX = zeros(8,1);
        % search for the maximal sigma
        for malpha = 2.6:0.2:3.4 % 2.3:0.2:4.4 % 4.1:1:4.1
            verbose = 0;
            for k=1:length(msigmaa)
                msigma = msigmaa(k);
                mlambda = msigma;
                [x obj iter eps info] = ...
                    SmoothTRWrap(A, b, D, mtheta, mphi, msigma, mkiu, malpha, mlambda);
                if verbose > 0
                    fprintf('objective=%.8f, #iterate=%d, tolerance=%.8f\n', ...
                        obj, iter, eps);
                end
                % for j=1:8
                %    fprintf('x(%d) = %.16f\n', j, x(j));
                % end
                % x = [0.655; 0.208; 1.24e-7; 6.89e-7; 0.232; 7.42e-13; 3.67e-12; 1.71e-8];
                % x = [0.655; 0.208; 0; 0; 0.232; 0; 0; 0];
                % x = [0.545; 0.237; 0; 0.098; 0.165; 0; 0; 0.059];
                RSE = 0;
                [m n] = size(F);
                B = F(:,1:n-1); c = F(:,n);
                
                %
                % tB = B(:,[1 2 5]);
                % z = tB\c;
                % v = sum((tB*z - c).^2)/30;
                
                v = B*x;
                MSE = sum((v-c).^2)/m;
                if MSE < maxMSE
                    maxMSE = MSE;
                    maxSigma = msigma;
                    maxAlpha = malpha;
                    maxX = x;
                end
                
                fprintf('MSE=%f, sigma=%f\n', ...
                    MSE, msigma);
                for j=1:8
                    fprintf('x(%d)=%.6f ', j, x(j));
                end
                fprintf('\n');
            end
        end
        for j=1:8
            fprintf('%.8f\n', maxX(j));
        end
        fprintf('maxMSE=%f, maxSigma=%f, maxAlpha=%f\n', ... 
            maxMSE, maxSigma, maxAlpha);
        
    elseif exps(k) == 4
        
        % 4. Normal distribution, linear regression
        mtheta = 1; mphi = 5; mkiu = 1.0;
        % malpha = 3.7
        % data size and noise scale
        n = 60; sig = 1;
        % repeat times
        ntimes = 100;
        fid = fopen(sprintf('g_solve_info_lin_phi_%d_%d_%.1f_n_%d_sig_%d-sng-5.txt', ... 
            mphi, mtheta, mkiu, n, sig), 'w');
        for malpha = 3.7:0.2:3.7 % 0.2:0.2:2.4 % 3.7:1:3.7 % 2.3:0.2:4.4
            for k = 1:20 % msigma = 20:0.5:23 % k=1:10 % msigma = 10:0.5:12 % 5:1:5% k=1:20 % msigma = 22.0:1:22.0 % 5:1:50
                msigma = 24.5;
                mlambda = msigma;
                sumC  = 0;
                sumIC = 0;
                arrRME = zeros(ntimes, 1);
                tic;
                for k = 1 : ntimes
                    [RME C IC] = run_gaussian_exper(mtheta,...
                        mphi, msigma, mkiu, malpha, mlambda, n, sig, 0, mtheta);
                    arrRME(k) = RME;
                    sumC = sumC + C;
                    sumIC = sumIC + IC;
                end
                sortRME = sort(arrRME);
                MRME = (sortRME(ntimes/2) + sortRME(ntimes/2+1))/2;
                
                fprintf(fid, 'MRME=%f, C=%f, IC=%f, lambda=%f, alpha=%f\n', MRME, ...
                    sumC/ntimes, sumIC/ntimes, mlambda, malpha);
                toc
            end
        end
        fclose(fid);
    elseif exps(k) == 5
        % A = [1.2 0.5 0.7; -0.8 0 -0.5];
        % b = [-1; 1];
        % D = eye(3); 
        % mtheta = 1; msigma = 10; mkiu = 1;
        % [x obj iter eps info] = ...
        %         SmoothTRWrap(A, b, D, mtheta, msigma, mkiu);
        
        % 5. Normal distribution, logistic regression
        mtheta = 2; mphi = 3; mkiu = 0.5;
        
        n = 200; sig = 1;
        
        fid = fopen(sprintf('g_solve_info_lr_phi_%d_%d_%.1f_n_%d_sig_%d-3.txt', ...
            mphi, mtheta, mkiu, n, sig), 'w');
        % for msigma = 2.4:0.5:8.4
        for malpha = 3.7:1:3.7 % 2.3:0.2:4.4 % 0.2:0.2:2.4 % 3.7:1:3.7 % 0.2:0.2:2.4
            for msigma = 5:0.5:50 % 8.4:0.5:28.4
                mlambda = msigma;
                ntimes = 100;
                sumC = 0;
                sumIC = 0;
                arrRME= zeros(ntimes, 1);
                tic
                for k=1:ntimes
                    [RME C IC] = run_gaussian_exper(mtheta,...
                        mphi, msigma, mkiu, malpha, mlambda, n, sig, 0, mtheta);
                    arrRME(k) = RME;
                    sumC = sumC + C;
                    sumIC = sumIC + IC;
                end
                toc
                sortRME = sort(arrRME);
                MRME = (sortRME(ntimes/2) + sortRME(ntimes/2+1))/2;
                
                fprintf(fid, 'MRME=%f, C=%f, IC=%f, lambda=%f, alpha=%f\n', MRME, ...
                    sumC/ntimes, sumIC/ntimes, mlambda, malpha);
            end
        end
        fclose(fid);
    elseif exps(k) == 6
        % counter example
        global g_proc_info;
        global g_nit;
        global g_fid; 
        
        g_fid = fopen('g_proc_info.txt', 'w');
        g_nit = 1;
        % g_proc_info = zeros(1024, 16);
        
        A = [0.25 -0.25 0; 0 0.5 0; 0 0 1]; 
        b = [0; 0.5*(6+1.0/sqrt(5)); 0];
        D = eye(3); 
        
        mtheta = 1; mkiu = 0.5;
        msigma = 1; mphi = 3;
        mlambda = 0; malpha = 0; 
        % H = [-1/8 -1/8; -1/8 5/8 - sqrt(5)/100];
        [x obj iter eps info] = ...
                SmoothTRWrap(A, b, D, mtheta, mphi, msigma, mkiu, malpha, mlambda);
        
        % %
        % adjust the optimal solution
        L = (msigma * mkiu * (1-mkiu) ./(2*sum(A.*A, 1))).^(1/(2-mkiu));
        L = L';
        Lg = ((msigma * mkiu) / (2 * norm(A,2) * norm(b)))^(1/(1-mkiu));
        Lg = Lg * ones(3, 1);
        
        I = find(L < Lg);
        L(I) = Lg(I);
        
        I = find(abs(x) <= L);
        if length(I) > 0
            x(I) = 0;
        end
        % end-adjust
        % % 
        
        % %
        % check the optimality 
        % first order 
        scaled_grad = 2*diag(x)*A'*(A*x-b) + msigma * mkiu * abs(x).^mkiu;
        nscaled_grad = norm(scaled_grad)
        % second order 
        I = find(x > 0);
        H = 2*diag(x)*A'*A*diag(x) + msigma*mkiu*(1-mkiu)*diag(abs(x).^msigma);
        Hh = H(I, I); 
        [R p] = chol(Hh);
        if p == 0
            fprintf('second order condition satisfied\n');
        end
        % end-check
        % % 
        fclose(g_fid);
        
        save('g_proc_info.txt', 'g_proc_info', '-ascii');
        
        x
        obj 
        iter
        
        % plot the figure
        % [X1 X2] = meshgrid(0.5:0.1:3.5, 4:0.1:6);
        % Z = (X1-X2).^2/16 + (X2 - 6 - 1/sqrt(5)).^2/4 + sqrt(abs(X1)) + sqrt(abs(X2));
        % surf(X1, X2, Z);
        % hold on; 
    end
end

% repeat run for checking 
function [RME C IC] = run_gaussian_exper(mtheta, mphi, msigma, mkiu, ... 
    malpha, mlambda, n, sig, verbose, logistic) 

if logistic == 1
    beta = [3; 1.5; 0; 0; 2; 0; 0; 0];
else
    beta = [3; 1.5; 0; 0; 2; 0; 0; 0];
end

dim = length(beta);

A = zeros(n, dim); 
b = zeros(dim, 1);
D = eye(dim);
R = zeros(dim, dim);

if logistic == 1
    for i=1:dim
        for j=1:dim
            R(i,j) = 0.5^(abs(i-j));
        end
    end
else
    for i=1:dim-2
        for j=1:dim-2
            R(i,j) = 0.5^(abs(i-j));
        end
    end
    R(dim-1, dim-1) = 0.25;
    R(dim, dim) = 0.25;
end

% gen rand matrix
mu = zeros(1, dim);
CR = chol(R);
A = repmat(mu, n, 1) + randn(n, dim)*CR;

if logistic == 2
    [nn m] = size(A);
    A(:,m-1:m) = rand(nn, 2) > 0.5;
end

for i=1:n
    if logistic == 2
        u = A(i,:)*beta;
        u = exp(u)/(1+exp(u));
        b(i) = rand() < u;
    else
        b(i) = A(i,:)*beta + sig * randn();
    end
end

[x obj iter eps info] = ...
    SmoothTRWrap(A, b, D, mtheta, mphi, msigma, mkiu, malpha, mlambda);
% x = -x;
dx = (x - beta);

xls = A\b;
dxls = (xls - beta);
mdx = dx'*R*dx;
mdx = mdx*mdx;
mdxls = dxls'*R*dxls;
mdxls = mdxls * mdxls;
RME = mdx/mdxls;

C = 0;
IC = 0;
verbose = 0;
for k=1:dim
    if(beta(k) == 0 && abs(x(k)) <= 1e-8)
        C = C + 1;
    elseif (beta(k) ~= 0 && abs(x(k)) <= 1e-8)
        IC = IC + 1;
    end
    % print the result
    if verbose == 1
        fprintf('estimated=%f, true=%f\n', x(k), beta(k));
    end
end
    
% generate the difference matrix
function D = gen_diff_mat(n) 
D1 = zeros(n-1, n);
for k=1:n-1
    D1(k,k) = 1;
    D1(k,k+1) = -1;
end
D = [kron(D1, eye(n)); kron(eye(n), D1)];

% generate the image
function I = gen_img(n) 
I = zeros(n,n);
% white 
for r=(n/2-5):(n/2-1)
  I = I + gen_circle(n, n/2, n/2, r, 1);
end

% a single cicle
function I = gen_circle(n, x0, y0, r, c)
I = zeros(n,n);
t = 0:0.1:2*pi;
len = length(t);
x = r*sin(t)+x0; 
y = r*cos(t)+y0;
x = round(x); y = round(y); 
for k=1:len
    I(x(k),y(k)) = c;
end