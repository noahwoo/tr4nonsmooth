function Y = get_data(n,R,mu,s)
% GET_DATA - Generate multivariate normal data
%
%   X = get_data(n,R)
%
%   X = get_data(n,R,mu,s)
%
% Generates a sample of n observations of variables having user defined
% correlations R (number of variables is determined by size of R). 
% Optionally, the user may also define the mean (mu) and variance (s) 
% of the variables. If unspecified, these two arguments will default to 
% a mean vector of zeros and variance vector of ones.
%
% For p variables, R is pxp, mu is px1, s is px1.
%
%
% Note that the parameters you supply are the expected values in
% population. Each row returned is a draw from the multivariate normal with mean mu and
% the implied covariance matrix from your supplied correlation matrix.
% Thus the actual sample of data that gets returned will have values
% that are slightly different than the population parameters. A larger sample size 
% will get you closer to the parameters you supply.
%
%
% Example 1: Create two highly correlated variables.
% 
% R = [1 .95; .95 1]
% 
% R =
% 
%     1.0000    0.9500
%     0.9500    1.0000
%     
% X = get_data(100,R);
% 
% 
% corrcoef(X)
% 
% ans =
% 
%     1.0000    0.9587
%     0.9587    1.0000
%
% 


% Example 2: Two highly correlated variables and a third variable negatively correlated to both. Also with user specified means and variances.
% 
% 
% 
% R =
% 
%     1.0000    0.9500   -0.9500
%     0.9500    1.0000   -0.9500
%    -0.9500   -0.9500    1.0000
% 
% 
% mu =
% 
%      5
%      6
%      7
% 
% 
% s =
% 
%     20
%     80
%      3
% 
% 
%  
% >> X = get_data(1000,R,mu,s);
% 
% 
% 
% >> corrcoef(X)
% 
% ans =
% 
%     1.0000    0.9543   -0.9503
%     0.9543    1.0000   -0.9531
%    -0.9503   -0.9531    1.0000
% 
% 
% 
% >> mean(X)
% 
% ans =
% 
%     5.1359    6.2795    6.9641
% 
% 
% 
% >> var(X)
% 
% ans =
% 
%    20.0260   82.9633    3.1008


% It is possible that the correlation matrix you supply can not be used, 
% since it may be impossible for any set of data to have that structure. If this
% happens, it will show up as an error in the Cholesky decomposition of
% your implied covariance matrix.
%


% -------------
% input check
% -------------

if nargin>4
    error('Too many inputs.')
end

if ((nargin ~= 4) && (nargin ~= 2))
    error('Incorrect number of inputs, see help.')
end


% define number of variables p; this can be inferred from either the number
% of rows or columns of matrix R
[p , p2] = size(R) ;

if p ~= p2
    error('Correlation matrix R must be square.')
end

if (nargin == 2) 
    % set defaults for when mu and s not specified
    mu = zeros(p,1) ;
    s = ones(p,1);
end


[n_mu, p_mu]=size(mu);
if ~(n_mu==p && p_mu==1)
    error('Vector mu must be px1.')
end

[n_s, p_s]=size(s);
if ~(n_s==p && p_s==1)
    error('Vector s must be px1.')
end

% -----------------
% end input check
% -----------------



% find the implied covariance matrix S 
DS = sqrt(diag(s));
S = DS*R*DS;

try
    L = chol(S,'lower');
catch
    error('Correlation specified is not possible. Matrix R must be positive definite.'); 
end



% create n random draws (going left to right) of the px1 vector X, where X ~ Np(0,Ip)
% (this means each of the p variables is in a row)
X = randn(p,n);

% now hit each draw (each column vector) with L and add mu 
% (each of the new variables is still in a row)
Y = L*X + repmat(mu,1,n);

% now put each of the p variables into p columns, instead of being in p rows
Y = Y' ;


% end of code

