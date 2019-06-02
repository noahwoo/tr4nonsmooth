function [x,redf,par,iter,z,info] = gqtpar(A,b,delta,rtol,atol,maxit,par)
% [x,redf,par,iter,z,info] = gqtpar(A,b,delta,rtol,atol,maxit,par)
% by More' and Sorensen from Minpack.
%
%     Given an n by n symmetric matrix A, an n-vector b, and a
%     positive number delta, this subroutine determines a vector
%     x which minimizes the quadratic function
%
%           f(x) = (x,Ax)/2 + (b,x)
%
%     subject to the constraint
%
%           norm(x) .le. delta.
%
%     The norm is the Euclidean norm and (u,v) is the inner
%     product between vectors u and v.
%
%       A is an n by n array. On input the full upper triangle must
%         contain the full upper triangle of the symmetric matrix A.
%         On output the array contains the matrix A.
%
%       b is an input array of length n which must contain the
%         elements of the vector b.
%       delta is a positive input variable which specifies an upper
%         bound on the euclidean norm of x.
%
%       rtol is a positive input variable which specifies the
%         relative accuracy desired in the solution. If xsol is the
%         solution and redf(x) is the reduction in the objective
%         function at x, then termination occurs when
%
%               norm(x) .le. (1 + rtol)*delta
%
%               redf(x)  .ge. ((1 - rtol)**2)*redf(xsol)
%
%       atol is a nonnegative input variable which specifies the
%         absolute accuracy desired in the solution. If xsol is the
%         solution and redf(x) is the reduction in the objective
%         function at x, then termination occurs when
%
%               norm(x) .le. (1 - rtol)*delta
%
%               max(redf(x),redf(xsol)) .le. atol
%
%       maxit is a positive integer input variable which specifies
%         the maximum number of iterations.
%
%       par is a nonnegative variable. On input par contains an
%         initial estimate of the Lagrange multiplier corresponding
%         to the constraint norm(x) .le. delta.  On output par
%         contains the final estimate of this multiplier.
%
%       info is an integer output variable set as follows.
%
%         info = 1   The reduction redf has the relative accuracy
%                    specified by rtol.
%
%         info = 2   The reduction redf has the absolute accuracy
%                    specified by atol.
%
%         info = 3   Rounding errors prevent further progress.
%
%         info = 4   Number of iterations has reached maxit.
%
%       x is an output array of length n which contains the final
%         estimate of the solution vector.
%
%       redf is an output variable which contains the reduction in
%         the objective function at the output x.
%
%       iter is a nonnegative integer output variable set to the
%         number of iterations.
%
%       z is an output array of length n which contains the final
%          estimate of a direction of negative curvature.
%

% This m-file only contains the help message. The source code is contained
% in the original Fortran subroutine gqtpar.f and in the gateway subroutine
% gqtparg.f.
