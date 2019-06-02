      subroutine gqtpar(n,a,lda,b,delta,rtol,atol,maxit,
     *                  par,info,x,redf,iter,z,wa1,wa2)
      integer n,lda,maxit,info,iter
      double precision delta,rtol,atol,par,redf
      double precision a(lda,n),b(n),x(n),z(n),wa1(n),wa2(n)
c     **********
c
c     Subroutine gqtpar
c
c     Given an n by n symmetric matrix A, an n-vector b, and a
c     positive number delta, this subroutine determines a vector
c     x which minimizes the quadratic function
c
c           f(x) = (x,Ax)/2 + (b,x)
c
c     subject to the constraint
c
c           norm(x) .le. delta.
c
c     The norm is the Euclidean norm and (u,v) is the inner
c     product between vectors u and v.
c
c     The subroutine statement is
c
c       subroutine gqtpar(n,a,lda,b,delta,rtol,atol,maxit,
c                         par,info,x,redf,iter,z,wa1,wa2)
c
c     where
c
c       n is a positive integer input variable set to the order of A.
c
c       a is an n by n array. On input the full upper triangle must
c         contain the full upper triangle of the symmetric matrix A.
c         On output the array contains the matrix A.
c
c       lda is a positive integer input variable not less than n
c         which specifies the leading dimension of the array a.
c
c       b is an input array of length n which must contain the
c         elements of the vector b.
c
c       delta is a positive input variable which specifies an upper
c         bound on the euclidean norm of x.
c
c       rtol is a positive input variable which specifies the
c         relative accuracy desired in the solution. If xsol is the
c         solution and redf(x) is the reduction in the objective
c         function at x, then termination occurs when
c
c               norm(x) .le. (1 + rtol)*delta
c
c               redf(x)  .ge. ((1 - rtol)**2)*redf(xsol)
c
c       atol is a nonnegative input variable which specifies the
c         absolute accuracy desired in the solution. If xsol is the
c         solution and redf(x) is the reduction in the objective
c         function at x, then termination occurs when
c
c               norm(x) .le. (1 - rtol)*delta
c
c               max(redf(x),redf(xsol)) .le. atol
c
c       maxit is a positive integer input variable which specifies
c         the maximum number of iterations.
c
c       par is a nonnegative variable. On input par contains an
c         initial estimate of the Lagrange multiplier corresponding
c         to the constraint norm(x) .le. delta.  On output par
c         contains the final estimate of this multiplier.
c
c       info is an integer output variable set as follows.
c
c         info = 1   The reduction redf has the relative accuracy
c                    specified by rtol.
c
c         info = 2   The reduction redf has the absolute accuracy
c                    specified by atol.
c
c         info = 3   Rounding errors prevent further progress.
c
c         info = 4   Number of iterations has reached maxit.
c
c       x is an output array of length n which contains the final
c         estimate of the solution vector.
c
c       redf is an output variable which contains the reduction in
c         the objective function at the output x.
c
c       iter is a nonnegative integer output variable set to the
c         number of iterations.
c
c       z is an output array of length n which contains the final
c          estimate of a direction of negative curvature.
c
c       wa1 and wa2 are work arrays of length n.
c
c     Subprograms called
c
c       Minpack-supplied ... chfac,enorm,estsv,rsol
c
c       Fortran-supplied ... dble
c                            abs,max,min,sign,sqrt
c
c     Argonne National Laboratory. Minpack Project. February 1984.
c
c     **********
      integer i,j,nrank
      logical rednc
      double precision alpha,anorm,bnorm,parc,parl,pars,paru,parf,
     *       prod,rxnorm,rznorm,temp,xnorm
      double precision enorm
c
c     Initialization.
c
      parf = 0.0
      xnorm = 0.0
      rxnorm = 0.0
      rednc = .false.
      do 20 j = 1, n
         x(j) = 0.0
         z(j) = 0.0
         wa1(j) = a(j,j)
         do 10 i = j+1, n
            a(i,j) = a(j,i)
   10       continue
   20    continue
c
c     Calculate the l1-norm of A and the l2-norm of b.
c
      anorm = 0.0
      do 40 j = 1, n
         temp = 0.0
         do 30 i = 1, n
            temp = temp + abs(a(i,j))
   30       continue
         wa2(j) = temp - abs(wa1(j))
         anorm = max(anorm,temp)
   40    continue
      bnorm = enorm(n,b)
c
c     Calculate a lower bound, pars, for the domain of the problem.
c     Also calculate an upper bound, paru, and a lower bound, parl,
c     for the Lagrange multiplier.
c
      pars = -anorm
      parl = -anorm
      paru = -anorm
      do 50 j = 1, n
         pars = max(pars,-wa1(j))
         parl = max(parl,wa1(j)+wa2(j))
         paru = max(paru,-wa1(j)+wa2(j))
   50    continue
      parl = max(dble(0.0),bnorm/delta-parl,pars)
      paru = max(dble(0.0),bnorm/delta+paru)
c
c     If the input par lies outside of the interval (parl,paru),
c     set par to the closer endpoint.
c
      par = max(par,parl)
      par = min(par,paru)
c
c     Special case: parl = paru.
c
      paru = max(paru,(1+rtol)*parl)
c
c     Beginning of an iteration.
c
      info = 0
      do 140 iter = 1, maxit
c
c        Safeguard par.
c
         if (par .le. pars .and. paru .gt. 0.0)
     *      par = max(dble(0.001),sqrt(parl/paru))*paru
c
c        Attempt to compute the Cholesky factorization.
c
         do 70 j = 1, n
            do 60 i = 1, j-1
               a(i,j) = a(j,i)
   60          continue
            a(j,j) = wa1(j) + par
   70       continue
         call chfac(n,a,lda,wa2,nrank)
c
c        Compute an approximate solution x.
c
         if (nrank .eq. n) then
            parf = par
            call rsol(n,a,lda,b,wa2,2)
            rxnorm = enorm(n,wa2)
            call rsol(n,a,lda,wa2,x,1)
            do 80 i = 1, n
               x(i) = -x(i)
   80          continue
            xnorm = enorm(n,x)
c
c           Convergence test.
c
            if (abs(xnorm-delta) .le. rtol*delta .or.
     *         par .eq. 0.0 .and. xnorm .le. (1+rtol)*delta) info = 1
c
c           Compute a direction of negative curvature and
c           use this information to improve pars.
c
            call estsv(n,a,lda,rznorm,z)
            pars = max(pars,par-rznorm**2)
c
c           Compute a negative curvature solution of the form
c           x + alpha*z where norm(x + alpha*z) = delta.
c
            rednc = .false.
            if (xnorm .le. delta) then
               prod = 0.0
               do 90 i = 1, n
                  prod = prod + z(i)*(x(i)/delta)
   90             continue
               temp = (delta - xnorm)*((delta + xnorm)/delta)
               alpha = temp/(abs(prod) + sqrt(prod**2 + temp/delta))
               rznorm = alpha*rznorm
               alpha = sign(alpha,prod)
               if ((rznorm/delta)**2 + par*(xnorm/delta)**2 .le. par)
     *            rednc = .true.
c
c              Convergence tests.
c
               if ((rznorm/delta)**2 .le.
     *            rtol*(2-rtol)*(par + (rxnorm/delta)**2)) then
                  info = 1
               else if ((par + (rxnorm/delta)**2)/2 .le. 
     *            (atol/delta)/delta .and. iter .ne. 1) then
                  if (info .ne. 1) info = 2
                  end if
               end if
c
c           Compute the Newton correction parc to par.
c
            if (xnorm .eq. 0.0) then
               parc = -par
            else
               do 100 j = 1, n
                  wa2(j) = x(j)/xnorm
  100             continue
               call rsol(n,a,lda,wa2,wa2,2)
               temp = enorm(n,wa2)
               parc = (((xnorm - delta)/delta)/temp)/temp
               end if
c
c           Update parl or paru.
c
            if (xnorm .gt. delta) parl = max(parl,par)
            if (xnorm .lt. delta) paru = min(paru,par)
c
c        Use the rank information from the Cholesky decomposition.
c
         else
            temp = enorm(nrank+1,wa2)
            parc = -(a(nrank+1,nrank+1)/temp)/temp
            pars = max(pars,par+parc)

c           Modify paru to guarantee that nrank = n.

            paru = max(paru,(1+2*rtol)*pars)
            end if
c
c        Use pars to update parl.
c
         parl = max(parl,pars)
c
c        Termination tests.
c
         if (info .eq. 0) then
            if (iter .eq. maxit) info = 4
            if (paru - pars .le. rtol*(2-rtol)*pars) info = 3
            if (paru .eq. 0.0) info = 2
            end if
c
c        If exiting, store the best approximation and
c        restore the upper triangle of A.
c
         if (info .ne. 0) then
            par = parf
            redf = (rxnorm**2 + (par*xnorm)*xnorm)/2
            if (rednc) then
               redf = ((rxnorm**2 + (par*delta)*delta) - rznorm**2)/2
               do 110 i = 1, n
                  x(i) = x(i) + alpha*z(i)
  110             continue
               end if
            do 130 j = 1, n
               do 120 i = 1, j-1
                  a(i,j) = a(j,i)
  120             continue
               a(j,j) = wa1(j)
  130          continue
            return
            end if
c
c        Compute an improved estimate for par.
c
         par = max(parl,par+parc)
c
c        End of an iteration.
c
  140    continue
c
c     Last card of subroutine gqtpar.
c
      end
      subroutine chfac(n,a,lda,z,info)
      integer n,lda,info
      double precision a(lda,n),z(n)
c     **********
c
c     Subroutine chfac
c
c     Given an n by n symmetric matrix A, this subroutine computes
c     a modified Cholesky decomposition for A.
c
c     If the Cholesky decompositon of A exists then this subroutine
c     computes an upper triangular matrix R (the Cholesky factor)
c     such that
c
c                T
c           A = R R
c
c     If the Cholesky decomposition of A does not exist, then there
c     is a largest leading submatrix of A which is not positive
c     definite. This subroutine computes a Cholesky decompositon of
c     this submatrix, and a direction z of negative curvature for A.
c     Thus, (z,Az) is not positive, where (x,y) is the inner product
c     between vectors x and y.
c
c     The subroutine statement is
c
c       subroutine chfac(n,a,lda,z,info)
c
c     where
c
c       n is a positive integer input variable set to the order of A.
c
c       a is an n by n array. On input the full upper triangle must
c         contain the full upper triangle of the matrix A. If the
c         Cholesky decomposition exists, then on output the full
c         upper triangle contains the full upper triangle of the
c         matrix R. Otherwise, on output the first info columns in
c         the full upper triangle contains the full upper triangle
c         of the Cholesky factor for the order info leading
c         submatrix of A.
c
c       lda is a positive integer input variable not less than n
c         which specifies the leading dimension of the array a.
c
c       z is an output array of length n. If A does not have a
c         Cholesky decomposition then z is a vector such that
c         (z,Az) is not positive. Otherwise, z is not referenced.
c
c       info is an integer output variable which specifies the
c         order of the largest leading submatrix of A which is
c         positive definite. If no leading submatrix is positive
c         definite, then info = 0.
c
c     Subprograms called
c
c       Fortran-supplied ... sqrt
c
c       Minpack-supplied ... rsol
c
c     Argonne National Laboratory. Minpack Project. January 1984.
c
c     **********
      integer i,j,k
      double precision sum1,sum2
c
c     Factorization of A.
c
      do 40 j = 1, n
         sum2 = 0.0
         do 20 i = 1, j-1
            sum1 = 0.0
            do 10 k = 1, i-1
               sum1 = sum1 + a(k,j)*a(k,i)
   10          continue
            a(i,j) = (a(i,j) - sum1)/a(i,i)
            sum2 = sum2 + a(i,j)**2
   20       continue
c
c        If it is not possible to complete the factorization then
c        compute a direction of negative curvature and terminate.
c
         if (a(j,j) .gt. sum2) then
            a(j,j) = sqrt(a(j,j)-sum2)
         else
            a(j,j) = a(j,j) - sum2
            info = j - 1
            call rsol(info,a,lda,a(1,j),z,1)
            z(j) = -1.0
            do 30 i = j+1, n
               z(i) = 0.0
   30          continue
            return
            end if
   40    continue
      info = n
      return
c
c     Last card of subroutine chfac.
c
      end
      subroutine rsol(n,r,ldr,b,x,job)
      integer n,ldr,job
      double precision r(ldr,n),b(n),x(n)
c     **********
c
c     Subroutine rsol
c
c     Given an n by n upper triangular matrix R and an n-vector b,
c     this subroutine solves the system Rx = b if job = 1 and the
c     system (R transpose)x = b if job = 2. It is assumed that R
c     is nonsingular.
c
c     The subroutine statement is
c
c       subroutine rsol(n,r,ldr,b,x,job)
c
c     where
c
c       n is a positive integer input variable set to the order of R.
c
c       r is an n by n input array. The full upper triangle must
c         contain the full upper triangle of the matrix R.
c         The strict lower triangle is not referenced.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       b is an input array of length n which must contain the
c         elements of the vector b.
c
c       x is an output array of length n. If job = 1 then x contains
c         the solution to the system Rx = b, while if job = 2 the x
c         contains the solution to the system (R transpose)x = b.
c
c       job is a positive input variable. The system Rx = b is solved
c         if job = 1, and (R transpose)x = b is solved if job = 2.
c
c     Argonne National Laboratory. Minpack Project. January 1984.
c
c     **********
      integer i,j
      double precision sum,temp
c
c     Solve Rx = b.
c
      if (job .eq. 1) then
         do 10 i = 1, n
            x(i) = b(i)
   10       continue
         do 30 j = n, 1, -1
            x(j) = x(j)/r(j,j)
            temp = x(j)
            do 20 i = 1, j-1
               x(i) = x(i) - r(i,j)*temp
   20          continue
   30       continue
         return
c
c     Solve (R transpose)x = b.
c
      else if (job .eq. 2) then
         do 50 j = 1, n
            sum = b(j)
            do 40 i = 1, j-1
               sum = sum - r(i,j)*x(i)
   40          continue
            x(j) = sum/r(j,j)
   50       continue
         return
         end if
c
c     Last card of subroutine rsol.
c
      end
      subroutine estsv(n,r,ldr,svmin,z)
      integer n,ldr
      double precision svmin
      double precision r(ldr,n),z(n)
c     **********
c
c     Subroutine estsv
c
c     Given an n by n upper triangular matrix R, this subroutine
c     estimates the smallest singular value and the associated
c     singular vector of R.
c
c     In the algorithm a vector e is selected so that the solution y
c     to the system (R transpose)y = e is large. The choice of sign
c     for the components of e cause maximal local growth in the
c     components of y as the forward substitution proceeds.
c     The vector z is the solution of the system Rz = y, and the
c     estimate svmin is enorm(y)/enorm(z).
c
c     The subroutine statement is
c
c       subroutine estsv(n,r,ldr,svmin,z)
c
c     where
c
c       n is a positive integer input variable set to the order of R.
c
c       r is an n by n input array. The full upper triangle must
c         contain the full upper triangle of the matrix R.
c         The strict lower triangle is not referenced.
c
c       ldr is a positive integer input variable not less than n
c         which specifies the leading dimension of the array r.
c
c       svmin is a nonnegative output variable which contains an
c         estimate for the smallest singular value of R.
c
c       z is an output array of length n which defines a singular
c         vector associated with the estimate svmin. On output
c         enorm(Rz) = svmin and enorm(z) = 1.
c
c     Subprograms called
c
c       Minpack-supplied ... enorm
c
c       Fortran-supplied ... abs,sign
c
c     Argonne National Laboratory. Minpack Project. January 1984.
c
c     **********
      integer i,j
      double precision e,s,sm,temp,w,wm,ynorm,znorm
      double precision enorm
      do 10 i = 1, n
         z(i) = 0.0
   10    continue
c
c     This choice of e makes the algorithm scale invariant.
c
      e = abs(r(1,1))
      if (e .eq. 0.0) then
         svmin = 0.0
         z(1) = 1.0
         return
         end if
c
c     Solve (R transpose)y = e.
c
      do 50 i = 1, n
         e = sign(e,-z(i))
c
c        Scale y. The factor of 0.01 reduces the number of scalings.
c
         if (abs(e-z(i)) .gt. abs(r(i,i))) then
            temp = min(dble(0.01),abs(r(i,i))/abs(e-z(i)))
            do 20 j = 1, n
               z(j) = temp*z(j)
   20          continue
            e = temp*e
            end if
c
c        Determine the two possible choices of y(i).
c
         if (r(i,i) .eq. 0.0) then
            w = 1.0
            wm = 1.0
         else
            w =  (e - z(i))/r(i,i)
            wm =  -(e + z(i))/r(i,i)
            end if
c
c       Choose y(i) based on the predicted value
c       of y(j) for j = i+1,...,n.
c
         s = abs(e-z(i))
         sm = abs(e+z(i))
         do 30 j = i+1, n
            sm = sm + abs(z(j)+wm*r(i,j))
            z(j) = z(j) + w*r(i,j)
            s = s + abs(z(j))
   30       continue
         if (s .lt. sm) then
            temp = wm - w
            w = wm
            do 40 j = i+1, n
               z(j) = z(j) + temp*r(i,j)
   40          continue
            end if
         z(i) = w
   50    continue
      ynorm = enorm(n,z)
c
c     Solve Rz = y.
c
      do 80 j = n, 1, -1
c
c        Scale z.
c
         if (abs(z(j)) .gt. abs(r(j,j))) then
            temp = min(dble(0.01),abs(r(j,j))/abs(z(j)))
            do 60 i = 1, n
               z(i) = temp*z(i)
   60          continue
            ynorm = temp*ynorm
            end if
         if (r(j,j) .eq. 0.0) then
             z(j) = 1.0
         else
             z(j) = z(j)/r(j,j)
             end if
         temp = z(j)
         do 70 i = 1, j-1
            z(i) = z(i) - temp*r(i,j)
   70       continue
   80    continue
c
c     Compute svmin and normalize z.
c
      znorm = enorm(n,z)
      svmin = ynorm/znorm
      do 90 i = 1, n
         z(i) = z(i)/znorm
   90    continue
      return
c
c     Last card of subroutine estsv.
c
      end
      double precision function enorm(n,x)
      integer n
      double precision x(n)
c     **********
c
c     function enorm
c
c     given an n-vector x, this function calculates the
c     euclidean norm of x.
c
c     the euclidean norm is computed by accumulating the sum of
c     squares in three different sums. the sums of squares for the
c     small and large components are scaled so that no overflows
c     occur. non-destructive underflows are permitted. underflows
c     and overflows do not occur in the computation of the unscaled
c     sum of squares for the intermediate components.
c     the definitions of small, intermediate and large components
c     depend on two constants, rdwarf and rgiant. the main
c     restrictions on these constants are that rdwarf**2 not
c     underflow and rgiant**2 not overflow. the constants
c     given here are suitable for every known computer.
c
c     the function statement is
c
c       double precision function enorm(n,x)
c
c     where
c
c       n is a positive integer input variable.
c
c       x is an input array of length n.
c
c     subprograms called
c
c       fortran-supplied ... dabs,dsqrt
c
c     argonne national laboratory. minpack project. march 1980.
c     burton s. garbow, kenneth e. hillstrom, jorge j. more
c
c     **********
      integer i
      double precision agiant,floatn,one,rdwarf,rgiant,s1,s2,s3,xabs,
     *                 x1max,x3max,zero
      data one,zero,rdwarf,rgiant /1.0d0,0.0d0,5.422d-20,1.304d19/
      s1 = zero
      s2 = zero
      s3 = zero
      x1max = zero
      x3max = zero
      floatn = n
      agiant = rgiant/floatn
      do 90 i = 1, n
         xabs = dabs(x(i))
         if (xabs .gt. rdwarf .and. xabs .lt. agiant) go to 70
            if (xabs .le. rdwarf) go to 30
c
c              sum for large components.
c
               if (xabs .le. x1max) go to 10
                  s1 = one + s1*(x1max/xabs)**2
                  x1max = xabs
                  go to 20
   10          continue
                  s1 = s1 + (xabs/x1max)**2
   20          continue
               go to 60
   30       continue
c
c              sum for small components.
c
               if (xabs .le. x3max) go to 40
                  s3 = one + s3*(x3max/xabs)**2
                  x3max = xabs
                  go to 50
   40          continue
                  if (xabs .ne. zero) s3 = s3 + (xabs/x3max)**2
   50          continue
   60       continue
            go to 80
   70    continue
c
c           sum for intermediate components.
c
            s2 = s2 + xabs**2
   80    continue
   90    continue
c
c     calculation of norm.
c
      if (s1 .eq. zero) go to 100
         enorm = x1max*dsqrt(s1+(s2/x1max)/x1max)
         go to 130
  100 continue
         if (s2 .eq. zero) go to 110
            if (s2 .ge. x3max)
     *         enorm = dsqrt(s2*(one+(x3max/s2)*(x3max*s3)))
            if (s2 .lt. x3max)
     *         enorm = dsqrt(x3max*((s2/x3max)+(x3max*s3)))
            go to 120
  110    continue
            enorm = x3max*dsqrt(s3)
  120    continue
  130 continue
      return
c
c     last card of function enorm.
c
      end
