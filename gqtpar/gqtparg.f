*
*     Gateway function for gqtpar.f from Minpack,
*     by More' and Sorensen.
*     
*     Marcelo Marazzi, Jan 1999.
*
      subroutine mexFunction(nlhs, plhs, nrhs, prhs)
*
*     mexFunction declarations (these don't change).
*
      integer plhs(*), prhs(*)
      integer nlhs, nrhs
*
*     computational function declarations
*
      integer n,maxit,info,iter,NMAX
      parameter (NMAX = 300000)
      double precision wa1(NMAX),wa2(NMAX)
*
*     Miscellaneous parameters declarations
*
      integer nrow,ncol,size,i
      double precision rmaxit, rinfo, riter
*
*     Matlab counterparts of the computational functions 
*     parameters declarations
*
      integer Ap,bp,deltap,rtolp,atolp,maxitp,
     &        infop,xp,redfp,iterp,zp,parp
*
*     Interface subroutines declarations
*
      integer mxGetM,mxGetN,mxIsNumeric,mxCreateFull,
     &        mxGetPr


*     Check for proper number of arguments.

      if (nrhs .ne. 7) then
         call mexErrMsgTxt('Only 7 input arguments allowed.')
      elseif (nlhs .ne. 6) then
         call mexErrMsgTxt('Only 6 output arguments allowed.')
      end if

* Check to ensure input parameters are numeric (not strings).

      if (mxIsNumeric(prhs(1)) .eq. 0) then
         call mexErrMsgTxt('Input A must be a numeric matrix.')
      else if (mxIsNumeric(prhs(2)) .eq. 0) then
         call mexErrMsgTxt('Input b must be a numeric matrix.')
      else if (mxIsNumeric(prhs(3)) .eq. 0) then
         call mexErrMsgTxt('Input delta must be a numeric matrix.')
      else if (mxIsNumeric(prhs(4)) .eq. 0) then
         call mexErrMsgTxt('Input rtol must be a numeric matrix.')
      else if (mxIsNumeric(prhs(5)) .eq. 0) then
         call mexErrMsgTxt('Input atol must be a numeric matrix.')
      else if (mxIsNumeric(prhs(6)) .eq. 0) then
         call mexErrMsgTxt('Input maxit must be a numeric matrix.')
      else if (mxIsNumeric(prhs(7)) .eq. 0) then
         call mexErrMsgTxt('Input par must be a numeric matrix.')
      end if

*     Get the size of the input parameters.

      nrow = mxGetM(prhs(1))
      ncol = mxGetN(prhs(1))
      n = nrow

      if (ncol .ne. nrow) then
         call mexErrMsgTxt('Input A must be square.')
      end if

      nrow = mxGetM(prhs(2))
      ncol = mxGetN(prhs(2))
      
      if ((nrow .ne. 1) .and. (ncol .ne. 1)) then
         call mexErrMsgTxt('Input b must be a vector.')
      end if

      do i = 3,7
         nrow = mxGetM(prhs(i))
         ncol = mxGetN(prhs(i))

         if ((nrow .ne. 1) .or. (ncol .ne. 1)) then
         call mexErrMsgTxt('Last 5 input parameters must be scalars.')
         end if
      end do

      if (n .gt. NMAX) then
         call mexErrMsgTxt('gqtpar: Increase parameter NMAX
     $   (currently too small)')
      end if

*     Get input matrix data

      Ap = mxGetPr(prhs(1))
      bp = mxGetPr(prhs(2))

      deltap = mxGetPr(prhs(3))
      rtolp  = mxGetPr(prhs(4))
      atolp  = mxGetPr(prhs(5))
      maxitp = mxGetPr(prhs(6))
      parp   = mxGetPr(prhs(7))

*     Fortran-Mex files don't handle integer scalars.
*     Thus, we read from Matlab the input parameter "maxit" 
*     and store as a double precision Fortran variable "rmaxit".
*
      size = 1
      call mxCopyPtrToReal8(maxitp, rmaxit, size)

*     Create matrices for the return arguments

      plhs(1) = mxCreateDoubleMatrix(n,1,0)
      plhs(2) = mxCreateDoubleMatrix(1,1,0)
      plhs(3) = mxCreateDoubleMatrix(1,1,0)
      plhs(4) = mxCreateDoubleMatrix(1,1,0)
      plhs(5) = mxCreateDoubleMatrix(n,1,0)
      plhs(6) = mxCreateDoubleMatrix(1,1,0)

      xp     = mxGetPr(plhs(1))
      redfp  = mxGetPr(plhs(2))    
      parp   = mxGetPr(plhs(3))
      iterp  = mxGetPr(plhs(4))    
      zp     = mxGetPr(plhs(5))    
      infop  = mxGetPr(plhs(6))    
*
*     Transform "rmaxit" to the integer Fortran variable "maxit" 
*     so as to pass it to the computational routine (which requires 
*     that this parameter be an integer).
*
      maxit = int(rmaxit)

*     Call the Fortran subroutine

      call gqtpar(n,%val(Ap),n,%val(bp),%val(deltap),
     &     %val(rtolp),%val(atolp),maxit,%val(parp),info,
     &     %val(xp),%val(redfp),iter,%val(zp),wa1,wa2)    
*
*     Transform the Fortran integer output variables "iter"
*     and "info" to the double precision Fortran variables
*     "riter" and "rinfo" respectively. Then transform them
*     into the Matlab variables (pointers) "iterp" and "infop".
*
      riter = dfloat(iter)
      rinfo = dfloat(info)

      size = 1
      call mxCopyReal8ToPtr(riter,iterp,size)
      call mxCopyReal8ToPtr(rinfo,infop,size)

      return
      end
