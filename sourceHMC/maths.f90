module mathsModule
implicit none

! Overloading the functions
interface sigmoid
  module procedure sigmoidVector, sigmoidScalar
end interface

interface invSigmoid
  module procedure invsigmoidVector, invSigmoidScalar
end interface

interface diffSigmoid
  module procedure diffSigmoidVector, diffSigmoidScalar
end interface

contains

!------------------------------------------------------------------
! Calculates log(Gamma(x))
! ALNGAM computes the logarithm of the gamma function.
!
!  Modified:
!
!    13 January 2008
!
!  Author:
!
!    Original FORTRAN77 version by Allan Macleod.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Allan Macleod,
!    Algorithm AS 245,
!    A Robust and Reliable Algorithm for the Logarithm of the Gamma Function,
!    Applied Statistics,
!    Volume 38, Number 2, 1989, pages 397-402.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) XVALUE, the argument of the Gamma function.
!
!    Output, integer ( kind = 4 ) IFAULT, error flag.
!    0, no error occurred.
!    1, XVALUE is less than or equal to 0.
!    2, XVALUE is too big.
!
!    Output, real ( kind = 8 ) ALNGAM, the logarithm of the gamma function of X.
! http://people.sc.fsu.edu/~jburkardt/f_src/f_src.html
!------------------------------------------------------------------
	function alngam ( xvalue, ifault )

	  real ( kind = 8 ) alngam
	  real ( kind = 8 ), parameter :: alr2pi = 0.918938533204673D+00
	  integer ( kind = 4 ) ifault
	  real ( kind = 8 ), dimension ( 9 ) :: r1 = (/ &
	    -2.66685511495D+00, &
	    -24.4387534237D+00, &
	    -21.9698958928D+00, &
	     11.1667541262D+00, &
	     3.13060547623D+00, &
	     0.607771387771D+00, &
	     11.9400905721D+00, &
	     31.4690115749D+00, &
	     15.2346874070D+00 /)
	  real ( kind = 8 ), dimension ( 9 ) :: r2 = (/ &
	    -78.3359299449D+00, &
	    -142.046296688D+00, &
	     137.519416416D+00, &
	     78.6994924154D+00, &
	     4.16438922228D+00, &
	     47.0668766060D+00, &
	     313.399215894D+00, &
	     263.505074721D+00, &
	     43.3400022514D+00 /)
	  real ( kind = 8 ), dimension ( 9 ) :: r3 = (/ &
	    -2.12159572323D+05, &
	     2.30661510616D+05, &
	     2.74647644705D+04, &
	    -4.02621119975D+04, &
	    -2.29660729780D+03, &
	    -1.16328495004D+05, &
	    -1.46025937511D+05, &
	    -2.42357409629D+04, &
	    -5.70691009324D+02 /)
	  real ( kind = 8 ), dimension ( 5 ) :: r4 = (/ &
	     0.279195317918525D+00, &
	     0.4917317610505968D+00, &
	     0.0692910599291889D+00, &
	     3.350343815022304D+00, &
	     6.012459259764103D+00 /)
	  real ( kind = 8 ) x
	  real ( kind = 8 ) x1
	  real ( kind = 8 ) x2
	  real ( kind = 8 ), parameter :: xlge = 5.10D+05
	  real ( kind = 8 ), parameter :: xlgst = 1.0D+30
	  real ( kind = 8 ) xvalue
	  real ( kind = 8 ) y

	  x = xvalue
	  alngam = 0.0D+00
	!
	!  Check the input.
	!
	  if ( xlgst <= x ) then
	    ifault = 2
	    return
	  end if

	  if ( x <= 0.0D+00 ) then
	    ifault = 1
	    return
	  end if

	  ifault = 0
	!
	!  Calculation for 0 < X < 0.5 and 0.5 <= X < 1.5 combined.
	!
	  if ( x < 1.5D+00 ) then

	    if ( x < 0.5D+00 ) then

	      alngam = - log ( x )
	      y = x + 1.0D+00
	!
	!  Test whether X < machine epsilon.
	!
	      if ( y == 1.0D+00 ) then
	        return
	      end if

	    else

	      alngam = 0.0D+00
	      y = x
	      x = ( x - 0.5D+00 ) - 0.5D+00

	    end if

	    alngam = alngam + x * (((( &
	        r1(5)   * y &
	      + r1(4) ) * y &
	      + r1(3) ) * y &
	      + r1(2) ) * y &
	      + r1(1) ) / (((( &
	                  y &
	      + r1(9) ) * y &
	      + r1(8) ) * y &
	      + r1(7) ) * y &
	      + r1(6) )

	    return

	  end if
	!
	!  Calculation for 1.5 <= X < 4.0.
	!
	  if ( x < 4.0D+00 ) then

	    y = ( x - 1.0D+00 ) - 1.0D+00

	    alngam = y * (((( &
	        r2(5)   * x &
	      + r2(4) ) * x &
	      + r2(3) ) * x &
	      + r2(2) ) * x &
	      + r2(1) ) / (((( &
	                  x &
	      + r2(9) ) * x &
	      + r2(8) ) * x &
	      + r2(7) ) * x &
	      + r2(6) )
	!
	!  Calculation for 4.0 <= X < 12.0.
	!
	  else if ( x < 12.0D+00 ) then

	    alngam = (((( &
	        r3(5)   * x &
	      + r3(4) ) * x &
	      + r3(3) ) * x &
	      + r3(2) ) * x &
	      + r3(1) ) / (((( &
	                  x &
	      + r3(9) ) * x &
	      + r3(8) ) * x &
	      + r3(7) ) * x &
	      + r3(6) )
	!
	!  Calculation for 12.0 <= X.
	!
	  else

	    y = log ( x )
	    alngam = x * ( y - 1.0D+00 ) - 0.5D+00 * y + alr2pi

	    if ( x <= xlge ) then

	      x1 = 1.0D+00 / x
	      x2 = x1 * x1

	      alngam = alngam + x1 * ( ( &
	             r4(3)   * &
	        x2 + r4(2) ) * &
	        x2 + r4(1) ) / ( ( &
	        x2 + r4(5) ) * &
	        x2 + r4(4) )

	    end if

	  end if

	  return
	end function alngam

!*****************************************************************************80
!
!! DIGAMA calculates DIGAMMA ( X ) = d ( LOG ( GAMMA ( X ) ) ) / dX
!
!  Modified:
!
!    03 June 2013
!
!  Author:
!
!    Original FORTRAN77 version by Jose Bernardo.
!    FORTRAN90 version by John Burkardt.
!
!  Reference:
!
!    Jose Bernardo,
!    Algorithm AS 103:
!    Psi ( Digamma ) Function,
!    Applied Statistics,
!    Volume 25, Number 3, 1976, pages 315-317.
!
!  Parameters:
!
!    Input, real ( kind = 8 ) X, the argument of the digamma function.
!    0 < X.
!
!    Output, integer IFAULT, error flag.
!    0, no error.
!    1, X <= 0.
!
!    Output, real ( kind = 8 ) DIGAMA, the value of the digamma function at X.
!
! http://people.sc.fsu.edu/~jburkardt/f_src/f_src.html
	function digama ( x, ifault )


	implicit none

	real ( kind = 8 ), parameter :: euler_mascheroni = 0.57721566490153286060D+00
	real ( kind = 8 ) digama
	integer ( kind = 4 ) ifault
	real ( kind = 8 ) r
	real ( kind = 8 ) x
	real ( kind = 8 ) x2
	!
	!  Check the input.
	!
		if ( x <= 0.0D+00 ) then
			digama = 0.0D+00
			ifault = 1
			return
		end if
		!
		!  Initialize.
		!
		ifault = 0
		x2 = x
		digama = 0.0D+00
		!
		!  Approximation for small argument.
		!
		if ( x2 <= 0.00001D+00 ) then
			digama = - euler_mascheroni - 1.0D+00 / x2
			return
		end if
		!
		!  Reduce to DIGAMA(X + N).
		!
		do while ( x2 < 8.5D+00 )
			digama = digama - 1.0D+00 / x2
			x2 = x2 + 1.0D+00
		end do
		!
		!  Use Stirling's (actually de Moivre's) expansion.
		!
		r = 1.0D+00 / x2
		digama = digama + log ( x2 ) - 0.5D+00 * r
		r = r * r
		digama = digama &
			- r * ( 1.0D+00 / 12.0D+00 &
			- r * ( 1.0D+00 / 120.0D+00 &
			- r *   1.0D+00 / 252.0D+00 ) )

		return
	end function digama

!-------------------------------------------------
! This function computes log(1+x)
!-------------------------------------------------
	function log1p(x)
	real(kind=8) :: x, log1p 

		if (1.d0 + x .eq. 1.d0) then
		   log1p = x
		else
			log1p = log(1.d0 + x) * x / ((1.d0 + x) - 1.d0)
		endif		
	end function log1p

!-------------------------------------------------
! This function computes log(sum(x)) given a vector of the logarithm of x
!-------------------------------------------------
	function logSum(logx, logy)
	real(kind=8) :: logx, logy, logsum, xlog, ylog, negDiff

		xlog = logx
		ylog = logy
		
! Make X the maximum
		if (ylog > xlog) then
			xlog = logy
			ylog = logx
		endif

! X is now bigger
		if (xlog == -tiny(1.d0)) then
			logSum = xlog
			return
		endif

! Test if the number we are adding is many orders of magnitude below
		negDiff = ylog - xlog
		if (negDiff < -30.d0) then
			logSum = xlog
			return
		else
			logSum = xlog + log1p(exp(negDiff))
			return
		endif
			
	end function logSum
	
!-------------------------------------------------
! This function computes log(sum(x)) given a vector of the logarithm of x
!-------------------------------------------------
	function logSumVector(logx)
	real(kind=8) :: logx(:), logSumVector, t
	integer :: i, n
		n = size(logx)
		
		logSumVector = logx(1)
		do i = 2, n			
			logSumVector = logSum(logSumVector, logx(i))
		enddo
		
	end function logSumVector
	
!-----------------------------------------------------------------------
! Return a sigmoid function
!-----------------------------------------------------------------------
	function sigmoidVector(x, lower, upper)
	real(kind=8) :: x(:), sigmoidVector(size(x)), lower, upper
		sigmoidVector = lower + (upper-lower) / (1.d0 + exp(-x))
	end function sigmoidVector
	
!-----------------------------------------------------------------------
! Return a sigmoid function
!-----------------------------------------------------------------------
	function sigmoidScalar(x, lower, upper)
	real(kind=8) :: x, sigmoidScalar, lower, upper
		sigmoidScalar = lower + (upper-lower) / (1.d0 + exp(-x))
	end function sigmoidScalar

!-----------------------------------------------------------------------
! Return the inverse of the sigmoid function
!-----------------------------------------------------------------------
	function invSigmoidVector(x, lower, upper)
	real(kind=8) :: x(:), invsigmoidVector(size(x)), lower, upper
		invSigmoidVector = log( (lower - x) / (x - upper) )
	end function invSigmoidVector
	
!-----------------------------------------------------------------------
! Return the inverse of the sigmoid function
!-----------------------------------------------------------------------
	function invSigmoidScalar(x, lower, upper)
	real(kind=8) :: x, invsigmoidScalar, lower, upper
		invSigmoidScalar = log( (lower - x) / (x - upper) )
	end function invSigmoidScalar

!-----------------------------------------------------------------------
! Return the derivative of the sigmoid function
!-----------------------------------------------------------------------
	function diffSigmoidVector(x, lower, upper)
	real(kind=8) :: x(:), diffSigmoidVector(size(x)), lower, upper
		diffSigmoidVector = (upper-lower) * exp(-x) / (1.d0+exp(-x))**2
	end function diffSigmoidVector
	
!-----------------------------------------------------------------------
! Return the derivative of the sigmoid function
!-----------------------------------------------------------------------
	function diffSigmoidScalar(x, lower, upper)
	real(kind=8) :: x, diffSigmoidScalar, lower, upper
		diffSigmoidScalar = (upper-lower) * exp(-x) / (1.d0+exp(-x))**2
	end function diffSigmoidScalar

end module mathsModule