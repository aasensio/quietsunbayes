module mcmc_class_hierarchical
use mcmc_class
implicit none

	type pixelType		
		real(kind=8) :: noise
	end type pixelType
	
	type, EXTENDS (mcmc) :: mcmc_hierarchical

		integer :: npixels
		type(pixelType), pointer :: pixel(:)
		real(kind=8) :: BMin, BMax, muMin, muMax, fMin, fMax, phiMin, phiMax
		
		real(kind=8), pointer, dimension(:) :: CV2, CV3, CQ2, CQ3, CU2, CU3
	
 		contains
 			procedure :: evalPosterior => logPosteriorWeakField
 			procedure :: initialValues => initialValuesWeakField
	end type mcmc_hierarchical
	
	real(kind=8), allocatable :: B_i(:), mu_i(:), phi_i(:), f_i(:), hyperB_i(:), hypermu_i(:), hyperphi_i(:), hyperf_i(:)
	real(kind=8), allocatable :: c2p(:), s2p(:), sqrtMu(:)
		
	contains
	
!------------------------------------------------------------------
! Calculates log(Gamma(x))
!! ALNGAM computes the logarithm of the gamma function.
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

!------------------------------------------------------------------
! Set initial values for the parameters
!------------------------------------------------------------------
	subroutine initialValuesWeakField(this)
	class(mcmc_hierarchical), intent(inout) :: this
	integer :: loop, i, j
	real(kind=8) :: value

		loop = 1
		
! B
		do i = 1, this%npixels
			this%pars(loop) = 200.d0
			loop = loop + 1
		enddo

! mu
		do i = 1, this%npixels
			this%pars(loop) = 0.2
			loop = loop + 1
		enddo

! phi
		do i = 1, this%npixels
			this%pars(loop) = 0.2d0
			loop = loop + 1
		enddo

! f
		do i = 1, this%npixels
			this%pars(loop) = 0.1d0
			loop = loop + 1
		enddo

! hyperB
		do i = 1, 2
			this%pars(loop) = 1.d0
			loop = loop + 1
		enddo
		
! hypermu
		do i = 1, 2
			this%pars(loop) = 1.d0
			loop = loop + 1
		enddo
		
! hyperf
		do i = 1, 2
			this%pars(loop) = 1.d0
			loop = loop + 1
		enddo
		
! hyperphi
		do i = 1, 2
			this%pars(loop) = 1.d0
			loop = loop + 1
		enddo

		this%ls = 1.d0
		
	end subroutine initialValuesWeakField

!------------------------------------------------------------------
! Evaluate the posterior. This function overrides the function in the parent
!------------------------------------------------------------------
	function logPosteriorWeakField(this, outBounds) result (logP)
	class(mcmc_hierarchical), intent(in) :: this
	integer, intent(inout) :: outBounds
	real(kind=8) :: logP, nu
	integer :: i, j, loop, ierr

! The vector of trial parameters enters with the following order
! B(1..npixels)
! mu(1..npixels)
! f(1..npixels)
! phi(1..npixels)

! hyperB(1..2)
! hypermu(1..2)
! hyperf(1..2)
! hyperphi(1..2)


! Fill the vectors to simplify the notation and test for the boundaries
		outBounds = 0
		loop = 1
		
		do i = 1, this%npixels
			B_i(i) = this%trial(loop)
			loop = loop + 1
			if (B_i(i) < this%BMin .or. B_i(i) > this%BMax) then
				outBounds = 1
				return
			endif
		enddo
		
		do i = 1, this%npixels
			mu_i(i) = this%trial(loop)
			loop = loop + 1
			if (mu_i(i) < this%muMin .or. mu_i(i) > this%muMax) then
				outBounds = 1
				return
			endif
		enddo
		
		do i = 1, this%npixels
			f_i(i) = this%trial(loop)
			loop = loop + 1
			if (f_i(i) < this%fMin .or. f_i(i) > this%fMax) then
				outBounds = 1
				return
			endif
		enddo
		
		do i = 1, this%npixels
			phi_i(i) = this%trial(loop)
			loop = loop + 1
			if (phi_i(i) < this%phiMin .or. phi_i(i) > this%phiMax) then
				outBounds = 1
				return
			endif
		enddo

! Hyperpriors
		hyperB_i(1) = this%trial(loop)
		loop = loop + 1
		if (hyperB_i(1) < 0) then
			outBounds = 1
			return
		endif
		
		hyperB_i(2) = this%trial(loop)
		loop = loop + 1
		if (hyperB_i(2) < 0) then
			outBounds = 1
			return
		endif
		
		hypermu_i(1) = this%trial(loop)
		loop = loop + 1
		if (hypermu_i(1) < 0) then
			outBounds = 1
			return
		endif
		
		hypermu_i(2) = this%trial(loop)
		loop = loop + 1
		if (hypermu_i(2) < 0) then
			outBounds = 1
			return
		endif
		
		hyperf_i(1) = this%trial(loop)
		loop = loop + 1
		if (hyperf_i(1) < 0) then
			outBounds = 1
			return
		endif
		
		hyperf_i(2) = this%trial(loop)
		loop = loop + 1
		if (hyperf_i(2) < 0) then
			outBounds = 1
			return
		endif
		
		hyperphi_i(1) = this%trial(loop)
		loop = loop + 1
		if (hyperphi_i(1) < 0) then
			outBounds = 1
			return
		endif
		
		hyperphi_i(2) = this%trial(loop)
		loop = loop + 1
		if (hyperphi_i(2) < 0) then
			outBounds = 1
			return
		endif
		
 		logP = 0.d0
 		nu = 0.001d0
 		
!-----------------
! LOG-PRIORS
!-----------------

! IG(B; a, b) for B
		logP = logP + sum( -(hyperB_i(1)+1.d0) * log(B_i) - hyperB_i(2) / B_i + hyperB_i(1) * log(hyperB_i(2)) - alngam(hyperB_i(1), ierr) )
		
! Beta(mu; a, b, 0, pi) for mu
		logP = logP + sum( (hypermu_i(1)-1.d0) * log(mu_i - this%muMin) + (hypermu_i(2)-1.d0) * log(this%muMax - mu_i) + (1.d0-hypermu_i(1)-hypermu_i(2)) * log(this%muMax - this%muMin) - &
			(alngam(hypermu_i(1),ierr) + alngam(hypermu_i(2),ierr) - alngam(hypermu_i(1)+hypermu_i(2),ierr)) )
			
! Beta(f; a, b, 0, pi) for f
		logP = logP + sum( (hyperf_i(1)-1.d0) * log(f_i - this%fMin) + (hyperf_i(2)-1.d0) * log(this%fMax-f_i) + (1.d0-hyperf_i(1)-hyperf_i(2)) * log(this%fMax - this%fMin) - &
			(alngam(hyperf_i(1),ierr) + alngam(hyperf_i(2),ierr) - alngam(hyperf_i(1)+hyperf_i(2),ierr)) )
			
! Beta(phi; a, b, 0, pi) for phi
		logP = logP + sum( (hyperphi_i(1)-1.d0) * log(phi_i - this%phiMin) + (hyperphi_i(2)-1.d0) * log(this%phiMax-phi_i) + (1.d0-hyperphi_i(1)-hyperphi_i(2)) * log(this%phiMax - this%phiMin) - &
			(alngam(hyperphi_i(1),ierr) + alngam(hyperphi_i(2),ierr) - alngam(hyperphi_i(1)+hyperphi_i(2),ierr)) )
			
! It is a scale parameter so we use a Jeffreys' prior
		logP = logP - (nu-1.d0) * log(hyperB_i(2)) - nu * hyperB_i(2)
		
!-----------------
! DATA LOG-LIKELIHOOD
!-----------------
		c2p = cos(2.d0 * phi_i)
		s2p = sin(2.d0 * phi_i)
		sqrtMu = sqrt(1.d0-mu_i**2)
		
  		logP = logP - 0.5d0 / 1.d-3**2 * (&
			sum( (B_i * mu_i * f_i)**2 * this%CV2 - (B_i * mu_i * f_i) * this%CV3 ) + &
			sum( (B_i * sqrtMu * f_i * c2p)**2 * this%CQ2 - (B_i * sqrtMu * f_i * c2p) * this%CQ3 ) + &
			sum( (B_i * sqrtMu * f_i * s2p)**2 * this%CU2 - (B_i * sqrtMu * f_i * s2p) * this%CU3 ) )
			
		return

	end function logPosteriorWeakField

end module mcmc_class_hierarchical