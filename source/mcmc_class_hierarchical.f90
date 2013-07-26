module mcmc_class_hierarchical
use mcmc_class
implicit none

	type pixelType		
		real(kind=8) :: noise
	end type pixelType
	
	type, EXTENDS (mcmc) :: mcmc_hierarchical

		integer :: npixels
		type(pixelType), pointer :: pixel(:)
		real(kind=8) :: BMin, BMax, BSigma2, muMin, muMax, muSigma2, fMin, fMax, fSigma2, phiMin, phiMax, phiSigma2
		integer :: BNodes, muNodes, fNodes, phiNodes		
		real(kind=8), pointer :: BLoc(:), muLoc(:), fLoc(:), phiLoc(:)
		
		real(kind=8), pointer, dimension(:) :: CV2, CV3, CQ2, CQ3, CU2, CU3
	
 		contains
 			procedure :: evalPosterior => logPosteriorWeakField
 			procedure :: initialValues => initialValuesWeakField
	end type mcmc_hierarchical
	
	real(kind=8), allocatable :: B_i(:), mu_i(:), phi_i(:), f_i(:), thetaB_i(:), thetamu_i(:), thetaphi_i(:), thetaf_i(:)
		
	contains
	
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

! thetaB
		do i = 1, this%BNodes
			this%pars(loop) = 1.d0
			loop = loop + 1
		enddo
		
! thetamu
		do i = 1, this%muNodes
			this%pars(loop) = 1.d0
			loop = loop + 1
		enddo
		
! thetaf
		do i = 1, this%fNodes
			this%pars(loop) = 1.d0
			loop = loop + 1
		enddo
		
! thetaphi
		do i = 1, this%phiNodes
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
	real(kind=8) :: logP
	integer :: i, j, loop, ierr

! The vector of trial parameters enters with the following order
! B(1..npixels)
! mu(1..npixels)
! f(1..npixels)
! phi(1..npixels)

! thetaB(1..BNodes)
! thetamu(1..muNodes)
! thetaf(1..fNodes)
! thetaphi(1..phiNodes)


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
		do i = 1, this%BNodes
			thetaB_i(i) = this%trial(loop)
			loop = loop + 1
			if (thetaB_i(i) < 0) then
				outBounds = 1
				return
			endif
		enddo
		
		do i = 1, this%muNodes
			thetamu_i(i) = this%trial(loop)
			loop = loop + 1
			if (thetamu_i(i) < 0) then
				outBounds = 1
				return
			endif
		enddo
		
		do i = 1, this%fNodes
			thetaf_i(i) = this%trial(loop)
			loop = loop + 1
			if (thetaf_i(i) < 0) then
				outBounds = 1
				return
			endif
		enddo
		
		do i = 1, this%phiNodes
			thetaphi_i(i) = this%trial(loop)
			loop = loop + 1
			if (thetaphi_i(i) < 0) then
				outBounds = 1
				return
			endif
		enddo


 		logP = 0.d0
 		
!-----------------
! LOG-PRIORS
!-----------------

! Normalize the weights to add to one
		thetaB_i = log(thetaB_i / sum(thetaB_i))
		thetamu_i = log(thetamu_i / sum(thetamu_i))
		thetaf_i = log(thetaf_i / sum(thetaf_i))
		thetaphi_i = log(thetaphi_i / sum(thetaphi_i))
		
		do i = 1, this%npixels
			
			logP = logP + logSumVector( thetaB_i - 0.5d0 * (B_i(i) - this%BLoc)**2 / this%BSigma2 )
			logP = logP + logSumVector( thetamu_i - 0.5d0 * (mu_i(i) - this%muLoc)**2 / this%muSigma2 )
			logP = logP + logSumVector( thetaf_i - 0.5d0 * (f_i(i) - this%fLoc)**2 / this%fSigma2 )
			logP = logP + logSumVector( thetaphi_i - 0.5d0 * (phi_i(i) - this%phiLoc)**2 / this%phiSigma2 )
		
		enddo
		
!-----------------
! DATA LOG-LIKELIHOOD
!-----------------
  		logP = sum((B_i * mu_i * f_i)**2 * this%CV2 - (B_i * mu_i * f_i) * this%CV3) + &
			sum((B_i * sqrt(1.d0-mu_i**2) * f_i * cos(phi_i))**2 * this%CQ2 - (B_i * sqrt(1.d0-mu_i**2) * f_i * cos(phi_i)) * this%CQ3) + &
			sum((B_i * sqrt(1.d0-mu_i**2) * f_i * sin(phi_i))**2 * this%CU2 - (B_i * sqrt(1.d0-mu_i**2) * f_i * sin(phi_i)) * this%CU3)
			
		return

	end function logPosteriorWeakField

end module mcmc_class_hierarchical