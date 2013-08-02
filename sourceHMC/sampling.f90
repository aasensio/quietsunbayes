module samplingModule
use varsModule
use mathsModule, only : alngam, digama
implicit none
contains

!------------------------------------------------
! Carry out the sampling using the HMC
!------------------------------------------------
	subroutine doSampling	
	real(kind=8), allocatable :: st(:), stepSize(:), logPGradient(:), logPGradientNew(:), st2(:)
	real(kind=8) :: scaleFactor, logP, logP2
	integer :: seed, fbInt, maxStep, resume, nburn, i
	character(len=128) :: flPfx
	
		scaleFactor = 0.5d0
 		seed = 1234
		fbInt = 10
		maxStep = 10
		resume = 0
		nburn = 0
		flPfx = 'test'
		
		allocate(st(nVariables))
		allocate(stepSize(nVariables))
				
		call initialValuesWeakField(st, stepSize)
				
! Test derivatives
! 		allocate(logPGradient(nVariables))
! 		allocate(logPGradientNew(nVariables))
! 		allocate(st2(nVariables))
! 				
! 		call negLogPosterior(nVariables,st,logP,logPGradient)
! 		do i = 1, nVariables
! 			st2 = st
! 			st2(i) = st(i) + 1.d-6
! 			call negLogPosterior(nVariables,st2,logP2,logPGradientNew)
! 			print *, (logP2 - logP) / 1.d-6, logPGradient(i)
! 		enddo
! 		stop
								
		open(unit=20,file= (trim(flPfx)//".extract"),action='write',status='replace',access='stream')
		
		call run_guided_hmc(nVariables,st,scaleFactor,maxStep,stepSize,flPfx,seed,resume,&
			fbInt, negLogPosterior, writeHMCProcess, nBurn, nSteps)
			
		close(20)		
		
	end subroutine doSampling
	
!------------------------------------------------------------------
! Set initial values for the parameters
!------------------------------------------------------------------
	subroutine initialValuesWeakField(pars, stepSize)
	real(kind=8) :: pars(:), stepSize(:)
	integer :: loop, i, j
	real(kind=8) :: value

		loop = 1
		
! B
		do i = 1, npixels
			pars(loop) = 200.d0
			stepSize(loop) = 5.d0
			loop = loop + 1
		enddo

! mu
		do i = 1, npixels
			pars(loop) = 0.2d0
			stepSize(loop) = 0.05d0
			loop = loop + 1
		enddo

! f
		do i = 1, npixels
			pars(loop) = 0.2d0
			stepSize(loop) = 0.05d0
			loop = loop + 1
		enddo
		
! phi
		do i = 1, npixels
			pars(loop) = 0.1d0
			stepSize(loop) = 0.03d0
			loop = loop + 1
		enddo


! hyperB
		pars(loop) = 2.d0
		stepSize(loop) = 0.2d0
		loop = loop + 1
		
		pars(loop) = 200.d0
		stepSize(loop) = 15.d0
		loop = loop + 1
				
! hypermu
		pars(loop) = 150.5d0
		stepSize(loop) = 15.d0
		loop = loop + 1
		
		pars(loop) = 200.5d0
		stepSize(loop) = 15.d0
		loop = loop + 1
		
! hyperf		
		pars(loop) = 2.5d0
		stepSize(loop) = 1.2d0
		loop = loop + 1
		
		pars(loop) = 3.5d0
		stepSize(loop) = 1.5d0
		loop = loop + 1
		
! hyperphi		
		pars(loop) = 11.5d0
		stepSize(loop) = 4.d0
		loop = loop + 1
		
		pars(loop) = 200.5d0
		stepSize(loop) = 50.d0
		loop = loop + 1
		
	end subroutine initialValuesWeakField

	
!------------------------------------------------
! Negative log-posterior and its gradient
!------------------------------------------------
 	subroutine negLogPosterior(nVariables,trial,logP,logPGradient)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: trial
   real(kind=8), dimension(nVariables) :: logPGradient
   real(kind=8) :: logP, nu, logP2
	integer :: i, outBounds, loop, ierr
		
! Fill the vectors to simplify the notation and test for the boundaries		
		loop = 1				

! B
		do i = 1, nPixels
			B_i(i) = trial(loop)
			loop = loop + 1
			if (B_i(i) < BMin .or. B_i(i) > BMax) then
				logP = -1.d15
				logPGradient = -1.d15
				return
			endif
		enddo

! mu
		do i = 1, npixels
			mu_i(i) = trial(loop)
			loop = loop + 1
			if (mu_i(i) < muMin .or. mu_i(i) > muMax) then
				logP = -1.d15
				logPGradient = -1.d15
				return
			endif
		enddo
	
! f
		do i = 1, npixels
			f_i(i) = trial(loop)
			loop = loop + 1
			if (f_i(i) < fMin .or. f_i(i) > fMax) then
				logP = -1.d15
				logPGradient = -1.d15
				return
			endif
		enddo
	
! phi
		do i = 1, npixels
			phi_i(i) = trial(loop)
			loop = loop + 1
			if (phi_i(i) < phiMin .or. phi_i(i) > phiMax) then
				logP = -1.d15
				logPGradient = -1.d15
				return
			endif
		enddo

! Hyperpriors
		hyperB_i(1) = trial(loop)
		loop = loop + 1
		if (hyperB_i(1) < 0) then
			logP = -1.d15
			logPGradient = -1.d15
			return
		endif
		
		hyperB_i(2) = trial(loop)
		loop = loop + 1
		if (hyperB_i(2) < 0) then
			logP = -1.d15
			logPGradient = -1.d15
			return
		endif
		
		hypermu_i(1) = trial(loop)
		loop = loop + 1
		if (hypermu_i(1) < 0) then
			logP = -1.d15
			logPGradient = -1.d15
			return
		endif
		
		hypermu_i(2) = trial(loop)
		loop = loop + 1
		if (hypermu_i(2) < 0) then
			logP = -1.d15
			logPGradient = -1.d15
			return
		endif
		
		hyperf_i(1) = trial(loop)
		loop = loop + 1
		if (hyperf_i(1) < 0) then
			logP = -1.d15
			logPGradient = -1.d15
			return
		endif
		
		hyperf_i(2) = trial(loop)
		loop = loop + 1
		if (hyperf_i(2) < 0) then
			logP = -1.d15
			logPGradient = -1.d15
			return
		endif
		
		hyperphi_i(1) = trial(loop)
		loop = loop + 1
		if (hyperphi_i(1) < 0) then
			logP = -1.d15
			logPGradient = -1.d15
			return
		endif
		
		hyperphi_i(2) = trial(loop)
		loop = loop + 1
		if (hyperphi_i(2) < 0) then
			logP = -1.d15
			logPGradient = -1.d15
			return
		endif
		
 		logP = 0.d0
		logPGradient = 0.d0
 		nu = 0.001d0
 		
!-----------------
! LOG-PRIORS
!-----------------

!-------------
! IG(B; a, b) for B
!-------------
 		logP = logP + sum( -(hyperB_i(1)+1.d0) * log(B_i) - hyperB_i(2) / B_i + hyperB_i(1) * log(hyperB_i(2)) - alngam(hyperB_i(1), ierr) )
 		
! dIG/dB_i
		logPGradient(1:nPixels) = logPGradient(1:nPixels) + hyperB_i(2) / B_i**2 - (1.d0+hyperB_i(1)) / B_i
! dIG/dalpha
		logPGradient(4*nPixels+1) = sum( log(hyperB_i(2)) - log(B_i) - digama(hyperB_i(1), ierr) )
! dIG/dbeta
		logPGradient(4*nPixels+2) = sum( -1.d0 / B_i + hyperB_i(1) / hyperB_i(2) )
		
!-------------		
! Beta(mu; a, b, 0, pi) for mu
!-------------
		logP = logP + sum( (hypermu_i(1)-1.d0) * log(mu_i - muMin) + (hypermu_i(2)-1.d0) * log(muMax - mu_i) + (1.d0-hypermu_i(1)-hypermu_i(2)) * log(muMax - muMin) - &
			(alngam(hypermu_i(1),ierr) + alngam(hypermu_i(2),ierr) - alngam(hypermu_i(1)+hypermu_i(2),ierr)) )
						
! dBeta/dmu_i
		logPGradient(nPixels+1:2*nPixels) = logPGradient(nPixels+1:2*nPixels) + (1.d0 - hypermu_i(2)) / (muMax - mu_i) + (hypermu_i(1) - 1.d0) / (mu_i - muMin)
! dBeta/dalpha
		logPGradient(4*nPixels+3) = sum( -log(muMax - muMin) + log(mu_i - muMin) - digama(hypermu_i(1), ierr) + digama(hypermu_i(1)+hypermu_i(2), ierr) )
! dBeta/dbeta
		logPGradient(4*nPixels+4) = sum( -log(muMax - muMin) + log(muMax - mu_i) - digama(hypermu_i(2), ierr) + digama(hypermu_i(1)+hypermu_i(2), ierr) )
		

!-------------
! Beta(f; a, b, 0, pi) for f
!-------------
		logP = logP + sum( (hyperf_i(1)-1.d0) * log(f_i - fMin) + (hyperf_i(2)-1.d0) * log(fMax-f_i) + (1.d0-hyperf_i(1)-hyperf_i(2)) * log(fMax - fMin) - &
			(alngam(hyperf_i(1),ierr) + alngam(hyperf_i(2),ierr) - alngam(hyperf_i(1)+hyperf_i(2),ierr)) )
! dBeta/df_i
		logPGradient(2*nPixels+1:3*nPixels) = logPGradient(2*nPixels+1:3*nPixels) + (1.d0 - hyperf_i(2)) / (fMax - f_i) + (hyperf_i(1) - 1.d0) / (f_i - fMin)
! dBeta/dalpha
		logPGradient(4*nPixels+5) = sum( -log(fMax - fMin) + log(f_i - fMin) - digama(hyperf_i(1), ierr) + digama(hyperf_i(1)+hyperf_i(2), ierr) )
! dBeta/dbeta
		logPGradient(4*nPixels+6) = sum( -log(fMax - fMin) + log(fMax - f_i) - digama(hyperf_i(2), ierr) + digama(hyperf_i(1)+hyperf_i(2), ierr) )
			
!-------------
! Beta(phi; a, b, 0, pi) for phi
!-------------
		logP = logP + sum( (hyperphi_i(1)-1.d0) * log(phi_i - phiMin) + (hyperphi_i(2)-1.d0) * log(phiMax-phi_i) + (1.d0-hyperphi_i(1)-hyperphi_i(2)) * log(phiMax - phiMin) - &
			(alngam(hyperphi_i(1),ierr) + alngam(hyperphi_i(2),ierr) - alngam(hyperphi_i(1)+hyperphi_i(2),ierr)) )
			
! dBeta/dphi_i
		logPGradient(3*nPixels+1:4*nPixels) = logPGradient(3*nPixels+1:4*nPixels) + (1.d0 - hyperphi_i(2)) / (phiMax - phi_i) + (hyperphi_i(1) - 1.d0) / (phi_i - phiMin)
! dBeta/dalpha
		logPGradient(4*nPixels+7) = sum( -log(phiMax - phiMin) + log(phi_i - phiMin) - digama(hyperphi_i(1), ierr) + digama(hyperphi_i(1)+hyperphi_i(2), ierr) )
! dBeta/dbeta
		logPGradient(4*nPixels+8) = sum( -log(phiMax - phiMin) + log(phiMax - phi_i) - digama(hyperphi_i(2), ierr) + digama(hyperphi_i(1)+hyperphi_i(2), ierr) )
			
!-------------
! It is a scale parameter so we use a Jeffreys' prior
!-------------
		logP = logP - (nu-1.d0) * log(hyperB_i(2)) - nu * hyperB_i(2)
		logPGradient(4*nPixels+2) = logPGradient(4*nPixels+2) + (1.d0 - nu) / hyperB_i(2) - nu
		
!-----------------
! DATA LOG-LIKELIHOOD
!-----------------
		c2p = cos(2.d0 * phi_i)
		s2p = sin(2.d0 * phi_i)
		sqrtMu = sqrt(1.d0-mu_i**2)
		dsqrtMudMu = -mu_i / sqrtMu
		
  		logP = logP - 0.5d0 / 1.d-3**2 * (&
			sum( CV1 + (B_i * mu_i * f_i)**2 * CV2 - (B_i * mu_i * f_i) * CV3 ) + &
			sum( CQ1 + (B_i * sqrtMu * f_i * c2p)**2 * CQ2 - (B_i * sqrtMu * f_i * c2p) * CQ3 ) + &
			sum( CU1 + (B_i * sqrtMu * f_i * s2p)**2 * CU2 - (B_i * sqrtMu * f_i * s2p) * CU3 ) )
			
! dlogL/dB
		logPGradient(1:nPixels) = logPGradient(1:nPixels) - 0.5d0 / 1.d-3**2 * (&
			( 2.d0 * B_i * (mu_i * f_i)**2 * CV2 - (mu_i * f_i) * CV3 ) + &
			( 2.d0 * B_i * (sqrtMu * f_i * c2p)**2 * CQ2 - (sqrtMu * f_i * c2p) * CQ3 ) + &
			( 2.d0 * B_i * (sqrtMu * f_i * s2p)**2 * CU2 - (sqrtMu * f_i * s2p) * CU3 ) )
			
! dlogL/dmu
		logPGradient(nPixels+1:2*nPixels) = logPGradient(nPixels+1:2*nPixels) - 0.5d0 / 1.d-3**2 * (&
			( 2.d0 * mu_i * (B_i * f_i)**2 * CV2 - (B_i * f_i) * CV3 ) + &
			dsqrtMudMu * ( 2.d0 * sqrtMu * (B_i * f_i * c2p)**2 * CQ2 - (B_i * f_i * c2p) * CQ3 ) + &
			dsqrtMudMu * ( 2.d0 * sqrtMu * (B_i * f_i * s2p)**2 * CU2 - (B_i * f_i * s2p) * CU3 ) )
			
! dlogL/df
		logPGradient(2*nPixels+1:3*nPixels) = logPGradient(2*nPixels+1:3*nPixels) - 0.5d0 / 1.d-3**2 * (&
			( 2.d0 * f_i * (B_i * mu_i)**2 * CV2 - (B_i * mu_i) * CV3 ) + &
			( 2.d0 * f_i * (B_i * sqrtMu * c2p)**2 * CQ2 - (B_i * sqrtMu * c2p) * CQ3 ) + &
			( 2.d0 * f_i * (B_i * sqrtMu * s2p)**2 * CU2 - (B_i * sqrtMu * s2p) * CU3 ) )
			
! dlogL/dphi
		logPGradient(3*nPixels+1:4*nPixels) = logPGradient(3*nPixels+1:4*nPixels) - 0.5d0 / 1.d-3**2 * (&			
			( -4.d0 * c2p * s2p * (B_i * sqrtMu * f_i)**2 * CQ2 + 2.d0 * s2p * (B_i * sqrtMu * f_i) * CQ3 ) + &
			( 4.d0 * s2p * c2p * (B_i * sqrtMu * f_i)**2 * CU2 - 2.d0 * c2p * (B_i * sqrtMu * f_i) * CU3 ) )
					
		logP = -logP
		logPGradient = -logPGradient
				  
	end subroutine negLogPosterior

!------------------------------------------------
! A subroutine to write the extract file
! I have assumed that the unit=20 is opened for
! writing (append) earlier. In general only write
! those parametes which are estimated. The files
! can be really big depending on the dimensionality
!------------------------------------------------
	subroutine writeHMCProcess(nVariables,x,v,g)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: x
   real(kind=8), dimension(nVariables) :: g
   real(kind=8) :: v, xwrite(8)
	integer i
		
	xwrite = x(4*nPixels+1:nVariables)
	write(20) xwrite
				
	end subroutine writeHMCProcess

end module samplingModule