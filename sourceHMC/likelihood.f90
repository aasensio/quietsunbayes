module likelihoodModule
use varsModule
use mathsModule, only : alngam, digama, sigmoid, diffSigmoid, invSigmoid
implicit none

contains
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
		B_i = trial(1:nPixels)
		if (any(B_i < BMin .or. B_i > BMax)) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif		

! mu
		mu_i = sigmoid(trial(nPixels+1:2*nPixels), muMin, muMax)
		if (any(mu_i < muMin .or. mu_i > muMax)) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif
	
! f
		f_i = sigmoid(trial(2*nPixels+1:3*nPixels), fMin, fMax)
		if (any(f_i < fMin .or. f_i > fMax)) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif
	
! phi
		phi_i = sigmoid(trial(3*nPixels+1:4*nPixels), phiMin, phiMax)
		if (any(phi_i < phiMin .or. phi_i > phiMax)) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif

! Hyperpriors
		loop = 4*nPixels+1
		hyperB_i(1) = trial(loop)
		loop = loop + 1
		if (hyperB_i(1) < 0) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif
		
		hyperB_i(2) = trial(loop)
		loop = loop + 1
		if (hyperB_i(2) < 0) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif
		
		hypermu_i(1) = trial(loop)
		loop = loop + 1
		if (hypermu_i(1) < 0) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif
		
		hypermu_i(2) = trial(loop)
		loop = loop + 1
		if (hypermu_i(2) < 0) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif
		
		hyperf_i(1) = trial(loop)
		loop = loop + 1
		if (hyperf_i(1) < 0) then
			logP = 1.d15
			logPGradient = 1.d15
			return
		endif
		
		hyperf_i(2) = trial(loop)
		loop = loop + 1
		if (hyperf_i(2) < 0) then
			logP = 1.d15
			logPGradient = 1.d15
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
! It is a scale parameter so we use a Jeffreys' prior
!-------------
		logP = logP - (nu-1.d0) * log(hyperB_i(2)) - nu * hyperB_i(2)
		logPGradient(4*nPixels+2) = logPGradient(4*nPixels+2) + (1.d0 - nu) / hyperB_i(2) - nu
		
!-------------
! Hyperpriors for alpha and beta for each Beta prior		
!-------------
		logP = logP - 2.5d0 * log(sum(hypermu_i))
		logPGradient(4*nPixels+3) = logPGradient(4*nPixels+3) - 2.5d0 / sum(hypermu_i)
		logPGradient(4*nPixels+4) = logPGradient(4*nPixels+4) - 2.5d0 / sum(hypermu_i)
		
		logP = logP - 2.5d0 * log(sum(hyperf_i))
		logPGradient(4*nPixels+5) = logPGradient(4*nPixels+5) - 2.5d0 / sum(hyperf_i)
		logPGradient(4*nPixels+6) = logPGradient(4*nPixels+6) - 2.5d0 / sum(hyperf_i)
		
		
!-----------------
! DATA LOG-LIKELIHOOD
!-----------------
		c2p = cos(2.d0 * phi_i)
		s2p = sin(2.d0 * phi_i)
		sqrtMu = sqrt(1.d0-mu_i**2)
		dsqrtMudMu = -mu_i / sqrtMu
		
  		logP = logP - 0.5d0 / sigma_n**2 * (&
			sum( CV1 + (B_i * mu_i * f_i)**2 * CV2 - (B_i * mu_i * f_i) * CV3 ) + &
			sum( CQ1 + (B_i**2 * (1.d0-mu_i**2) * f_i * c2p)**2 * CQ2 - (B_i**2 * (1.d0-mu_i**2) * f_i * c2p) * CQ3 ) + &
			sum( CU1 + (B_i**2 * (1.d0-mu_i**2) * f_i * s2p)**2 * CU2 - (B_i**2 * (1.d0-mu_i**2) * f_i * s2p) * CU3 ) )
			
! dlogL/dB
		logPGradient(1:nPixels) = logPGradient(1:nPixels) - 0.5d0 / sigma_n**2 * (&
			( 2.d0 * B_i * (mu_i * f_i)**2 * CV2 - (mu_i * f_i) * CV3 ) + &
			( 4.d0 * B_i**3 * ((1.d0-mu_i**2) * f_i * c2p)**2 * CQ2 - 2.d0 * B_i * ((1.d0-mu_i**2) * f_i * c2p) * CQ3 ) + &
			( 4.d0 * B_i**3 * ((1.d0-mu_i**2) * f_i * s2p)**2 * CU2 - 2.d0 * B_i * ((1.d0-mu_i**2) * f_i * s2p) * CU3 ) )
			
! dlogL/dmu
		logPGradient(nPixels+1:2*nPixels) = logPGradient(nPixels+1:2*nPixels) - 0.5d0 / sigma_n**2 * (&
			( 2.d0 * mu_i * (B_i * f_i)**2 * CV2 - (B_i * f_i) * CV3 ) + &
			( 2.d0 * (1.d0-mu_i**2) * (-2.d0*mu_i) * (B_i**2 * f_i * c2p)**2 * CQ2 - (-2.d0*mu_i) * (B_i**2 * f_i * c2p) * CQ3 ) + &
			( 2.d0 * (1.d0-mu_i**2) * (-2.d0*mu_i) * (B_i**2 * f_i * s2p)**2 * CU2 - (-2.d0*mu_i) * (B_i**2 * f_i * s2p) * CU3 ) )
			
! dlogL/df
		logPGradient(2*nPixels+1:3*nPixels) = logPGradient(2*nPixels+1:3*nPixels) - 0.5d0 / sigma_n**2 * (&
			( 2.d0 * f_i * (B_i * mu_i)**2 * CV2 - (B_i * mu_i) * CV3 ) + &
			( 2.d0 * f_i * (B_i**2 * (1.d0-mu_i**2) * c2p)**2 * CQ2 - (B_i**2 * (1.d0-mu_i**2) * c2p) * CQ3 ) + &
			( 2.d0 * f_i * (B_i**2 * (1.d0-mu_i**2) * s2p)**2 * CU2 - (B_i**2 * (1.d0-mu_i**2) * s2p) * CU3 ) )
			
! dlogL/dphi
		logPGradient(3*nPixels+1:4*nPixels) = logPGradient(3*nPixels+1:4*nPixels) - 0.5d0 / sigma_n**2 * (&			
			( -4.d0 * c2p * s2p * (B_i**2 * (1.d0-mu_i**2) * f_i)**2 * CQ2 + 2.d0 * s2p * (B_i**2 * (1.d0-mu_i**2) * f_i) * CQ3 ) + &
			( 4.d0 * s2p * c2p * (B_i**2 * (1.d0-mu_i**2) * f_i)**2 * CU2 - 2.d0 * c2p * (B_i**2 * (1.d0-mu_i**2) * f_i) * CU3 ) )
			
!-----------------
! SIGMOID TRANSFORMATION FOR BOUNDED PARAMETERS
!-----------------		
		logPGradient(nPixels+1:2*nPixels) = logPGradient(nPixels+1:2*nPixels) * diffSigmoid(trial(nPixels+1:2*nPixels), muMin, muMax)
		logPGradient(2*nPixels+1:3*nPixels) = logPGradient(2*nPixels+1:3*nPixels) * diffSigmoid(trial(2*nPixels+1:3*nPixels), fMin, fMax)
		logPGradient(3*nPixels+1:4*nPixels) = logPGradient(3*nPixels+1:4*nPixels) * diffSigmoid(trial(3*nPixels+1:4*nPixels), phiMin, phiMax)
					
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
	real(kind=8), dimension(nVariables) :: x, x2
   real(kind=8), dimension(nVariables) :: g, meanOld
   real(kind=8) :: v, xwrite(6)
	integer i
			
		meanOld = parsMean
		parsMean = meanOld + (x - meanOld) / (nStep + 1.d0)		
		parsVariance = (nStep - 1.d0) / nStep * parsVariance + (x - meanOld)**2 / (nStep+1.d0)**2 + (x - meanOld)**2 / nStep
		
		xwrite = x(4*nPixels+1:nVariables)
		write(20) xwrite
		
		x2 = x
		x2(nPixels+1:2*nPixels) = sigmoid(x(nPixels+1:2*nPixels), muMin, muMax)
		x2(2*nPixels+1:3*nPixels) = sigmoid(x(2*nPixels+1:3*nPixels), fMin, fMax)
		x2(3*nPixels+1:4*nPixels) = sigmoid(x(3*nPixels+1:4*nPixels), phiMin, phiMax)
		
		write(21) x2
		
		open(unit=25,file= "test.stddev",action='write',status='replace')
		write(25,*) sqrt(parsVariance)
		close(25)
		
		open(unit=25,file= "test.mean",action='write',status='replace')
		write(25,*) parsMean
		close(25)
		
		nStep = nStep + 1.d0
				
	end subroutine writeHMCProcess
	
!------------------------------------------------
! A subroutine to write the extract file
! I have assumed that the unit=20 is opened for
! writing (append) earlier. In general only write
! those parametes which are estimated. The files
! can be really big depending on the dimensionality
!------------------------------------------------
	subroutine writeHMCProcessBurnin(nVariables,x,v,g)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: x
   real(kind=8), dimension(nVariables) :: g, meanOld
   real(kind=8) :: v, xwrite(8)
	integer i
							
		write(20) x		
						
	end subroutine writeHMCProcessBurnin
		
end module likelihoodModule