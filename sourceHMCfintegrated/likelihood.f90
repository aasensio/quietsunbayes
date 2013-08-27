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
		B_i = sigmoid(trial(1:nPixels), BMin, BMax)

! mu
		mu_i = sigmoid(trial(nPixels+1:2*nPixels), muMin, muMax)
		
! phi
		phi_i = sigmoid(trial(2*nPixels+1:3*nPixels), phiMin, phiMax)			

! Hyperpriors
		hyperB_i(1) = sigmoid(trial(3*nPixels+1), 0.d0, hyperparRanges(1,1))
		hyperB_i(2) = sigmoid(trial(3*nPixels+2), 0.d0, hyperparRanges(2,1))
		
		hypermu_i(1) = sigmoid(trial(3*nPixels+3), 0.d0, hyperparRanges(1,2))
		hypermu_i(2) = sigmoid(trial(3*nPixels+4), 0.d0, hyperparRanges(2,2))
												
 		logP = 0.d0
		logPGradient = 0.d0
 		nu = 0.001d0
 		
!-----------------
! LOG-PRIORS
!-----------------

!-------------
! IG(B; a, b) for B
!-------------
 		logP = logP + sum( -(hyperB_i(1)+1.d0) * log(B_i) - hyperB_i(2) / B_i ) + nPixels * hyperB_i(1) * log(hyperB_i(2)) - nPixels * alngam(hyperB_i(1), ierr)
 		
! dIG/dB_i
		logPGradient(1:nPixels) = logPGradient(1:nPixels) + hyperB_i(2) / B_i**2 - (1.d0+hyperB_i(1)) / B_i
! dIG/dalpha
		logPGradient(3*nPixels+1) = nPixels * log(hyperB_i(2)) - sum(log(B_i)) - nPixels * digama(hyperB_i(1), ierr)
! dIG/dbeta
		logPGradient(3*nPixels+2) = sum(-1.d0 / B_i) + nPixels * hyperB_i(1) / hyperB_i(2)
		
!-------------		
! Beta(mu; a, b, 0, pi) for mu
!-------------
		logP = logP + sum( (hypermu_i(1)-1.d0) * log(mu_i - muMin) + (hypermu_i(2)-1.d0) * log(muMax - mu_i) ) + nPixels * (1.d0-hypermu_i(1)-hypermu_i(2)) * log(muMax - muMin) - &
			nPixels * (alngam(hypermu_i(1),ierr) + alngam(hypermu_i(2),ierr) - alngam(hypermu_i(1)+hypermu_i(2),ierr))
						
! dBeta/dmu_i
		logPGradient(nPixels+1:2*nPixels) = logPGradient(nPixels+1:2*nPixels) + (1.d0 - hypermu_i(2)) / (muMax - mu_i) + (hypermu_i(1) - 1.d0) / (mu_i - muMin)
! dBeta/dalpha
		logPGradient(3*nPixels+3) = -nPixels * log(muMax - muMin) + sum( log(mu_i - muMin) ) - nPixels * digama(hypermu_i(1), ierr) + nPixels * digama(hypermu_i(1)+hypermu_i(2), ierr)
! dBeta/dbeta
		logPGradient(3*nPixels+4) = -nPixels * log(muMax - muMin) + sum( log(muMax - mu_i) ) - nPixels * digama(hypermu_i(2), ierr) + nPixels * digama(hypermu_i(1)+hypermu_i(2), ierr)
		

!-------------
! It is a scale parameter so we use a Jeffreys' prior
!-------------
		logP = logP - (nu-1.d0) * log(hyperB_i(2)) - nu * hyperB_i(2)
		logPGradient(3*nPixels+2) = logPGradient(3*nPixels+2) + (1.d0 - nu) / hyperB_i(2) - nu
		
!-------------
! Hyperpriors for alpha and beta for each Beta prior		
!-------------
		logP = logP - 2.5d0 * log(sum(hypermu_i))
		logPGradient(3*nPixels+3) = logPGradient(3*nPixels+3) - 2.5d0 / sum(hypermu_i)
		logPGradient(3*nPixels+4) = logPGradient(3*nPixels+4) - 2.5d0 / sum(hypermu_i)
					
!-----------------
! DATA LOG-LIKELIHOOD
!-----------------
		c2p = cos(2.d0 * phi_i)
		s2p = sin(2.d0 * phi_i)
		A = 0.5d0 / sigma_n**2 * (CV1 + CQ1 + CU1)
		B = 0.5d0 / sigma_n**2 * (B_i * mu_i * CV3 + B_i**2 * (1.d0-mu_i**2) * (c2p*CQ3+s2p*CU3) )
		C = 0.5d0 / sigma_n**2 * ( (B_i * mu_i)**2 * CV2 + (B_i**2 * (1.d0-mu_i**2))**2 * (c2p**2 * CQ2 + s2p**2 * CU2) )
						
		dlogPdB = B / (2.d0*C) + (exp(-0.25d0 * B**2/C) - exp(-0.25d0*(B-2.d0*C)**2 / C) ) / (erf(0.5d0*B/sqrt(C)) - erf(0.5d0*(B-2.d0*C) / sqrt(C) )) / (sqrt(PI*C))
		dlogPdC = -0.5d0/C - 0.25d0 * (B/C)**2 - ( exp(-0.25d0*B**2/C) / (2.d0*sqrt(PI)*C**1.5) + 2.d0/sqrt(PI) * exp(-0.25d0*(B-2.d0*C)**2/C) * (-1.d0/sqrt(C)-0.25d0*(B-2.d0*C)/C**1.5) )/ &
			(erf(0.5d0*B/sqrt(C)) - erf(0.5d0*(B-2.d0*C) / sqrt(C) ))
				
		dBdB_i = 0.5d0 / sigma_n**2 * (mu_i*CV3 + 2.d0*B_i*(1.d0-mu_i**2) * (c2p*CQ3 + s2p*CU3))
		dCdB_i = 0.5d0 / sigma_n**2 * (2.d0*B_i*mu_i**2*CV2 + 4.d0*B_i**3*(1.d0-mu_i**2)**2* (c2p**2*CQ2 + s2p**2*CU2))
				
		dBdmu_i = 0.5d0 / sigma_n**2 * (B_i*CV3 - 2.d0*mu_i*B_i**2 * (c2p*CQ3 + s2p*CU3))
		dCdmu_i = 0.5d0 / sigma_n**2 * (2.d0*mu_i*B_i**2*CV2 - 4.d0*mu_i*(1.d0-mu_i**2)*B_i**4 * (c2p**2*CQ2 + s2p**2*CU2))
		
		dBdphi_i = 0.5d0 / sigma_n**2 * (2.d0*B_i**2*(1.d0-mu_i**2) * (-s2p*CQ3 + c2p*CU3))
		dCdphi_i = 0.5d0 / sigma_n**2 * (4.d0*B_i**4*(1.d0-mu_i**2)**2 * c2p * s2p * (-CQ2 + CU2))
						
  		logP = logP - 0.5d0 * sum(log(C)) - sum(A - B**2 / (4.d0*C)) + sum( log( erf(0.5d0*B/sqrt(C)) - erf(0.5d0*(B-2.d0*C) / sqrt(C) ) ) )
			
! dlogL/dB
		logPGradient(1:nPixels) = logPGradient(1:nPixels) + dlogPdB * dBdB_i + dlogPdC * dCdB_i
					
! dlogL/dmu
		logPGradient(nPixels+1:2*nPixels) = logPGradient(nPixels+1:2*nPixels) +  + dlogPdB * dBdmu_i + dlogPdC * dCdmu_i
						
! dlogL/dphi
		logPGradient(2*nPixels+1:3*nPixels) = logPGradient(2*nPixels+1:3*nPixels) + dlogPdB * dBdphi_i + dlogPdC * dCdphi_i
			
!-----------------
! SIGMOID TRANSFORMATION FOR BOUNDED PARAMETERS
!-----------------			
		logPGradient(1:nPixels) = logPGradient(1:nPixels) * diffSigmoid(trial(1:nPixels), BMin, BMax)
		logPGradient(nPixels+1:2*nPixels) = logPGradient(nPixels+1:2*nPixels) * diffSigmoid(trial(nPixels+1:2*nPixels), muMin, muMax)		
		logPGradient(2*nPixels+1:3*nPixels) = logPGradient(2*nPixels+1:3*nPixels) * diffSigmoid(trial(2*nPixels+1:3*nPixels), phiMin, phiMax)
		
		logPGradient(3*nPixels+1) = logPGradient(3*nPixels+1) * diffSigmoid(trial(3*nPixels+1), 0.d0, hyperparRanges(1,1))
		logPGradient(3*nPixels+2) = logPGradient(3*nPixels+2) * diffSigmoid(trial(3*nPixels+2), 0.d0, hyperparRanges(2,1))
		
		logPGradient(3*nPixels+3) = logPGradient(3*nPixels+3) * diffSigmoid(trial(3*nPixels+3), 0.d0, hyperparRanges(1,2))
		logPGradient(3*nPixels+4) = logPGradient(3*nPixels+4) * diffSigmoid(trial(3*nPixels+4), 0.d0, hyperparRanges(2,2))
					
					
		logP = -logP
		logPGradient = -logPGradient
		
!   		print *, f_i(1:2), trial(2*nPixels+1:2*nPixels+2)
								  
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
   real(kind=8) :: v, xwrite(4)
	integer i
			
		meanOld = parsMean
		parsMean = meanOld + (x - meanOld) / (nStep + 1.d0)		
		parsVariance = (nStep - 1.d0) / nStep * parsVariance + (x - meanOld)**2 / (nStep+1.d0)**2 + (x - meanOld)**2 / nStep
					
		x2 = x
		x2(1:nPixels) = sigmoid(x(1:nPixels), BMin, BMax)
		x2(nPixels+1:2*nPixels) = sigmoid(x(nPixels+1:2*nPixels), muMin, muMax)		
		x2(2*nPixels+1:3*nPixels) = sigmoid(x(2*nPixels+1:3*nPixels), phiMin, phiMax)
		
		x2(3*nPixels+1) = sigmoid(x(3*nPixels+1), 0.d0, hyperparRanges(1,1))
		x2(3*nPixels+2) = sigmoid(x(3*nPixels+2), 0.d0, hyperparRanges(2,1))
		x2(3*nPixels+3) = sigmoid(x(3*nPixels+3), 0.d0, hyperparRanges(1,2))
		x2(3*nPixels+4) = sigmoid(x(3*nPixels+4), 0.d0, hyperparRanges(2,2))
				
		xwrite = x2(3*nPixels+1:nVariables)
		
		write(20) xwrite
		
! If we are in the last steps of the chain, write them to a file for continuation of the chain
		if ( (nSteps - nStep) < nStepsBurn ) then
			write(21) x
		endif
				
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