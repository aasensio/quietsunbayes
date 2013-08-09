module samplingModule
use varsModule
use likelihoodModule, only : negLogPosterior, writeHMCProcess, writeHMCProcessBurnin
use maximumlikelihoodModule, only : getMaximumLikelihood
implicit none
contains

!------------------------------------------------
! Carry out the sampling using the HMC
!------------------------------------------------
	subroutine doSampling(action)
	character(len=100) :: action
	real(kind=8), allocatable :: st(:), stepSize(:), logPGradient(:), logPGradientNew(:), st2(:), meanOld(:)
	real(kind=8) :: scaleFactor, logP, logP2
	integer :: seed, fbInt, maxStep, resume, nburn, i, nStepsBurn
	character(len=128) :: flPfx
	
		scaleFactor = 0.5d0
		seed = 1234
		fbInt = 10
		maxStep = 1
		if (action(1:4) == 'CONT') then
			resume = 1
		else
			resume = 0
		endif
		nburn = 0
		flPfx = 'test            '
		
		allocate(st(nVariables))
		allocate(stepSize(nVariables))
				
		if (resume == 0) then
			call initialValuesWeakField(st, stepSize)
!  			call getMaximumLikelihood(st)
								
			parsOld = st		
			nStep = 1
			parsVariance = 0.d0
			parsMean = 0.d0
					
! Test derivatives
! 	 		allocate(logPGradient(nVariables))
! 	 		allocate(logPGradientNew(nVariables))
! 	 		allocate(st2(nVariables))
! 	 				
! 	 		call negLogPosterior(nVariables,st,logP,logPGradient)
! 			do i = 1, nVariables
! 				st2 = st
! 				st2(i) = st(i) + 1.d-6
! 				call negLogPosterior(nVariables,st2,logP2,logPGradientNew)
! 				print *, (logP2 - logP) / 1.d-6, logPGradient(i)
! 			enddo
! 			stop

			open(unit=20,file= flPfx(1:len_trim(flPfx))//".burnin",action='write',status='replace',access='stream')		
			nStepsBurn = 200
			call run_guided_hmc(nVariables,st,scaleFactor,maxStep,stepSize,flPfx(1:len_trim(flPfx)),seed,resume,&
				fbInt, negLogPosterior, writeHMCProcessBurnin, nBurn, nStepsBurn)			
			close(20)
			
	! Estimate width of distributions
			open(unit=20,file=flPfx(1:len_trim(flPfx))//".burnin",action='read',status='old',access='stream')
			parsMean = 0.d0
			parsVariance = 0.d0
			nStep = 1
			allocate(meanOld(nVariables))
			do i = 1, nStepsBurn
				read(20) st
				if (i > 100) then
					meanOld = parsMean
					parsMean = meanOld + (st - meanOld) / (nStep + 1.d0)		
					parsVariance = (nStep - 1.d0) / nStep * parsVariance + (st - meanOld)**2 / (nStep+1.d0)**2 + (st - meanOld)**2 / nStep
					nStep = nStep + 1
				endif
			enddo		
			deallocate(meanOld)
			
			close(20)
			
			stepSize = sqrt(parsVariance)
			st = parsMean
			
			open(unit=20,file='variances',action='write',status='replace')
			do i = 1, nVariables
				write(20,*) parsMean(i), stepSize(i)
			enddo
			close(20)
		
		endif
		
		maxStep = 10
		
		if (resume == 0) then					
			open(unit=20,file= (trim(flPfx)//".extract"),action='write',status='replace',access='stream')
			open(unit=21,file= (trim(flPfx)//".extract2"),action='write',status='replace',access='stream')
		else
			open(unit=20,file= (trim(flPfx)//".extract"),action='write',position='append',access='stream')
		endif
		
		call run_guided_hmc(nVariables,st,scaleFactor,maxStep,stepSize,flPfx(1:len_trim(flPfx)),seed,resume,&
			fbInt, negLogPosterior, writeHMCProcess, nBurn, nSteps)
		close(20)
		
		open(unit=20,file='posterior.sizes',action='write',status='replace')
		write(20,*) 6, nSteps
		close(20)
		
	end subroutine doSampling
	
!------------------------------------------------------------------
! Set initial values for the parameters
!------------------------------------------------------------------
	subroutine initialValuesWeakField(pars, stepSize)
	real(kind=8) :: pars(:), stepSize(:)
	integer :: loop, i, j
	real(kind=8) :: value
	real(kind=8), allocatable :: Bpar(:), Bperp(:), azimuth(:)

		loop = 1
		
		allocate(Bpar(npixels))
		allocate(Bperp(npixels))
		allocate(azimuth(npixels))
		
! Compute maximum-likelihood solution
		Bpar = 0.5d0 * CV3 / CV2
		Bperp = sqrt(0.5d0 * sqrt(CQ3**2 + CU3**2) / CQ2)
		
		do i = 1, npixels
			azimuth(i) = 0.5d0 * atan2(CU3(i), CQ3(i))
 			if (azimuth(i) < 0.d0) azimuth(i) = azimuth(i) + PI
		enddo
					
! B
		do i = 1, npixels
			pars(loop) = sqrt(Bpar(i)**2 + Bperp(i)**2) / 0.5
			if (pars(loop) > Bmax) pars(loop) = Bmax - 200.d0
			if (pars(loop) < BMin) pars(loop) = BMin + 200.d0
			stepSize(loop) = 40.d0			
			loop = loop + 1
		enddo		

! mu
		do i = 1, npixels
			pars(loop) = cos(atan2(Bperp(i),Bpar(i)))
			stepSize(loop) = 0.1d0
			loop = loop + 1
		enddo				
		
! f
		do i = 1, npixels
			pars(loop) = 0.5d0
			stepSize(loop) = 0.3d0
			loop = loop + 1
		enddo
		
! phi
		do i = 1, npixels
			pars(loop) = azimuth(i)
			stepSize(loop) = 2.0d0
			loop = loop + 1
		enddo
								
	end subroutine initialValuesWeakField
	
end module samplingModule