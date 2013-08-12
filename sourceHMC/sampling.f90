module samplingModule
use varsModule
use likelihoodModule, only : negLogPosterior, writeHMCProcess, writeHMCProcessBurnin
use maximumlikelihoodModule, only : getMaximumLikelihood
use mathsModule, only : invSigmoid
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
	
		scaleFactor = 0.2d0
		seed = 1234
		fbInt = 10
		maxStep = 10
		if (action(1:4) == 'CONT') then
			resume = 1
		else
			resume = 0
		endif
		nburn = 0
		flPfx = 'test'
		
		allocate(st(nVariables))
		allocate(stepSize(nVariables))
				
		if (resume == 0) then
			call initialValuesWeakField(st, stepSize)			
! 			call getMaximumLikelihood(st)
								
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
		scaleFactor = 0.5d0
		
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
	real(kind=8), allocatable :: Bpar(:), Bperp(:), azimuth(:), BModulus(:), fillFactor(:)

		loop = 1
		
		allocate(Bpar(npixels))
		allocate(Bperp(npixels))
		allocate(azimuth(npixels))
		allocate(BModulus(npixels))
		allocate(fillFactor(npixels))
		
! Compute maximum-likelihood solution
		Bpar = 0.5d0 * CV3 / CV2
		Bperp = sqrt(0.5d0 * sqrt(CQ3**2 + CU3**2) / CQ2)		
		azimuth = 0.5d0 * atan2(CU3, CQ3)
		fillFactor = 0.3d0
		BModulus = sqrt(Bpar**2 + Bperp**2) / fillFactor
		
		where (BModulus > BMax)
			BModulus = BMax - 200.d0
		endwhere
		where (BModulus < BMin)
			BModulus = BMin + 200.d0
		endwhere
		
		where (azimuth < 0.d0)
			azimuth = azimuth + PI
		endwhere
							
! B
		pars(1:nPixels) = BModulus
		stepSize = 100.d0		

! mu
		pars(nPixels+1:2*nPixels) = invSigmoid(cos(atan2(Bperp,Bpar)), muMin, muMax)
		stepSize = 2.d0
		
! f
		pars(2*nPixels+1:3*nPixels) = invSigmoid(fillFactor, fMin, fMax)
		stepSize = 2.d0
		
! phi
		pars(3*nPixels+1:4*nPixels) = invSigmoid(azimuth, phiMin, phiMax)
		stepSize = 2.d0
		
		loop = 4*nPixels + 1
! hyperB	
		pars(loop) = 1.d0
		stepSize(loop) = 2.2d0
		loop = loop + 1
		
		pars(loop) = 140.d0
		stepSize(loop) = 20.d0
		loop = loop + 1
				
! hypermu
		pars(loop) = 20.5d0
		stepSize(loop) = 10.d0
		loop = loop + 1
		
		pars(loop) = 20.5d0
		stepSize(loop) = 10.d0
		loop = loop + 1
		
! hyperf		
		pars(loop) = 8.d0
		stepSize(loop) = 5.d0
		loop = loop + 1
		
		pars(loop) = 20.5d0
		stepSize(loop) = 10.d0
		loop = loop + 1
		
		deallocate(Bpar)
		deallocate(Bperp)
		deallocate(azimuth)
		deallocate(BModulus)
		deallocate(fillFactor)
						
	end subroutine initialValuesWeakField
	
end module samplingModule