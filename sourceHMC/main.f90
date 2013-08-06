program main
use ioModule, only : readObservations
use samplingModule, only : doSampling
use varsModule
use mathsModule, only : digama
implicit none
	
	integer :: i, nargs
	character(len=100) :: action, arg
	character(len=128) :: fileObs
	
! Find the parameters passed to the program
	nargs = iargc()
	
	if (nargs <= 1) then
		print *, 'Parameters missing. The options are:'
		print *, '  - START nsteps -> start a MCMC with nsteps'
		print *, '  - CONT nsteps -> continue a previous MCMC and proceed another nsteps steps'		
		stop
	endif
	
	if (nargs == 2) then
		call getarg(1,action)
		call getarg(2,arg)
		read(arg,*) nSteps
	endif
	
! Read ranges of parameters
	open(unit=12,file='conf.dat',action='read',status='old')
	read(12,*)
	read(12,*) fileObs
	read(12,*)
	read(12,*)
	read(12,*) BMin, BMax
	read(12,*) muMin, muMax
	read(12,*) fMin, fMax
	read(12,*) phiMin, phiMax
	close(12)
	
	print *, 'B range : ', BMin, BMax
 	print *, 'mu range : ', muMin, muMax
 	print *, 'f range : ', fMin, fMax
 	print *, 'phi range : ', phiMin, phiMax
				
! Read the observations
	print *, 'Reading observations'
	call readObservations(fileObs)
	
	print *, 'Data read'
	
	sigma_n = 1.d-3


! Allocate memory for the model parameters
	allocate(B_i(npixels))
	allocate(mu_i(npixels))
	allocate(f_i(npixels))
	allocate(phi_i(npixels))
	
	allocate(sqrtMu(npixels))
	allocate(dsqrtMudMu(npixels))
	allocate(c2p(npixels))
	allocate(s2p(npixels))
	
! Hyperparameters
	allocate(hyperB_i(2))
	allocate(hypermu_i(2))
	allocate(hyperf_i(2))
	allocate(hyperphi_i(2))
	
! Number of parameters
! - B, mu, f and phi for each pixel plus the hyperparameters
	nVariables = nPixels * 4 + 8
	
	allocate(parsOld(nVariables))
	allocate(parsMean(nVariables))
	allocate(parsVariance(nVariables))
	
	allocate(parsToSave(8))
	do i = 1, 8
		parsToSave(i) = nPixels*4+i
	enddo

	print *, 'Number of parameters = ', nVariables
	print *, 'Saving only hyperparameters '

	call doSampling
	
	print *, 'Finished.'
	
	close(12)
	close(13)
				
end program main