program main
use ioModule, only : readObservations
use samplingModule, only : doSampling
use varsModule
use mathsModule, only : digama
implicit none
	
	integer :: i, j, nargs
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
	
	allocate(hyperparRanges(2,2))
	
! Read ranges of parameters
	open(unit=12,file='conf_fmarginalized.dat',action='read',status='old')
	read(12,*)
	read(12,*) fileObs
	read(12,*)
	read(12,*)
	read(12,*) BMin, BMax
	read(12,*) muMin, muMax	
	read(12,*) phiMin, phiMax
	read(12,*)
	read(12,*)
	do i = 1, 2
		read(12,*) (hyperparRanges(j,i),j=1,2)
	enddo
	close(12)
	
	print *, 'B range : ', BMin, BMax
 	print *, 'mu range : ', muMin, muMax 	
 	print *, 'phi range : ', phiMin, phiMax
				
! Read the observations
	print *, 'Reading observations'
	call readObservations(fileObs)
	
	print *, 'Data read'
	
	sigma_n = 1.d-3


! Allocate memory for the model parameters
	allocate(B_i(npixels))
	allocate(mu_i(npixels))	
	allocate(phi_i(npixels))
	
	allocate(sqrtMu(npixels))
	allocate(dsqrtMudMu(npixels))
	allocate(c2p(npixels))
	allocate(s2p(npixels))
	
	allocate(A(npixels))
	allocate(B(npixels))
	allocate(C(npixels))	
	allocate(dlogpdB(npixels))
	allocate(dlogpdC(npixels))
	allocate(dBdB_i(npixels), dBdmu_i(npixels), dBdphi_i(npixels), dCdB_i(npixels), dCdmu_i(npixels), dCdphi_i(npixels))
	
! Hyperparameters
	allocate(hyperB_i(2))
	allocate(hypermu_i(2))
	
! Number of parameters
! - B, mu, f and phi for each pixel plus the hyperparameters
	nVariables = nPixels * 3 + 4
	
	allocate(parsOld(nVariables))
	allocate(parsInitial(nVariables))
	allocate(parsMean(nVariables))
	allocate(parsVariance(nVariables))
	
	allocate(parsToSave(4))
	do i = 1, 4
		parsToSave(i) = nPixels*3+i
	enddo

	print *, 'Number of parameters = ', nVariables
	print *, 'Saving only hyperparameters '
 	
	call doSampling(action)
	
	print *, 'Finished.'
	
	close(12)
	close(13)
				
end program main