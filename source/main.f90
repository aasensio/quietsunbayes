program main
use ioModule, only : readObservations
! use samplingModule, only : doSampling
use mcmc_class_hierarchical
implicit none

	type(mcmc_hierarchical) :: chain
	integer :: i, nargs, nsteps, npar
	character(len=100) :: action, arg
	
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
		read(arg,*) nsteps
	endif

! Read the observations
	print *, 'Reading observations'
	call readObservations(chain)
	
! Read ranges of parameters
	open(unit=12,file='conf.dat',action='read',status='old')
	read(12,*)
	read(12,*) chain%BMin, chain%BMax, chain%BNodes, chain%BSigma2
	read(12,*) chain%muMin, chain%muMax, chain%muNodes, chain%muSigma2
	read(12,*) chain%fMin, chain%fMax, chain%fNodes, chain%fSigma2
	read(12,*) chain%phiMin, chain%phiMax, chain%phiNodes	, chain%phiSigma2
	close(12)
	
	chain%BSigma2 = chain%BSigma2**2
	chain%muSigma2 = chain%muSigma2**2
	chain%fSigma2 = chain%fSigma2**2
	chain%phiSigma2 = chain%phiSigma2**2
	
	allocate(chain%BLoc(chain%BNodes))
	allocate(chain%muLoc(chain%muNodes))
	allocate(chain%fLoc(chain%fNodes))
	allocate(chain%phiLoc(chain%phiNodes))
	do i = 1, chain%BNodes
		chain%BLoc(i) = (i-1.d0) * (chain%BMax - chain%BMin) / (chain%BNodes-1.d0) + chain%BMin
	enddo
	print *, 'B hyperparameter positions : ', chain%BLoc
	do i = 1, chain%BNodes
		chain%muLoc(i) = (i-1.d0) * (chain%muMax - chain%muMin) / (chain%muNodes-1.d0) + chain%muMin
	enddo
	print *, 'mu hyperparameter positions : ', chain%muLoc
	do i = 1, chain%BNodes
		chain%fLoc(i) = (i-1.d0) * (chain%fMax - chain%fMin) / (chain%fNodes-1.d0) + chain%fMin
	enddo
	print *, 'f hyperparameter positions : ', chain%fLoc
	do i = 1, chain%BNodes
		chain%phiLoc(i) = (i-1.d0) * (chain%phiMax - chain%phiMin) / (chain%phiNodes-1.d0) + chain%phiMin
	enddo	
	print *, 'phiB hyperparameter positions : ', chain%phiLoc
		
	print *, 'Data read'

! Allocate memory for the model parameters
	allocate(B_i(chain%npixels))
	allocate(mu_i(chain%npixels))
	allocate(f_i(chain%npixels))
	allocate(phi_i(chain%npixels))
	
! Hyperparameters
	allocate(thetaB_i(chain%BNodes))
	allocate(thetamu_i(chain%muNodes))
	allocate(thetaf_i(chain%fNodes))
	allocate(thetaphi_i(chain%phiNodes))
	
! Number of parameters
! - B, mu, f and phi for each pixel plus the hyperparameters
	npar = chain%npixels * 4 + (chain%BNodes + chain%muNodes + chain%fNodes + chain%phiNodes)

	print *, 'Number of parameters = ', npar

! If this is the start of a chain, open the file destroying the previous one
	if (index(action,'START')) then
		open(unit=12,file='posterior',action='write',status='replace',form='unformatted',access='stream')
		open(unit=13,file='proposal',action='write',status='replace',form='unformatted',access='stream')
	endif
	
! If it is a continuation, open the file without destroying it
	if (index(action,'CONT')) then
		open(unit=12,file='posterior',action='write',status='old',form='unformatted', access='stream')		 		
		open(unit=13,file='proposal',action='write',status='old',form='unformatted', access='stream')		
	endif

	call chain%initChain(npar, 'METROPOLIS_GIBBS', action)

	do i = 1, nsteps
		if (modulo(i, 1000) == 0) then
			write(*,FMT='(A,I7,A,I7,A,F5.1,A,E)') 'Iteration : ', i,'/',nsteps,' (', 100.d0*i/nsteps,'%) - logP(max) = ', chain%maxlogP
		endif
		call chain%stepMGChain()
	enddo

	print *, 'Finished.'

! Finalize chain
	call chain%finalizeChain()
	
	close(12)
	close(13)
		
		
	open(unit=12,file='bestpars',action='write',status='replace')
	write(12,*) chain%bestpars
	close(12)
	
	open(unit=12,file='ls.out',action='write',status='replace')
	write(12,*) chain%pars
	write(12,*) chain%ls
	write(12,*) chain%lsaccepted
	write(12,*) chain%logp
	write(12,*) chain%step, chain%steps_in_batch
	close(12)
	
end program main