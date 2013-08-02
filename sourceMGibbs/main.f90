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
	read(12,*) chain%BMin, chain%BMax
	read(12,*) chain%muMin, chain%muMax
	read(12,*) chain%fMin, chain%fMax
	read(12,*) chain%phiMin, chain%phiMax
	close(12)
	
	print *, 'B range : ', chain%BMin, chain%BMax
 	print *, 'mu range : ', chain%muMin, chain%muMax
 	print *, 'f range : ', chain%fMin, chain%fMax
 	print *, 'phi range : ', chain%phiMin, chain%phiMax
				
	print *, 'Data read'

! Allocate memory for the model parameters
	allocate(B_i(chain%npixels))
	allocate(mu_i(chain%npixels))
	allocate(f_i(chain%npixels))
	allocate(phi_i(chain%npixels))
	
	allocate(sqrtMu(chain%npixels))
	allocate(c2p(chain%npixels))
	allocate(s2p(chain%npixels))
	
! Hyperparameters
	allocate(hyperB_i(2))
	allocate(hypermu_i(2))
	allocate(hyperf_i(2))
	allocate(hyperphi_i(2))
	
! Number of parameters
! - B, mu, f and phi for each pixel plus the hyperparameters
	npar = chain%npixels * 4 + 8
	
	allocate(chain%parsToSave(8))
	do i = 1, 8
		chain%parsToSave(i) = chain%npixels*4+i
	enddo

	print *, 'Number of parameters = ', npar
	print *, 'Saving only hyperparameters '

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