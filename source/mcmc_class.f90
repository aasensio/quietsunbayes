!*****************************************************************
! This module defines a class to do MCMC sampling using several algorithms:
!   - Adaptive Metropolis
!   - Metropolis-within-Gibbs
!
! Everything is general except for the function eval_posterior, which evaluates
! the log-posterior distribution of your problem for a given set of parameters
! Since you might need additional data to evaluate the posterior (priors, data, etc.)
! this function is usually overriden by a child class defined outwards
!*****************************************************************
module mcmc_class
	implicit none
! 	private
	real(kind=8) :: pi = 3.1415926535897931d0 ! Class-wide private constant

	type, public :: mcmc
		integer :: npar, step, step_accepted, alpha, steps_in_batch, algorithmID, iteration
		integer :: actual_chain_length
		real(kind=8), pointer :: pars(:), chain(:,:), mean(:), covariance(:,:), identity(:,:), ls_chain(:,:)
		real(kind=8), pointer :: trial(:), ls(:), lsaccepted(:), bestPars(:)
		real(kind=8) :: logp, logptrial, maxlogP
		character(len=40) :: algorithm
	contains
		procedure :: initChain
		procedure :: stepChain
		procedure :: stepMGChain
		procedure :: finalizeChain
		procedure :: getChain
		procedure :: getLSChain
		procedure :: thinChain
		procedure :: burninChain
		procedure :: updateStats
		procedure :: selectTrial
		procedure :: evalPosterior => logPosteriorExample
		procedure :: initialValues => initialValuesExample
	end type mcmc

	contains

!-------------------------------------------------------------
! Initialize the random number generator
!-------------------------------------------------------------
	subroutine init_random_seed()
	integer :: i, n, clock
   integer, dimension(:), allocatable :: seed
          
   	call random_seed(size = n)
      allocate(seed(n))
          
      call system_clock(count=clock)
          
      seed = clock + 37 * (/ (i - 1, i = 1, n) /)
      call random_seed(put = seed)
          
      deallocate(seed)
   end subroutine init_random_seed
          
!-------------------------------------------------------------
! Generates a random number following an uniform distribution in the interval [0,1]
! Call it with idum<0 for initialization
!-------------------------------------------------------------
	function randomu()
	real(kind=8) :: randomu
	
		call random_number(randomu)
							
	end function randomu
	
!-------------------------------------------------------------
! Generates a random number following an normal distribution with zero mean
! and unit variance
!-------------------------------------------------------------
	function randomn()

	real(kind=8) :: randomn

	real(kind=8) :: u, v

  		call random_number( u )
  		call random_number( v )

  		randomn = sqrt ( -2.0d0 * log ( u ) ) * cos ( 2.0d0 * pi * v )

		return
	end function randomn
!------------------------------------------------------------------
! Markov chain initialization
! npar: number of parameters of the chain
! algorithm: 'ADAPTIVE' / 'METROPOLIS_GIBBS'
! action: 'START' / 'CONT'
!------------------------------------------------------------------
	subroutine initChain(this, npar, algorithm, action)
	class(mcmc), intent(inout) :: this
	integer, intent(in) :: npar
	character(len=*), intent(in) :: algorithm, action
	integer :: i, outBounds

! Set number of parameters
		this%npar = npar

! Set sampling algorithm
		this%algorithm = algorithm

! Length of a chunk of chain
 		this%actual_chain_length = 10000

! Allocate memory for all variables
		allocate(this%pars(npar))
		allocate(this%bestPars(npar))
		allocate(this%chain(npar,this%actual_chain_length))
		allocate(this%mean(npar))
		allocate(this%covariance(npar,npar))
		allocate(this%trial(npar))
		allocate(this%identity(npar,npar))		

! Fill up the identity matrix
		this%identity = 0.d0
		do i = 1, npar
			this%identity(i,i) = 1.d0
		enddo

! Initialization of a few variables
		this%step = 0
		this%iteration = 0
		this%step_accepted = 0
		this%logp = 0.d0
		this%logptrial = 0
		
! Initialize random seed
		call init_random_seed()

! Adaptive Montecarlo
! algorithmID is set to avoid the comparison with the strings everytime we do a proposal
		if (index(this%algorithm,'METROPOLIS_GIBBS')) then
			this%algorithmID = 1
		endif

! If doing Metropolis-within-Gibbs sampling, allocate a few more variables
		if (index(this%algorithm,'METROPOLIS_GIBBS')) then
			allocate(this%ls(npar))
			allocate(this%lsaccepted(npar))
			allocate(this%ls_chain(npar,this%actual_chain_length))
			this%steps_in_batch = 0
			this%algorithmID = 2
			
! Read log-standard deviations from file
			if (index(action, 'CONT')) then
				print *, 'Reading continuation files...'
				open(unit=15,file='ls.out',action='read',status='old')
				read(15,*) (this%pars(i),i=1,npar)
				read(15,*) (this%ls(i),i=1,npar)
				read(15,*) (this%lsaccepted(i),i=1,npar)
				read(15,*) this%logp
				read(15,*) this%step, this%steps_in_batch
				close(15)
			endif
		endif

! First random sample inside the range of parameters and evaluate posterior
		if (index(action, 'START')) then
			call this%initialValues()
		endif

		
		this%trial = this%pars		
		this%logp = this%evalPosterior(outBounds)

		this%maxlogP = -1.d100
		
	end subroutine initChain

!------------------------------------------------------------------
! Set initial values for the parameters
!------------------------------------------------------------------
	subroutine initialValuesExample(this)
	class(mcmc), intent(inout) :: this

		this%pars = 1.d0

	end subroutine initialValuesExample

!------------------------------------------------------------------
! Return the chain
!------------------------------------------------------------------
	function getChain(this) result(out)
	class(mcmc), intent(in) :: this
	real(kind=8) :: out(this%iteration,this%npar)
	
		out = this%chain(:,1:this%iteration)

	end function getChain

	!------------------------------------------------------------------
! Return the chain
!------------------------------------------------------------------
	function getLSChain(this) result(out)
	class(mcmc), intent(in) :: this
	real(kind=8) :: out(this%iteration,this%npar)
	
		out = this%ls_chain(:,1:this%iteration)

	end function getLSChain

!------------------------------------------------------------------
! Markov chain step
!------------------------------------------------------------------
	subroutine stepChain(this)
	class(mcmc), intent(in) :: this
	real(kind=8) :: area
		area = pi * this%npar
	end subroutine stepChain

!------------------------------------------------------------------
! Markov chain Metropolis-within-Gibbs step
!------------------------------------------------------------------
	subroutine stepMGChain(this)
	class(mcmc), intent(inout) :: this
	real(kind=8) :: logp, old_pars(this%npar), trial(this%npar), acceptance, n_batch, delta_n
	real(kind=8) :: logptrial, alpha, r, ran
	real(kind=8), allocatable :: temp(:,:)
	logical :: reset_batch
	integer :: i, outBounds, filepos
		
! Add one step to the chain counters
		this%step = this%step + 1
		this%iteration = this%iteration + 1
		this%steps_in_batch = this%steps_in_batch + 1
		reset_batch = .FALSE.
		
		logp = this%logp

! Save the original set of parameters	
		this%trial = this%pars
		old_pars = this%pars

! Propose a new set of parameters
		trial = this%selectTrial()

! Do the Gibbs updating
		do i = 1, this%npar

! Set the i-th parameter to its new value and do a Metropolis step
			this%trial(i) = trial(i)

! Evaluate the posterior
			logptrial = this%evalPosterior(outBounds)
			
! If it is inside the boundaries
			if (outBounds == 0) then

! Metropolis-Hastings step verification
				r = exp(logptrial - logp)
				alpha = min(1.d0,r)
				ran = randomu()

! If the MH step is successful, we keep the modification of the i-th parameter
				if (ran < alpha) then
					this%pars(i) = this%trial(i)
					this%lsaccepted(i) = this%lsaccepted(i) + 1
					logp = logptrial
				else
! If not, return the original value of the parameter
					this%trial(i) = old_pars(i)
				endif
			else
! If not, return the original value of the parameter
				this%trial(i) = old_pars(i)
			endif

! Modify the logarithm of the variances to adapt the acceptance rate
			if ( modulo(this%step+1, 50) == 0) then
				acceptance = 1.d0 * this%lsaccepted(i) / this%steps_in_batch
				n_batch = this%step / 50.d0

				delta_n = min(0.01d0, 1.d0 / sqrt(n_batch))				
				
				if (acceptance >= 0.44d0) then
					this%ls(i) = this%ls(i) + delta_n
				endif
				if (acceptance < 0.44d0) then
					this%ls(i) = this%ls(i) - delta_n
				endif
				
				reset_batch = .TRUE.
			endif

		enddo

		this%logp = logp
		
! Update the chain in the output file
		if (this%iteration > this%actual_chain_length) then
			
			inquire(unit=12, size=filepos)
			write(12,pos=filepos+1) this%chain
			
			inquire(unit=13, size=filepos)
			write(13,pos=filepos+1) this%ls_chain

			open(unit=14,file='posterior.sizes',action='write',status='replace')
			write(14,*) this%npar, this%step
			close(14)
											
! Reset the iteration number
			this%iteration = 1
		
		endif
		
		this%chain(:,this%iteration) = this%pars
		this%ls_chain(:,this%iteration) = this%ls
		
! The next iteration starts a new batch of 50 steps
		if (reset_batch) then
			this%lsaccepted = 0
			this%steps_in_batch = 0
		endif

! Set the maximum value of the log-posterior
		if (this%logP > this%maxlogP) then
			this%maxlogP = this%logP
			this%bestPars = this%pars
		endif
		
	end subroutine stepMGChain

!------------------------------------------------------------------
! Finalize the chain writing all samples that haven't been written
!------------------------------------------------------------------
	subroutine finalizeChain(this)
	class(mcmc), intent(inout) :: this
	integer :: filepos
	
		inquire(unit=12, size=filepos)
		write(12,pos=filepos+1) this%chain(:,1:this%iteration)
			
		inquire(unit=13, size=filepos)
		write(13,pos=filepos+1) this%ls_chain(:,1:this%iteration)

		open(unit=14,file='posterior.sizes',action='write',status='replace')
		write(14,*) this%npar, this%step
		close(14)
					
	end subroutine finalizeChain

!------------------------------------------------------------------
! Markov chain thinning by a factor
!------------------------------------------------------------------
	subroutine thinChain(this, factor)
	class(mcmc), intent(inout) :: this
	integer :: factor
	real(kind=8), allocatable :: temp(:,:)
	integer :: i, loop
	
		allocate(temp(this%step / factor,this%npar))
		loop = 1
		do i = 1, this%step, factor
			temp(loop,:) = this%chain(:,i)
		enddo
		
! Now reallocate memory for the chain
		deallocate(this%chain)
		allocate(this%chain(this%npar,loop))
		this%chain = temp(:,1:loop)
		deallocate(temp)
		
! Change the length of the chain
		this%step = loop		
		
	end subroutine thinChain

!------------------------------------------------------------------
! Markov chain burn-in
!------------------------------------------------------------------
	subroutine burninChain(this, initial_steps)
	class(mcmc), intent(inout) :: this
	integer :: initial_steps
	real(kind=8), allocatable :: temp(:,:)
	integer :: i, loop
	
		allocate(temp(this%step - initial_steps,this%npar))
		temp = this%chain(:,initial_steps:this%step)
		
! Now reallocate memory for the chain
		deallocate(this%chain)
		allocate(this%chain(this%npar,this%step - initial_steps))	
		this%chain = temp
		deallocate(temp)
		
! Change the length of the chain
		this%step = this%step - initial_steps
	end subroutine burninChain

!------------------------------------------------------------------
! Update mean and covariance matrix
!------------------------------------------------------------------
	subroutine updateStats(this)
	class(mcmc), intent(in) :: this
	real(kind=8) :: mean_old(this%npar)
	integer :: i, j

! Update the mean
		mean_old = this%mean
		this%mean = mean_old + (this%pars - mean_old) / (this%step+1.d0)
 
! Update the covariance matrix
		do i = 1, this%npar
			do j = 1, this%npar
				this%covariance(i,j) = (this%step-1.d0)/this%step * this%covariance(i,j) + &
                                         (this%pars(i)-mean_old(i))*(this%pars(j)-mean_old(j)) / &
                                         (this%step+1.d0)**2 + (this%pars(i)-this%mean(i))*&
                                         (this%pars(j)-this%mean(j)) / (this%step*1.d0)
			enddo
		enddo
	end subroutine updateStats

!------------------------------------------------------------------
! Select a new trial point following a certain proposal density
!------------------------------------------------------------------
	function selectTrial(this) result (out)
	class(mcmc), intent(in) :: this
	real(kind=8) :: out(this%npar)
	integer :: i

! Adaptive Metropolis-within-Gibbs
! Draw npar random numbers from N(0,sigma^2), with sigma appropriate for each parameter
		if (this%algorithmID == 2) then
			do i = 1, this%npar
				out(i) = this%pars(i) + exp(this%ls(i)) * randomn()
			enddo
		endif
		
		return
		
	end function selectTrial

!------------------------------------------------------------------
! Evaluate the posterior
!------------------------------------------------------------------
	function logPosteriorExample(this, outBounds) result(logP)
	class(mcmc), intent(in) :: this
	integer, intent(inout) :: outBounds
	real(kind=8) :: logP
		logP = 0.d0
	end function logPosteriorExample

end module mcmc_class