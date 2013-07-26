module samplingModule
use varsModule, only : observationType
implicit none
contains

!------------------------------------------------
! Negative log-posterior and its gradient
!------------------------------------------------
! 	subroutine doSampling(obs)
! 	type(observationType) :: obs
! 	real(kind=8), allocatable :: st(:), stepSize(:)
! 	real(kind=8) :: scaleFactor
! 	integer :: seed, fbInt, maxStep, resume
! 	character(len=128) :: flPfx
! 	
! 		scaleFactor = 1.d0
!  		seed = 1234
! 		fbInt = 10
! 		maxStep = 10
! 		resume = 0
! 		flPfx = 'test'
! 		
! 		allocate(st(sdim))
! 		allocate(stepSize(sdim))
! 		
! 		
! 	end subroutine doSampling
! 	
! !------------------------------------------------
! ! Negative log-posterior and its gradient
! !------------------------------------------------
! 	subroutine negLogPosterior(ndim,x,v,g)
! 	integer ndim
! 	real(kind=8), dimension(ndim) :: x
!    real(kind=8), dimension(ndim) :: g
!    real(kind=8) :: v
! 	integer i				
! 		
! 		v = -v		
! 		g = -g
!   
! 	end subroutine negLogPosterior
! 
! !------------------------------------------------
! ! A subroutie to write the extract file
! ! I have assumed that the unit=20 is opened for
! ! writing (append) earlier. In general only write
! ! those parametes which are estimated. The files
! ! can be really big depending on the dimensionality
! !------------------------------------------------
! 	subroutine writeHMCProcess(ndim,x,v,g)
! 	integer ndim
! 	real(kind=8), dimension(ndim) :: x, x2
!    real(kind=8), dimension(ndim) :: g
!    real(kind=8) :: v
! 	integer i	
! 	
! 	do i=1,ndim-1				
! 		write(20,'(E18.10,a)',advance='no') x(i), ' '		
! 	enddo
! 		
! 	write(20,'(E18.10)') x(ndim)
! 		
! 	end subroutine writeHMCProcess

end module samplingModule