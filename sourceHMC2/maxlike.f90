module maximumlikelihoodModule
use varsModule
use likelihoodModule, only : negLogPosterior
use lbfgsbModule, only : setulb
implicit none
	
contains

!-----------------------------------------------------------------------
! Do the actual nested sampling to estimate the evidence and the posterior distribution
!-----------------------------------------------------------------------
	subroutine getMaximumLikelihood(x)
	real(kind=8) :: x(:)
   integer :: i, loop
   integer :: nbin, mmax, ndata, iprint, k
	character(len=60) :: task, csave
	logical :: lsave(4)
	integer, allocatable :: nbd(:), iwa(:), isave(:)
	real(kind=8) :: factr, pgtol, f, f2
	real(kind=8), allocatable :: l(:), u(:), g(:), dsave(:), wa(:), x2(:), g2(:)

		allocate(l(nVariables),u(nVariables),g(nVariables))

! Boundaries
		loop = 1
		l(1:nPixels) = BMin
		u(1:nPixels) = BMax
		
		l(nPixels+1:2*nPixels) = muMin + 1.d-3
		u(nPixels+1:2*nPixels) = muMax - 1.d-3
		
		l(2*nPixels+1:3*nPixels) = fMin - 1.d-3
		u(2*nPixels+1:3*nPixels) = fMax + 1.d-3
		
		l(3*nPixels+1:4*nPixels) = phiMin - 1.d-3
		u(3*nPixels+1:4*nPixels) = phiMax + 1.d-3
								
! Geodesic Levenberg-Marquardt
		mmax = 10

		allocate(nbd(nVariables),iwa(3*nVariables),isave(44))
		allocate(dsave(29),wa(2*mmax*nVariables+5*nVariables+11*mmax*mmax+8*mmax),x2(nVariables),g2(nVariables))
		
		nbd = 2
		
		iprint = 1
		factr = 1.d7
		pgtol = 1.d-8
		
		task = 'START'

		do while (task(1:2) == 'FG' .or. task(1:5) == 'NEW_X' .or. task(1:5) == 'START')			
			call setulb(nVariables,mmax,x,l,u,nbd,f,g,factr,pgtol,wa,iwa,task,iprint, csave,lsave,isave,dsave)			

! Compute function and gradient
			if (task(1:2) == 'FG') then				
				call negLogPosterior(nVariables,x,f,g)
			endif

! New iterate and continue the iteration
			if (task(1:5) == 'NEW_X') then
				open(unit=15,file='bestPars',status='replace',action='write')
				write(15,*) x
				close(15)
			endif
		enddo
		
		deallocate(nbd,iwa,isave)
		deallocate(g,dsave,wa,x2,g2)
												
	end subroutine getMaximumLikelihood
	
end module maximumlikelihoodModule