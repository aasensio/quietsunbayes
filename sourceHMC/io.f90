module ioModule
use varsModule, only : CV1, CV2, CV3, CQ1, CQ2, CQ3, CU1, CU2, CU3, nPixels
implicit none
contains

!------------------------------------------------
! Read the file with the observations
!------------------------------------------------
	subroutine readObservations
			
		open(unit=12,file='/home/aasensio/Dropbox/weakField.bin',status='old',action='read',form='unformatted')
		read(12) nPixels
		write(*,*) 'Number of pixels ', nPixels
		allocate(CV1(nPixels))
		allocate(CV2(nPixels))
		allocate(CV3(nPixels))
		allocate(CQ1(nPixels))
		allocate(CQ2(nPixels))
		allocate(CQ3(nPixels))
		allocate(CU1(nPixels))
		allocate(CU2(nPixels))
		allocate(CU3(nPixels))
		
		read(12) CV1
		read(12) CV2
		read(12) CV3
		
		read(12) CQ1
		read(12) CQ2
		read(12) CQ3
		
		read(12) CU1
		read(12) CU2
		read(12) CU3		
		close(12)
		
	end subroutine readObservations

end module ioModule