module ioModule
use mcmc_class_hierarchical, only : mcmc_hierarchical
implicit none
contains

!------------------------------------------------
! Read the file with the observations
!------------------------------------------------
	subroutine readObservations(chain)
	type(mcmc_hierarchical) :: chain
		
		open(unit=12,file='/scratch1/aasensio/HINODE/MAP_LITES/weakField.bin',status='old',action='read',form='unformatted')
		read(12) chain%nPixels
		write(*,*) 'Number of pixels ', chain%nPixels
		allocate(chain%CV2(chain%nPixels))
		allocate(chain%CV3(chain%nPixels))
		allocate(chain%CQ2(chain%nPixels))
		allocate(chain%CQ3(chain%nPixels))
		allocate(chain%CU2(chain%nPixels))
		allocate(chain%CU3(chain%nPixels))
		read(12) chain%CV2
		read(12) chain%CV3
		read(12) chain%CQ2
		read(12) chain%CQ3
		read(12) chain%CU2
		read(12) chain%CU3		
		close(12)
		
	end subroutine readObservations

end module ioModule