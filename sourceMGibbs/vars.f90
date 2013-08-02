module varsModule
implicit none

	type observationType
		integer :: nPixels
		real(kind=8), pointer, dimension(:) :: CV2, CV3, CQ2, CQ3, CU2, CU3
	end type observationType
	
	type(observationType) :: obs

end module varsModule