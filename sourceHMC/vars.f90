module varsModule
implicit none

	integer :: nPixels, nVariables, nSteps
		
	real(kind=8) :: BMin, BMax, muMin, muMax, fMin, fMax, phiMin, phiMax
		
	real(kind=8), allocatable, dimension(:) :: CV1, CV2, CV3, CQ1, CQ2, CQ3, CU1, CU2, CU3
		
	real(kind=8), allocatable, dimension(:) :: B_i, mu_i, phi_i, f_i, hyperB_i, hypermu_i, hyperphi_i, hyperf_i, parsToSave
	
	real(kind=8), allocatable, dimension(:) :: c2p, s2p, sqrtMu, dsqrtMudMu

end module varsModule