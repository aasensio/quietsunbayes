module likelihoodModule
use varsModule
use mathsModule, only : alngam, digama
implicit none

contains
!------------------------------------------------
! Negative log-posterior and its gradient
!------------------------------------------------
 	subroutine negLogPosterior(nVariables,trial,logP,logPGradient)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: trial
   real(kind=8), dimension(nVariables) :: logPGradient
   real(kind=8) :: logP, nu, logP2
	integer :: i, outBounds, loop, ierr
		
! Fill the vectors to simplify the notation and test for the boundaries		
		loop = 1
		
! B
		do i = 1, nPixels
			B_i(i) = trial(loop)
			loop = loop + 1
			if (B_i(i) < BMin .or. B_i(i) > BMax) then
				logP = 1.d15
				logPGradient = 1.d15				
				return
			endif
		enddo

! mu
		do i = 1, npixels
			mu_i(i) = trial(loop)
			loop = loop + 1
			if (mu_i(i) < muMin .or. mu_i(i) > muMax) then
				logP = 1.d15
				logPGradient = 1.d15				
				return
			endif
		enddo
	
! f
		do i = 1, npixels
			f_i(i) = trial(loop)
			loop = loop + 1
			if (f_i(i) < fMin .or. f_i(i) > fMax) then
				logP = 1.d15
				logPGradient = 1.d15				
				return
			endif
		enddo
	
! phi
		do i = 1, npixels
			phi_i(i) = trial(loop)
			loop = loop + 1
			if (phi_i(i) < phiMin .or. phi_i(i) > phiMax) then
				logP = 1.d15
				logPGradient = 1.d15				
				return
			endif
		enddo

				
 		logP = 0.d0
		logPGradient = 0.d0
 		nu = 0.001d0
 				
!-----------------
! DATA LOG-LIKELIHOOD
!-----------------
		c2p = cos(2.d0 * phi_i)
		s2p = sin(2.d0 * phi_i)		
		
  		logP = logP - 0.5d0 / sigma_n**2 * (&
			sum( CV1 + (B_i * mu_i * f_i)**2 * CV2 - (B_i * mu_i * f_i) * CV3 ) + &
			sum( CQ1 + (B_i**2 * (1.d0-mu_i**2) * f_i * c2p)**2 * CQ2 - (B_i**2 * (1.d0-mu_i**2) * f_i * c2p) * CQ3 ) + &
			sum( CU1 + (B_i**2 * (1.d0-mu_i**2) * f_i * s2p)**2 * CU2 - (B_i**2 * (1.d0-mu_i**2) * f_i * s2p) * CU3 ) )
			
! dlogL/dB
		logPGradient(1:nPixels) = logPGradient(1:nPixels) - 0.5d0 / sigma_n**2 * (&
			( 2.d0 * B_i * (mu_i * f_i)**2 * CV2 - (mu_i * f_i) * CV3 ) + &
			( 4.d0 * B_i**3 * ((1.d0-mu_i**2) * f_i * c2p)**2 * CQ2 - 2.d0 * B_i * ((1.d0-mu_i**2) * f_i * c2p) * CQ3 ) + &
			( 4.d0 * B_i**3 * ((1.d0-mu_i**2) * f_i * s2p)**2 * CU2 - 2.d0 * B_i * ((1.d0-mu_i**2) * f_i * s2p) * CU3 ) )
			
! dlogL/dmu
		logPGradient(nPixels+1:2*nPixels) = logPGradient(nPixels+1:2*nPixels) - 0.5d0 / sigma_n**2 * (&
			( 2.d0 * mu_i * (B_i * f_i)**2 * CV2 - (B_i * f_i) * CV3 ) + &
			( 2.d0 * (1.d0-mu_i**2) * (-2.d0*mu_i) * (B_i**2 * f_i * c2p)**2 * CQ2 - (-2.d0*mu_i) * (B_i**2 * f_i * c2p) * CQ3 ) + &
			( 2.d0 * (1.d0-mu_i**2) * (-2.d0*mu_i) * (B_i**2 * f_i * s2p)**2 * CU2 - (-2.d0*mu_i) * (B_i**2 * f_i * s2p) * CU3 ) )
			
! dlogL/df
		logPGradient(2*nPixels+1:3*nPixels) = logPGradient(2*nPixels+1:3*nPixels) - 0.5d0 / sigma_n**2 * (&
			( 2.d0 * f_i * (B_i * mu_i)**2 * CV2 - (B_i * mu_i) * CV3 ) + &
			( 2.d0 * f_i * (B_i**2 * (1.d0-mu_i**2) * c2p)**2 * CQ2 - (B_i**2 * (1.d0-mu_i**2) * c2p) * CQ3 ) + &
			( 2.d0 * f_i * (B_i**2 * (1.d0-mu_i**2) * s2p)**2 * CU2 - (B_i**2 * (1.d0-mu_i**2) * s2p) * CU3 ) )
			
! dlogL/dphi
		logPGradient(3*nPixels+1:4*nPixels) = logPGradient(3*nPixels+1:4*nPixels) - 0.5d0 / sigma_n**2 * (&			
			( -4.d0 * c2p * s2p * (B_i**2 * (1.d0-mu_i**2) * f_i)**2 * CQ2 + 2.d0 * s2p * (B_i**2 * (1.d0-mu_i**2) * f_i) * CQ3 ) + &
			( 4.d0 * s2p * c2p * (B_i**2 * (1.d0-mu_i**2) * f_i)**2 * CU2 - 2.d0 * c2p * (B_i**2 * (1.d0-mu_i**2) * f_i) * CU3 ) )
			
					
		logP = -logP
		logPGradient = -logPGradient
						  
	end subroutine negLogPosterior

!------------------------------------------------
! A subroutine to write the extract file
! I have assumed that the unit=20 is opened for
! writing (append) earlier. In general only write
! those parametes which are estimated. The files
! can be really big depending on the dimensionality
!------------------------------------------------
	subroutine writeHMCProcess(nVariables,x,v,g)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: x
   real(kind=8), dimension(nVariables) :: g, meanOld
   real(kind=8) :: v, xwrite(6)
	integer i
			
		meanOld = parsMean
		parsMean = meanOld + (x - meanOld) / (nStep + 1.d0)		
		parsVariance = (nStep - 1.d0) / nStep * parsVariance + (x - meanOld)**2 / (nStep+1.d0)**2 + (x - meanOld)**2 / nStep
		
		xwrite = x(4*nPixels+1:nVariables)
		write(20) xwrite
		
		write(21) x
		
		open(unit=25,file= "test.stddev",action='write',status='replace')
		write(25,*) sqrt(parsVariance)
		close(25)
		
		open(unit=25,file= "test.mean",action='write',status='replace')
		write(25,*) parsMean
		close(25)
		
		nStep = nStep + 1.d0
				
	end subroutine writeHMCProcess
	
!------------------------------------------------
! A subroutine to write the extract file
! I have assumed that the unit=20 is opened for
! writing (append) earlier. In general only write
! those parametes which are estimated. The files
! can be really big depending on the dimensionality
!------------------------------------------------
	subroutine writeHMCProcessBurnin(nVariables,x,v,g)
	integer :: nVariables
	real(kind=8), dimension(nVariables) :: x
   real(kind=8), dimension(nVariables) :: g, meanOld
   real(kind=8) :: v, xwrite(8)
	integer i
							
		write(20) x		
						
	end subroutine writeHMCProcessBurnin
		
end module likelihoodModule