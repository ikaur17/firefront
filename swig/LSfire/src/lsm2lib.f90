MODULE LSM2
    
    ! Compilazione:
    ! f2py -c -m libhypef NT_parameters.f90 NT_lib_mixture3_r1.f90 NT_lib_mixture3_r2.f90 NT_v06.f90 --f90flags="-cpp -fopenmp " -lgomp

    !USE NT_parameters
    USE iso_c_binding
    
#ifdef _OPENMP
    USE omp_lib
#endif

    implicit none

#define FLAG_FORALL

    
    integer, parameter :: dp = kind(1.0d0)
    
    character(len=3), parameter :: VERSION = "0.1"
    
    logical, parameter :: flag_debug = .false.
    
    
CONTAINS
 
	!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

	SUBROUTINE lsmtest(N, c, i1, i2, Nthreads) BIND(C, name = "lsm2_test")
            
		integer, intent(in) :: N, Nthreads
		integer, intent(inout) :: c, i1, i2 
		integer :: i
        				
#ifdef _OPENMP
		call omp_set_num_threads(Nthreads)
		i1 = 1
		i2 = -1
		!$OMP PARALLEL DO
		do i = 1, N
			c = c + omp_get_thread_num() + 1
		end do
		!$OMP END PARALLEL DO
#else
		i1 = -1
		i2 = 1
		do i = 1, N
			c = c + 1
		end do
		c = -c
#endif

#ifdef FLAG_FORALL
		c = 0
#endif

	END SUBROUTINE lsmtest
		
		
!***********************************************************************
!
!  File:        lsm_level_set_evolution2d.f
!  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
!                   Regents of the University of Texas.  All rights reserved.
!               (c) 2009 Kevin T. Chu.  All rights reserved.
!  Revision:    $Revision$
!  Modified:    $Date$
!  Description: F77 subroutines for 2D level set evolution equation
!
!***********************************************************************
! toolbox/lsm_level_set_evolution/lsm_level_set_evolution2d.f
!***********************************************************************
	SUBROUTINE lsm2dZeroOutLevelSetEqnRHS( &
		lse_rhs, &
		ilo_lse_rhs_gb, ihi_lse_rhs_gb, &
		jlo_lse_rhs_gb, jhi_lse_rhs_gb, Nthreads) &
		BIND(C, name = "LSM2_LSM2D_ZERO_OUT_LEVEL_SET_EQN_RHS")
		
		integer :: ilo_lse_rhs_gb, ihi_lse_rhs_gb
		integer :: jlo_lse_rhs_gb, jhi_lse_rhs_gb
		real(dp) :: lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb, &
			jlo_lse_rhs_gb:jhi_lse_rhs_gb)
		integer :: i, j, Nthreads
		
#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif
		
		
#ifdef FLAG_FORALL
		forall (j=jlo_lse_rhs_gb:jhi_lse_rhs_gb, i=ilo_lse_rhs_gb:ihi_lse_rhs_gb)
			lse_rhs(i,j) = 0.d0
		end forall
#else
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_lse_rhs_gb, jhi_lse_rhs_gb
			forall (i = ilo_lse_rhs_gb : ihi_lse_rhs_gb)
				lse_rhs(i,j) = 0.d0
			end forall
		end do
#endif 
		
	END SUBROUTINE lsm2dZeroOutLevelSetEqnRHS
	

!***********************************************************************
!
!  lsm2dComputeDn() computes the n-th undivided differences in the
!  specified direction given the (n-1)-th undivided differences.  The 
!  undivided differences in cells with insufficient data is set to a 
!  large number.
!
!  Arguments:
!    Dn (out):           n-th undivided differences 
!    Dn_minus_one (in):  (n-1)-th undivided differences 
!    n (in):             order of undivided differences to compute
!    *_gb (in):          index range for ghostbox
!    *_fb (in):          index range for fillbox
!
!  NOTES:
!   - The index ranges for all ghostboxes and the fillbox should 
!     correspond to the index range for cell-centered data.
!   - The undivided differences for odd n are face-centered (i.e.
!     indices are of the form (i+1/2)).  In this situation, the array
!     index corresponding to the (i+1/2)-th undivided difference is
!     i (i.e. the index shifted down to the nearest integer index). 
!   - When n is odd, Dn is computed on the faces of the grid cells
!     specified by the fillbox indices.  The index range for the 
!     undivided differences to be computed is ilo_fb to (ihi_fb+1); 
!     that is, the number of undivided difference computed is equal
!     to the number of faces associated with the fillbox grid cells
!     (ihi_fb - ilo_fb + 2).
!   - The ghostbox for Dn_minus_one MUST be at least one ghostcell width
!     larger than the fillbox.
!
!***********************************************************************
! toolbox/lsm_spatial_derivatives/lsm_spatial_derivatives2d.f
!***********************************************************************
	SUBROUTINE lsm2dComputeDn( &
		Dn, &
		ilo_Dn_gb, ihi_Dn_gb, &
		jlo_Dn_gb, jhi_Dn_gb, &
		Dn_minus_one, &
		ilo_Dn_minus_one_gb, ihi_Dn_minus_one_gb, &
		jlo_Dn_minus_one_gb, jhi_Dn_minus_one_gb, &
		ilo_fb, ihi_fb, jlo_fb, jhi_fb, &
		n, &
		dir, &
		Nthreads)

		! _gb refers to ghostbox 
		! _fb refers to fillbox 
		integer :: ilo_Dn_gb, ihi_Dn_gb, jlo_Dn_gb, jhi_Dn_gb
		integer :: ilo_Dn_minus_one_gb, ihi_Dn_minus_one_gb
		integer :: jlo_Dn_minus_one_gb, jhi_Dn_minus_one_gb
		integer :: ilo_fb, ihi_fb, jlo_fb, jhi_fb
		real(dp) :: Dn(ilo_Dn_gb:ihi_Dn_gb,jlo_Dn_gb:jhi_Dn_gb)
		real(dp) :: Dn_minus_one(ilo_Dn_minus_one_gb:ihi_Dn_minus_one_gb, &
			jlo_Dn_minus_one_gb:jhi_Dn_minus_one_gb)
		integer :: n, dir, i, j
		integer :: offset(1:2)
		integer :: fillbox_shift(1:2)
		real(dp) :: sign_multiplier
		real(dp), parameter :: big = 1.d10
		integer :: Nthreads

#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif

		! calculate offsets, fillbox shifts, and sign_multiplier used 
		! when computing undivided differences.
		! NOTE:  even and odd undivided differences are taken in
		!        opposite order because of the discrepancy between
		!        face- and cell-centered data.  the sign discrepancy 
		!        is taken into account by sign_multiplier
		do i = 1, 2
			offset(i) = 0
			fillbox_shift(i) = 0
		end do
		
		if (mod(n,2) == 1) then
			offset(dir) = 1
			sign_multiplier = 1.0
			fillbox_shift(dir) = 1
		else
			offset(dir) = -1
			sign_multiplier = -1.0
			fillbox_shift(dir) = 0
		endif

		! loop over cells with sufficient data {
#ifdef FLAG_FORALL
		forall(i=ilo_fb: ihi_fb+fillbox_shift(1), j=jlo_fb : jhi_fb+fillbox_shift(2) )
			Dn(i,j) = sign_multiplier * ( Dn_minus_one(i,j) - &
			Dn_minus_one(i-offset(1),j-offset(2)))
		end forall
#else
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_fb, jhi_fb+fillbox_shift(2)
			!do i = ilo_fb, ihi_fb+fillbox_shift(1)
			forall (i = ilo_fb : ihi_fb+fillbox_shift(1) )
				Dn(i,j) = sign_multiplier * ( Dn_minus_one(i,j) - &
				Dn_minus_one(i-offset(1),j-offset(2)))
			end forall
			!end do
		end do 
#endif
		! } end loop over grid

		! set undivided differences for cells with insufficient data to big {
#ifdef FLAG_FORALL
		forall(i = ilo_Dn_gb : ilo_fb-1, j = jlo_Dn_gb : jhi_Dn_gb)
			Dn(i,j) = big
		end forall
		
		forall(i = ihi_fb+fillbox_shift(1)+1 : ihi_Dn_gb, j = jlo_Dn_gb : jhi_Dn_gb)
			Dn(i,j) = big
		end forall
		
		forall(i = ilo_Dn_gb : ihi_Dn_gb, j = jlo_Dn_gb : jlo_fb-1)
			Dn(i,j) = big
		end forall
		
		forall(i = ilo_Dn_gb : ihi_Dn_gb, j = jhi_fb+fillbox_shift(2)+1 : jhi_Dn_gb)
			Dn(i,j) = big
		end forall
#else
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_Dn_gb, jhi_Dn_gb
		
			!do i = ilo_Dn_gb, ilo_fb-1
			forall(i = ilo_Dn_gb : ilo_fb-1)
				Dn(i,j) = big
			end forall
			!end do
		
			!do i = ihi_fb+fillbox_shift(1)+1, ihi_Dn_gb
			forall(i = ihi_fb+fillbox_shift(1)+1 : ihi_Dn_gb)
				Dn(i,j) = big
			end forall
			!end do
			
		end do
		
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_Dn_gb, jlo_fb-1
			!do i = ilo_Dn_gb, ihi_Dn_gb
			forall(i = ilo_Dn_gb : ihi_Dn_gb)
				Dn(i,j) = big
			end forall
			!end do
		end do

		!$OMP PARALLEL DO PRIVATE(i)
		do j = jhi_fb+fillbox_shift(2)+1, jhi_Dn_gb
			!do i = ilo_Dn_gb, ihi_Dn_gb
			forall(i = ilo_Dn_gb : ihi_Dn_gb)
				Dn(i,j) = big
			end forall
			!end do
		end do
#endif
		! } end setting big value for cells near boundary of ghostcell box

	END SUBROUTINE lsm2dComputeDn
      
    
    
!***********************************************************************
!
!  lsm2dHJENO1() computes the forward (plus) and backward (minus)
!  first-order Hamilton-Jacobi ENO approximations to the gradient of 
!  phi.
!
!  Arguments:
!    phi_*_plus (out):   components of grad(phi) in plus direction
!    phi_*_minus (out):  components of grad(phi) in minus direction
!    phi (in):           phi
!    D1 (in):            scratch space for holding undivided first-differences
!    dx, dy (in):        grid spacing
!    *_gb (in):          index range for ghostbox
!    *_fb (in):          index range for fillbox
!
!  NOTES:
!   - it is assumed that BOTH the plus AND minus derivatives have
!     the same fillbox
!
!***********************************************************************
! toolbox/lsm_spatial_derivatives/lsm_spatial_derivatives2d.f
!***********************************************************************
	SUBROUTINE lsm2dHJENO1( &
		phi_x_plus, phi_y_plus, &
		ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, &
		jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb, &
		phi_x_minus, phi_y_minus, &
		ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, &
		jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, &
		phi, &
		ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb, &
		D1, &
		ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb, &
		ilo_fb, ihi_fb, jlo_fb, jhi_fb, &
		dx, dy, Nthreads) &
		BIND(C, name = "LSM2_LSM2D_HJ_ENO1")

		! _grad_phi_plus_gb refers to ghostbox for grad_phi plus data
		! _grad_phi_minus_gb refers to ghostbox for grad_phi minus data
		! _phi_gb refers to ghostbox for phi data
		! _fb refers to fill-box for grad_phi data
		integer :: ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
		integer :: jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
		integer :: ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
		integer :: jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
		integer :: ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
		integer :: ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
		integer :: ilo_fb, ihi_fb, jlo_fb, jhi_fb
		real(dp) :: phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb, &
			jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
		real(dp) :: phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb, &
			jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
		real(dp) :: phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb, &
			jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
		real(dp) :: phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb, &
			jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
		real(dp) :: phi(ilo_phi_gb:ihi_phi_gb, jlo_phi_gb:jhi_phi_gb)
		real(dp) :: D1(ilo_D1_gb:ihi_D1_gb, jlo_D1_gb:jhi_D1_gb)
		real(dp) :: dx, dy, inv_dx, inv_dy
		integer :: i, j, Nthreads
		
		integer, parameter :: order = 1
		integer, parameter :: x_dir = 1
		integer, parameter :: y_dir = 2
		
#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif

		! compute inv_dx and inv_dy
		inv_dx = 1.0d0/dx
		inv_dy = 1.0d0/dy

		!----------------------------------------------------
		! compute phi_x_plus and phi_x_minus
		!----------------------------------------------------

		! compute first undivided differences in x-direction
		call lsm2dComputeDn(D1, &
			ilo_D1_gb, ihi_D1_gb, &
			jlo_D1_gb, jhi_D1_gb, &
			phi, &
			ilo_phi_gb, ihi_phi_gb, &
			jlo_phi_gb, jhi_phi_gb, &
			ilo_fb, ihi_fb, &
			jlo_fb, jhi_fb, &
			order, x_dir, Nthreads)

		! { begin loop over grid
#ifdef FLAG_FORALL
		forall(i = ilo_fb : ihi_fb, j = jlo_fb : jhi_fb)
			phi_x_plus(i,j) = D1(i+1,j)*inv_dx
			phi_x_minus(i,j) = D1(i,j)*inv_dx
		end forall
#else
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_fb, jhi_fb
			!do i = ilo_fb, ihi_fb
			forall(i = ilo_fb : ihi_fb)
				phi_x_plus(i,j) = D1(i+1,j)*inv_dx
				phi_x_minus(i,j) = D1(i,j)*inv_dx
			end forall
			!end do
		end do 
#endif
		! } end loop over grid 

		!----------------------------------------------------
		! compute phi_y_plus and phi_y_minus
		!----------------------------------------------------

		! compute first undivided differences in y-direction
		call lsm2dComputeDn(D1, &
			ilo_D1_gb, ihi_D1_gb, &
			jlo_D1_gb, jhi_D1_gb, &
			phi, &
			ilo_phi_gb, ihi_phi_gb, &
			jlo_phi_gb, jhi_phi_gb, &
			ilo_fb, ihi_fb, &
			jlo_fb, jhi_fb, &
			order, y_dir, Nthreads)

		! { begin loop over grid
#ifdef FLAG_FORALL
		forall(i = ilo_fb : ihi_fb, j = jlo_fb : jhi_fb)
			phi_y_plus(i,j) = D1(i,j+1)*inv_dy
			phi_y_minus(i,j) = D1(i,j)*inv_dy
		end forall
#else
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_fb, jhi_fb
			!do i = ilo_fb, ihi_fb
			forall(i = ilo_fb : ihi_fb)
				phi_y_plus(i,j) = D1(i,j+1)*inv_dy
				phi_y_minus(i,j) = D1(i,j)*inv_dy
			end forall
			!end do
		end do
#endif
		! } end loop over grid 

      END SUBROUTINE lsm2dHJENO1
      
 
!***********************************************************************
!
!  lsm2dAverageGradPhi() computes the average of the plus and minus 
!  derivatives:
!
!    phi_* = (phi_*_plus + phi_*_minus) / 2
!
!  Arguments:
!    phi_* (out):       components of average grad(phi)
!    phi_*_plus (in):   components of grad(phi) in plus direction
!    phi_*_minus (in):  components of grad(phi) in minus direction
!    *_gb (in):         index range for ghostbox
!    *_fb (in):         index range for fillbox
!
!***********************************************************************
! toolbox/lsm_spatial_derivatives/lsm_spatial_derivatives2d.f
!***********************************************************************
	SUBROUTINE lsm2dAverageGradPhi( &
		phi_x, phi_y, &
		ilo_grad_phi_gb, ihi_grad_phi_gb, &
		jlo_grad_phi_gb, jhi_grad_phi_gb, &
		phi_x_plus, phi_y_plus, &
		ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, &
		jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb, &
		phi_x_minus, phi_y_minus, &
		ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, &
		jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, &
		ilo_fb, ihi_fb, &
		jlo_fb, jhi_fb, &
		Nthreads) &
		BIND(C, name = "LSM2_LSM2D_AVERAGE_GRAD_PHI")
		
!      _gb refers to ghostbox 
!      _fb refers to fill-box

		integer :: ilo_grad_phi_gb, ihi_grad_phi_gb
		integer :: jlo_grad_phi_gb, jhi_grad_phi_gb
		integer :: ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
		integer :: jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
		integer :: ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
		integer :: jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
		integer :: ilo_fb, ihi_fb
		integer :: jlo_fb, jhi_fb
		real(dp) :: phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb, &
			jlo_grad_phi_gb:jhi_grad_phi_gb)
		real(dp) :: phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb, &
			jlo_grad_phi_gb:jhi_grad_phi_gb)
		real(dp) :: phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb, &
			jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
		real(dp) :: phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb, &
			jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
		real(dp) :: phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb, &
			jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
		real(dp) :: phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb, &
			jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
		integer :: i, j
		real(dp), parameter :: half = 0.5d0
		integer :: Nthreads
		
#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif
		
		! compute "phi-upwind" derivatives
		
		! { begin loop over grid
#ifdef FLAG_FORALL
		forall(i = ilo_fb : ihi_fb, j = jlo_fb : jhi_fb)
			phi_x(i,j) = half * ( phi_x_plus(i,j) + phi_x_minus(i,j) )
			phi_y(i,j) = half * ( phi_y_plus(i,j) + phi_y_minus(i,j) )
		end forall
#else
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_fb, jhi_fb
			!do i = ilo_fb, ihi_fb
			forall(i = ilo_fb : ihi_fb)
				phi_x(i,j) = half * ( phi_x_plus(i,j) + phi_x_minus(i,j) )
				phi_y(i,j) = half * ( phi_y_plus(i,j) + phi_y_minus(i,j) )
			end forall
			!end do
		end do
#endif
		! } end loop over grid
		
	END SUBROUTINE lsm2dAverageGradPhi
	
	
!***********************************************************************
!
!  lsm2dComputeUnitNormal() computes the unit normal vector to the 
!  interface from grad(phi). 
!
!  Arguments:
!    normal_* (out):  components of unit normal vector
!    phi_* (in):      components of grad(phi) 
!    dx, dy (in):     grid spacing
!    *_gb (in):       index range for ghostbox
!    *_fb (in):       index range for fillbox
!
!***********************************************************************
! toolbox/geometry/lsm_geometry2d.f
!***********************************************************************
	SUBROUTINE lsm2dComputeUnitNormal( &
		normal_x, normal_y, &
		ilo_normal_gb, ihi_normal_gb, &
		jlo_normal_gb, jhi_normal_gb, &
		phi_x, phi_y, &
		ilo_grad_phi_gb, ihi_grad_phi_gb, &
		jlo_grad_phi_gb, jhi_grad_phi_gb, &
		ilo_fb, ihi_fb, &
		jlo_fb, jhi_fb, &
		Nthreads) &
		BIND(C, name = "LSM2_LSM2D_COMPUTE_UNIT_NORMAL")
		
		!     _gb refers to ghostboxes 
		!     _fb refers to fill-box for normal data
		
		integer :: ilo_normal_gb, ihi_normal_gb
		integer :: jlo_normal_gb, jhi_normal_gb
		integer :: ilo_grad_phi_gb, ihi_grad_phi_gb
		integer :: jlo_grad_phi_gb, jhi_grad_phi_gb
		integer :: ilo_fb, ihi_fb
		integer :: jlo_fb, jhi_fb
		real(dp) :: normal_x(ilo_normal_gb:ihi_normal_gb, &
			jlo_normal_gb:jhi_normal_gb)
		real(dp) :: normal_y(ilo_normal_gb:ihi_normal_gb, &
			jlo_normal_gb:jhi_normal_gb)
		real(dp) :: phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb, &
			jlo_grad_phi_gb:jhi_grad_phi_gb)
		real(dp) :: phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb, &
			jlo_grad_phi_gb:jhi_grad_phi_gb)
		real(dp) :: norm_grad_phi, inv_norm_grad_phi
		integer :: i, j
		real(dp), parameter :: half = 0.5d0
		real(dp), parameter :: zero_tol = 1.d-11
		integer :: Nthreads
		
#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif
		
		!$OMP PARALLEL DO PRIVATE(i, norm_grad_phi, inv_norm_grad_phi)
		! { begin loop over grid
		do j = jlo_fb, jhi_fb
			do i = ilo_fb, ihi_fb
				! compute unit normal 
				norm_grad_phi = sqrt( phi_x(i,j)*phi_x(i,j) + phi_y(i,j)*phi_y(i,j) )
				if (norm_grad_phi .ge. zero_tol) then
					inv_norm_grad_phi = 1.0d0/norm_grad_phi
					normal_x(i,j) = phi_x(i,j)*inv_norm_grad_phi
					normal_y(i,j) = phi_y(i,j)*inv_norm_grad_phi
				else
					normal_x(i,j) = 1.0d0
					normal_y(i,j) = 0.0d0
				end if
			end do
		end do
		! } end loop over grid
		
	END SUBROUTINE lsm2dComputeUnitNormal
	
	

!***********************************************************************
! toolbox/lsm_level_set_evolution/lsm_level_set_evolution2d.f
!***********************************************************************
	SUBROUTINE lsm2dAddNormalVelTermToLSERHS( &
		lse_rhs, &
		ilo_lse_rhs_gb, ihi_lse_rhs_gb, &
		jlo_lse_rhs_gb, jhi_lse_rhs_gb, &
		phi_x_plus, phi_y_plus, &
		ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, &
		jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb, &
		phi_x_minus, phi_y_minus, &
		ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, &
		jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, &
		vel_n, &
		ilo_vel_gb, ihi_vel_gb, &
		jlo_vel_gb, jhi_vel_gb, &
		ilo_fb, ihi_fb, &
		jlo_fb, jhi_fb, &
		Nthreads ) &
		BIND(C, name = "LSM2_LSM2D_ADD_NORMAL_VEL_TERM_TO_LSE_RHS")
		
		integer :: ilo_lse_rhs_gb, ihi_lse_rhs_gb
		integer :: jlo_lse_rhs_gb, jhi_lse_rhs_gb
		integer :: ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
		integer :: jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
        integer :: ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
        integer :: jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
        integer :: ilo_vel_gb, ihi_vel_gb
        integer :: jlo_vel_gb, jhi_vel_gb
        integer :: ilo_fb, ihi_fb
        integer :: jlo_fb, jhi_fb
        real(dp) :: lse_rhs(ilo_lse_rhs_gb:ihi_lse_rhs_gb, &
			jlo_lse_rhs_gb:jhi_lse_rhs_gb)
        real(dp) :: phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb, &
			jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
        real(dp) :: phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb, &
			jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
        real(dp) :: phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb, &
			jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
		real(dp) :: phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb, &
			jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
		real(dp) :: vel_n(ilo_vel_gb:ihi_vel_gb, jlo_vel_gb:jhi_vel_gb)
		integer :: i, j
		real(dp) :: vel_n_cur, norm_grad_phi_sq
		real(dp), parameter :: zero_tol = 1.d-11
		integer :: Nthreads

#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif
		
		! { begin loop over grid
		!$OMP PARALLEL DO PRIVATE(i, vel_n_cur, norm_grad_phi_sq)	
		do j = jlo_fb, jhi_fb
			do i = ilo_fb, ihi_fb
				vel_n_cur = vel_n(i,j)
				if (abs(vel_n_cur) >= zero_tol) then
					! { begin Godunov selection of grad_phi
					if (vel_n_cur >= 0.d0) then
						norm_grad_phi_sq = max(max(phi_x_minus(i,j),0.d0)**2, &
							min(phi_x_plus(i,j),0.d0)**2 ) + &
							max(max(phi_y_minus(i,j),0.d0)**2, &
							min(phi_y_plus(i,j),0.d0)**2 )
					else
						norm_grad_phi_sq = max(min(phi_x_minus(i,j),0.d0)**2, &
							max(phi_x_plus(i,j),0.d0)**2 ) + &
							max(min(phi_y_minus(i,j),0.d0)**2, &
							max(phi_y_plus(i,j),0.d0)**2 )
					end if
					! } end Godunov selection of grad_phi
					! compute contribution to lse_rhs(i,j) 
					lse_rhs(i,j) = lse_rhs(i,j) - vel_n_cur*sqrt(norm_grad_phi_sq)
				end if
			end do 
		end do 
		! } end loop over grid
	
	END SUBROUTINE lsm2dAddNormalVelTermToLSERHS
	
	
!***********************************************************************
!
!  lsm2dRK1Step() takes a single first-order Runge-Kutta (i.e. Forward 
!  Euler) step.
!  
!  Arguments:
!    u_next (out):  u(t_cur+dt)
!    u_cur (in):    u(t_cur)
!    rhs (in):      right-hand side of time evolution equation
!    dt (in):       step size
!    *_gb (in):     index range for ghostbox
!    *_fb (in):     index range for fillbox
!
!***********************************************************************
! toolbox/time_integration/lsm_tvd_runge_kutta2d.f
!***********************************************************************
	SUBROUTINE lsm2dRK1Step( &
		u_next, &
		ilo_u_next_gb, ihi_u_next_gb, &
		jlo_u_next_gb, jhi_u_next_gb, &
		u_cur, &
		ilo_u_cur_gb, ihi_u_cur_gb, &
		jlo_u_cur_gb, jhi_u_cur_gb, &
		rhs, &
		ilo_rhs_gb, ihi_rhs_gb, &
		jlo_rhs_gb, jhi_rhs_gb, &
		ilo_fb, ihi_fb, jlo_fb, jhi_fb, &
		dt, &
		Nthreads)
		
		integer :: ilo_u_next_gb, ihi_u_next_gb
		integer :: jlo_u_next_gb, jhi_u_next_gb
		integer :: ilo_u_cur_gb, ihi_u_cur_gb
		integer :: jlo_u_cur_gb, jhi_u_cur_gb
		integer :: ilo_rhs_gb, ihi_rhs_gb
		integer :: jlo_rhs_gb, jhi_rhs_gb
		integer :: ilo_fb, ihi_fb, jlo_fb, jhi_fb
		real(dp) :: u_next(ilo_u_next_gb:ihi_u_next_gb, &
			jlo_u_next_gb:jhi_u_next_gb)
		real(dp) :: u_cur(ilo_u_cur_gb:ihi_u_cur_gb, jlo_u_cur_gb:jhi_u_cur_gb)
		real(dp) :: rhs(ilo_rhs_gb:ihi_rhs_gb, jlo_rhs_gb:jhi_rhs_gb)
		integer :: i, j
		real(dp) :: dt
		integer :: Nthreads
		
#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif
		
		! { begin loop over grid
#ifdef FLAG_FORALL
		forall(i = ilo_fb : ihi_fb, j = jlo_fb : jhi_fb)
			u_next(i,j) = u_cur(i,j) + dt*rhs(i,j)
		end forall
#else
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_fb, jhi_fb
			!do i = ilo_fb, ihi_fb
			forall(i = ilo_fb : ihi_fb)
				u_next(i,j) = u_cur(i,j) + dt*rhs(i,j)
			end forall
			!end do
		end do
#endif
		! } end loop over grid
	
	END SUBROUTINE lsm2dRK1Step
      
      
!***********************************************************************
!
!  lsm2dTVDRK2Stage1() advances the solution through first stage of the 
!  second-order TVD Runge-Kutta step.
!  
!  Arguments:
!    u_stage1 (out):  u_approx(t_cur+dt)
!    u_cur (in):      u(t_cur)
!    rhs (in):        right-hand side of time evolution equation
!    dt (in):         step size
!    *_gb (in):       index range for ghostbox
!    *_fb (in):       index range for fillbox
!
!  NOTES:
!   - the first stage of TVD RK2 is identical to a single RK1 step
!
!***********************************************************************
! toolbox/time_integration/lsm_tvd_runge_kutta2d.f
!***********************************************************************
	SUBROUTINE lsm2dTVDRK2Stage1( &
		u_stage1, &
		ilo_u_stage1_gb, ihi_u_stage1_gb, &
		jlo_u_stage1_gb, jhi_u_stage1_gb, &
		u_cur, &
		ilo_u_cur_gb, ihi_u_cur_gb, &
		jlo_u_cur_gb, jhi_u_cur_gb, &
		rhs, &
		ilo_rhs_gb, ihi_rhs_gb, &
		jlo_rhs_gb, jhi_rhs_gb, &
		ilo_fb, ihi_fb, jlo_fb, jhi_fb, &
		dt, &
		Nthreads) &
		BIND(C, name = "LSM2_LSM2D_TVD_RK2_STAGE1")
		
		integer :: ilo_u_stage1_gb, ihi_u_stage1_gb
		integer :: jlo_u_stage1_gb, jhi_u_stage1_gb
		integer :: ilo_u_cur_gb, ihi_u_cur_gb
		integer :: jlo_u_cur_gb, jhi_u_cur_gb
		integer :: ilo_rhs_gb, ihi_rhs_gb
		integer :: jlo_rhs_gb, jhi_rhs_gb
		integer :: ilo_fb, ihi_fb, jlo_fb, jhi_fb
		real(dp) :: u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb, &
			jlo_u_stage1_gb:jhi_u_stage1_gb)
		real(dp) :: u_cur(ilo_u_cur_gb:ihi_u_cur_gb, jlo_u_cur_gb:jhi_u_cur_gb)
		real(dp) :: rhs(ilo_rhs_gb:ihi_rhs_gb, jlo_rhs_gb:jhi_rhs_gb)
		real(dp) :: dt
		integer :: Nthreads
		
		! use lsm2dRK1Step() to compute first stage
		call lsm2dRK1Step(u_stage1, &
			ilo_u_stage1_gb, ihi_u_stage1_gb, &
			jlo_u_stage1_gb, jhi_u_stage1_gb, &
			u_cur, &
			ilo_u_cur_gb, ihi_u_cur_gb, &
			jlo_u_cur_gb, jhi_u_cur_gb, &
			rhs, &
			ilo_rhs_gb, ihi_rhs_gb, &
			jlo_rhs_gb, jhi_rhs_gb, &
			ilo_fb, ihi_fb, jlo_fb, jhi_fb, &
			dt, &
			Nthreads)
	
	END SUBROUTINE lsm2dTVDRK2Stage1



!***********************************************************************
!
!  lsm2dTVDRK2Stage2() completes advancing the solution through a 
!  single step of the second-order TVD Runge-Kutta step.
!  
!  Arguments:
!    u_next (out):   u(t_cur+dt)
!    u_stage1 (in):  u_approx(t_cur+dt)
!    u_cur (in):     u(t_cur)
!    rhs (in):       right-hand side of time evolution equation
!    dt (in):        step size
!    *_gb (in):      index range for ghostbox
!    *_fb (in):      index range for fillbox
!
!***********************************************************************
! toolbox/time_integration/lsm_tvd_runge_kutta2d.f
!***********************************************************************
	SUBROUTINE lsm2dTVDRK2Stage2( &
		u_next, &
		ilo_u_next_gb, ihi_u_next_gb, &
		jlo_u_next_gb, jhi_u_next_gb, &
		u_stage1, &
		ilo_u_stage1_gb, ihi_u_stage1_gb, &
		jlo_u_stage1_gb, jhi_u_stage1_gb, &
		u_cur, &
		ilo_u_cur_gb, ihi_u_cur_gb, &
		jlo_u_cur_gb, jhi_u_cur_gb, &
		rhs, &
		ilo_rhs_gb, ihi_rhs_gb, &
		jlo_rhs_gb, jhi_rhs_gb, &
		ilo_fb, ihi_fb, jlo_fb, jhi_fb, &
		dt, &
		Nthreads) &
		BIND(C, name = "LSM2_LSM2D_TVD_RK2_STAGE2")
		
		integer :: ilo_u_next_gb, ihi_u_next_gb
		integer :: jlo_u_next_gb, jhi_u_next_gb
		integer :: ilo_u_stage1_gb, ihi_u_stage1_gb
		integer :: jlo_u_stage1_gb, jhi_u_stage1_gb
		integer :: ilo_u_cur_gb, ihi_u_cur_gb
		integer :: jlo_u_cur_gb, jhi_u_cur_gb
		integer :: ilo_rhs_gb, ihi_rhs_gb
		integer :: jlo_rhs_gb, jhi_rhs_gb
		integer :: ilo_fb, ihi_fb, jlo_fb, jhi_fb
		real(dp) :: u_next(ilo_u_next_gb:ihi_u_next_gb, &
			jlo_u_next_gb:jhi_u_next_gb)
		real(dp) :: u_stage1(ilo_u_stage1_gb:ihi_u_stage1_gb, &
			jlo_u_stage1_gb:jhi_u_stage1_gb)
		real(dp) :: u_cur(ilo_u_cur_gb:ihi_u_cur_gb, &
			jlo_u_cur_gb:jhi_u_cur_gb)
		real(dp) :: rhs(ilo_rhs_gb:ihi_rhs_gb, jlo_rhs_gb:jhi_rhs_gb)
		integer :: i, j
		real(dp) :: dt
		integer :: Nthreads

#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif
		
		! { begin loop over grid
#ifdef FLAG_FORALL
		forall(i = ilo_fb : ihi_fb, j = jlo_fb : jhi_fb)
			u_next(i,j) = 0.5d0*( u_cur(i,j) + u_stage1(i,j) + dt*rhs(i,j) )
		end forall
#else
		!$OMP PARALLEL DO PRIVATE(i)
		do j = jlo_fb, jhi_fb
			!do i = ilo_fb, ihi_fb
			forall(i = ilo_fb : ihi_fb)
				u_next(i,j) = 0.5d0*( u_cur(i,j) + u_stage1(i,j) + dt*rhs(i,j) )
			end forall
			!end do
		end do
#endif
		! } end loop over grid
		
	END SUBROUTINE lsm2dTVDRK2Stage2
	
	
!***********************************************************************
!
! lsm2dSignedLinearExtrapolation() extrapolates 2D data from the index 
! range of the fillbox into cells in ghostbox at the specified boundary 
! location using signed linear extrapolation.  Extrapolation is "away" 
! from the zero level set, i.e. the zero level set will not be 
! artifically created at the boundary.
!
! Arguments:
!   phi (in/out):            phi
!   bdry_location_idx (in):  boundary location index
!   *_gb (in):               index range for ghostbox
!   *_fb (in):               index range for fillbox
! 
! NOTES:
!  - fillbox indices must be a subset of ghostbox indices.
!  - if bdry_location_idx is out of the range for 2D, then no 
!    ghostcell values are set 
!
!***********************************************************************
! serial/boundary_conditions/lsm_boundary_conditions.c
! toolbox/boundary_conditions/lsm_boundary_conditions2d.f
!***********************************************************************
	SUBROUTINE lsm2dSignedLinearExtrapolation( &
		phi, &
		ilo_gb, ihi_gb, jlo_gb, jhi_gb, &
		ilo_fb, ihi_fb, jlo_fb, jhi_fb, &
		bdry_location_idx, &
		Nthreads)
		
		integer :: ilo_gb, ihi_gb, jlo_gb, jhi_gb
		integer :: ilo_fb, ihi_fb, jlo_fb, jhi_fb
		integer :: bdry_location_idx
		real :: phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
		integer :: Nthreads
				
		! local variables
		integer :: i, j
		integer, parameter :: zero = 0
		real :: s, abs_diff, dist, slope
		real, parameter :: one = 1.0d0
		
		! { extrapolate data in x-direction at lower end
		if (bdry_location_idx .eq. 0) then
			! { begin j loop
			do j = jlo_gb, jhi_gb
				s = sign(one,phi(ilo_fb,j))
				abs_diff = abs(phi(ilo_fb,j) - phi(ilo_fb+1,j))
				slope = s*abs_diff
				do i = ilo_gb, ilo_fb-1
					dist = ilo_fb - i
					phi(i,j) = phi(ilo_fb,j) + slope*dist
				end do
			end do  
			! } end j loop
		! } end extrapolate data in x-direction at lower end
			
		! { extrapolate data in x-direction at upper end
		elseif (bdry_location_idx .eq. 1) then	
			!  { begin j loop
			do j = jlo_gb, jhi_gb
				s = sign(one,phi(ihi_fb,j))
				abs_diff = abs(phi(ihi_fb,j) - phi(ihi_fb-1,j))
				slope = s*abs_diff
				do i = ihi_fb+1, ihi_gb
					dist = i - ihi_fb
					phi(i,j) = phi(ihi_fb,j) + slope*dist
				end do 
			end do  
			! } end j loop
		! } end extrapolate data in x-direction at upper end
		
		! { extrapolate data in y-direction at lower end
		elseif (bdry_location_idx .eq. 2) then
			! { begin i loop
			do i = ilo_gb, ihi_gb
				s = sign(one,phi(i,jlo_fb))
				abs_diff = abs(phi(i,jlo_fb) - phi(i,jlo_fb+1))
				slope = s*abs_diff
				do j = jlo_gb, jlo_fb-1
					dist = jlo_fb - j 
					phi(i,j) = phi(i,jlo_fb) + slope*dist
				end do
			end do 
			! } end i loop
		! } end extrapolate data in y-direction at lower end
		
		! { extrapolate data in y-direction at upper end
		elseif (bdry_location_idx .eq. 3) then
			! { begin i loop
			do i = ilo_gb, ihi_gb
				s = sign(one,phi(i,jhi_fb))
				abs_diff = abs(phi(i,jhi_fb) - phi(i,jhi_fb-1))
				slope = s*abs_diff
				do j = jhi_fb+1, jhi_gb
					dist = j - jhi_fb
					phi(i,j) = phi(i,jhi_fb) + slope*dist
				end do 
			end do 
			! } end i loop
		! } end extrapolate data in y-direction at lower end
		
		end if
		
	END SUBROUTINE lsm2dSignedLinearExtrapolation


!***********************************************************************
! versione Undy (26.10.2013)
!***********************************************************************

SUBROUTINE lsm2dSignedLinearExtrapolation9( &
		phi, &
		ilo_gb, ihi_gb, jlo_gb, jhi_gb, &
		ilo_fb, ihi_fb, jlo_fb, jhi_fb, &
		Nthreads) &
		BIND(C, name = "LSM2_signedLinearExtrapolationBC")
		
		integer :: ilo_gb, ihi_gb, jlo_gb, jhi_gb
		integer :: ilo_fb, ihi_fb, jlo_fb, jhi_fb
		real :: phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
		integer :: Nthreads
				
		! local variables
		integer :: i, j
		integer, parameter :: zero = 0
		real :: s, abs_diff, dist, slope
		real, parameter :: one = 1.0d0
		
		
#ifdef _OPENMP
        call omp_set_num_threads(Nthreads)
#endif
		
		! { begin j loop
		!$OMP PARALLEL DO PRIVATE(i, s, abs_diff, slope, dist)
		do j = jlo_gb, jhi_gb
		
			! { extrapolate data in x-direction at lower end
			s = sign(one,phi(ilo_fb,j))
			abs_diff = abs(phi(ilo_fb,j) - phi(ilo_fb+1,j))
			slope = s*abs_diff
			!do i = ilo_gb, ilo_fb-1
			forall(i = ilo_gb : ilo_fb-1)
				!dist = ilo_fb - i
				phi(i,j) = phi(ilo_fb,j) + slope * (ilo_fb - i) !dist
			end forall
			!end do
			! } end extrapolate data in x-direction at lower end
			
			! { extrapolate data in x-direction at upper end
			s = sign(one,phi(ihi_fb,j))
			abs_diff = abs(phi(ihi_fb,j) - phi(ihi_fb-1,j))
			slope = s*abs_diff
			!do i = ihi_fb+1, ihi_gb
			forall(i = ihi_fb+1 : ihi_gb)
				!dist = i - ihi_fb
				phi(i,j) = phi(ihi_fb,j) + slope * (i - ihi_fb) !dist
			end forall
			!end do 
			! } end extrapolate data in x-direction at upper end
			
		end do 
		! } end j loop
		
		! { begin i loop
		!$OMP PARALLEL DO PRIVATE(j, s, abs_diff, slope, dist)
		do i = ilo_gb, ihi_gb
		
			! { extrapolate data in y-direction at lower end
			s = sign(one,phi(i,jlo_fb))
			abs_diff = abs(phi(i,jlo_fb) - phi(i,jlo_fb+1))
			slope = s*abs_diff
			!do j = jlo_gb, jlo_fb-1
			forall(j = jlo_gb : jlo_fb-1)
				!dist = jlo_fb - j 
				phi(i,j) = phi(i,jlo_fb) + slope * (jlo_fb - j) !dist
			end forall
			!end do
			! } end extrapolate data in y-direction at lower end

			! { extrapolate data in y-direction at upper end
			s = sign(one,phi(i,jhi_fb))
			abs_diff = abs(phi(i,jhi_fb) - phi(i,jhi_fb-1))
			slope = s*abs_diff
			!do j = jhi_fb+1, jhi_gb
			forall(j = jhi_fb+1 : jhi_gb)
				!dist = j - jhi_fb
				phi(i,j) = phi(i,jhi_fb) + slope * (j - jhi_fb) !dist
			end forall
			!end do 
			! } end extrapolate data in y-direction at lower end
		
		end do 
		! } end i loop
			
	END SUBROUTINE lsm2dSignedLinearExtrapolation9
    
END MODULE LSM2
