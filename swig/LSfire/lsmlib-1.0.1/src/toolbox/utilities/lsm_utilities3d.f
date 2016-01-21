c***********************************************************************
c
c  File:        lsm_utilities3d.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 routines for 3D level set method utility subroutines
c
c***********************************************************************



c***********************************************************************
      subroutine lsm3dMaxNormDiff(
     &  max_norm_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  klo_field1_gb, khi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  klo_field2_gb, khi_field2_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer klo_field1_gb, khi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer klo_field2_gb, khi_field2_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb,
     &            klo_field1_gb:khi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb,
     &            klo_field2_gb:khi_field2_gb)
      real max_norm_diff
      real next_diff
      integer i,j,k


c     initialize max_norm_diff
      max_norm_diff = abs( field1(ilo_ib,jlo_ib,klo_ib) 
     &                   - field2(ilo_ib,jlo_ib,klo_ib))

c       loop over included cells { 
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                next_diff = abs(field1(i,j,k) - field2(i,j,k))
                if (next_diff .gt. max_norm_diff) then
                  max_norm_diff = next_diff
                endif
 
            enddo
          enddo
        enddo
c       } end loop over grid 
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dAveAbsDiff(
     &  ave_abs_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  klo_field1_gb, khi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  klo_field2_gb, khi_field2_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      real ave_abs_diff
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer klo_field1_gb, khi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer klo_field2_gb, khi_field2_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib   
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb,
     &            klo_field1_gb:khi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb,
     &            klo_field2_gb:khi_field2_gb)
      
c     local variables      
      real sum_abs_diff, num_pts, next_diff
      real zero, one
      parameter (zero=0.d0, one=1.d0)
      integer i,j,k

c     initialize max_norm_diff
      sum_abs_diff = zero
      num_pts = zero

c       loop over grid { 
        do k=klo_ib,khi_ib
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	      next_diff = abs(field1(i,j,k) - field2(i,j,k))	      	      
              sum_abs_diff = sum_abs_diff + next_diff
              num_pts = num_pts + one
	      
          enddo
	 enddo 
        enddo
c       } end loop over grid 
      
      if( num_pts .gt. zero) then
         ave_abs_diff = sum_abs_diff / num_pts
      else 
         ave_abs_diff = zero
      endif
      	 
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dVoxelCountGreaterThanZero(
     &  count,
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      integer count
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
c     local variables      
      integer i,j,k

      count = 0
c       loop over grid { 
        do k=klo_ib,khi_ib
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j,k) .gt. 0.d0) then
	       count = count + 1
	     endif
	     
	  enddo       
         enddo
        enddo
c       } end loop over grid     
      	 
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dVoxelCountLessThanZero(
     &  count,
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      integer count
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      
c     local variables      
      integer i,j,k
      
      count = 0
c       loop over grid { 
        do k=klo_ib,khi_ib
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j,k) .lt. 0.d0) then
	       count = count + 1
	     endif
	        
          enddo
	 enddo 
        enddo
c       } end loop over grid     
      	 
      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm3dComputeStableAdvectionDt(
     &  dt,
     &  vel_x, vel_y, vel_z,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  klo_vel_gb, khi_vel_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer klo_vel_gb, khi_vel_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb,
     &           klo_vel_gb:khi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb,
     &           klo_vel_gb:khi_vel_gb)
      real vel_z(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb,
     &           klo_vel_gb:khi_vel_gb)
      real dx, dy, dz
      real inv_dx, inv_dy, inv_dz
      real cfl_number
      integer i, j, k
      real max_U_over_dX
      real U_over_dX_cur
      real small_number
      parameter (small_number = 1.d-99)

c     initialize max_U_over_dX to -1
      max_U_over_dX = -1.0d0

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy
      inv_dz = 1.d0/dz
  
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                U_over_dX_cur = abs(vel_x(i,j,k))*inv_dx
     &                        + abs(vel_y(i,j,k))*inv_dy
     &                        + abs(vel_z(i,j,k))*inv_dz

                if (U_over_dX_cur .gt. max_U_over_dX) then
                  max_U_over_dX = U_over_dX_cur  
                endif

            enddo
          enddo
        enddo
c       } end loop over grid

c     set dt
      dt = cfl_number / (max_U_over_dX + small_number);

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dComputeStableNormalVelDt(
     &  dt,
     &  vel_n,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  klo_vel_gb, khi_vel_gb,
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz, 
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer klo_vel_gb, khi_vel_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real vel_n(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb,
     &           klo_vel_gb:khi_vel_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_z_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_z_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real dx,dy,dz
      real inv_dx, inv_dy, inv_dz
      real max_dx_sq
      real cfl_number
      integer i,j,k
      real max_H_over_dX
      real H_over_dX_cur
      real phi_x_cur, phi_y_cur, phi_z_cur
      real norm_grad_phi
      real small_number
      parameter (small_number = 1.d-99)

c     compute max_dx_sq
      max_dx_sq = max(dx,dy,dz)
      max_dx_sq = max(dx,dy,dz) * max(dx,dy,dz)

c     initialize max_H_over_dX to -1
      max_H_over_dX = -1.0d0

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy
      inv_dz = 1.d0/dz
      
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                  phi_x_cur = max(abs(phi_x_plus(i,j,k)),
     &                          abs(phi_x_minus(i,j,k)))
                  phi_y_cur = max(abs(phi_y_plus(i,j,k)),
     &                          abs(phi_y_minus(i,j,k)))
                  phi_z_cur = max(abs(phi_z_plus(i,j,k)),
     &                          abs(phi_z_minus(i,j,k)))
                  norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                              + phi_y_cur*phi_y_cur 
     &                              + phi_z_cur*phi_z_cur + max_dx_sq )

                  H_over_dX_cur = abs(vel_n(i,j,k)) / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy 
     &                        + phi_z_cur*inv_dz )

                if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
                endif
	      
            enddo
          enddo
        enddo
c       } end loop over grid
      
c     set dt
      dt = cfl_number / (max_H_over_dX + small_number);

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dComputeStableConstNormalVelDt(
     &  dt,
     &  vel_n,
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real vel_n
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_z_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_z_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real dx,dy,dz
      real inv_dx,inv_dy,inv_dz
      real max_dx_sq
      real cfl_number
      integer i,j,k
      real max_H_over_dX, abs_vel_n, norm_grad_phi
      real H_over_dX_cur
      real phi_x_cur, phi_y_cur, phi_z_cur
      real small_number
      parameter (small_number = 1.d-99)

c     compute max_dx_sq
      max_dx_sq = max(dx,dy,dz)
      max_dx_sq = max_dx_sq * max_dx_sq

c     initialize max_H_over_dX to -1
      max_H_over_dX = -1

c     compute inv_dx, inv_dy, inv_dz
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy
      inv_dz = 1.d0/dz
      
      abs_vel_n = abs(vel_n)

c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                phi_x_cur = max(abs(phi_x_plus(i,j,k)),
     &         	              abs(phi_x_minus(i,j,k)))
                phi_y_cur = max(abs(phi_y_plus(i,j,k)),
     &	                        abs(phi_y_minus(i,j,k)))
                phi_z_cur = max(abs(phi_z_plus(i,j,k)),
     &	                        abs(phi_z_minus(i,j,k)))

                norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                              + phi_y_cur*phi_y_cur 
     &                              + phi_z_cur*phi_z_cur + max_dx_sq )

                H_over_dX_cur = abs_vel_n / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                          + phi_y_cur*inv_dy 
     &                          + phi_z_cur*inv_dz )
     
  
                if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
                endif

            enddo
          enddo
        enddo
c       } end loop over grid
      
c     set dt
      dt = cfl_number / (max_H_over_dX + small_number);

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dVolumeIntegralPhiLessThanZero(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  klo_F_gb, khi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer klo_F_gb, khi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb,
     &       klo_F_gb:khi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real dx,dy,dz
      real epsilon
      integer i,j,k
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     initialize int_F to zero
      int_F = 0.0d0

c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .lt. -epsilon) then
                  int_F = int_F + F(i,j,k)*dV
                elseif (phi_cur .lt. epsilon) then
                  one_minus_H = 
     &		    0.5d0*(1.d0-phi_cur_over_epsilon
     &                   -one_over_pi*sin(pi*phi_cur_over_epsilon))
                  int_F = int_F + one_minus_H*F(i,j,k)*dV
                endif
    
            enddo
          enddo
        enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dVolumeIntegralPhiGreaterThanZero(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  klo_F_gb, khi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer klo_F_gb, khi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb,
     &       klo_F_gb:khi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real dx,dy,dz
      real epsilon
      integer i,j,k
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     initialize int_F to zero
      int_F = 0.0d0
     
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
  
               phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .gt. epsilon) then
                  int_F = int_F + F(i,j,k)*dV
                elseif (phi_cur .gt. -epsilon) then
                  H = 0.5d0*( 1.d0 + phi_cur_over_epsilon 
     &                             + one_over_pi
     &                             * sin(pi*phi_cur_over_epsilon) )
                  int_F = int_F + H*F(i,j,k)*dV
                endif
          
            enddo
          enddo
        enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dSurfaceIntegral(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  klo_F_gb, khi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  phi_x, phi_y, phi_z,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  klo_grad_phi_gb, khi_grad_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer klo_F_gb, khi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb,
     &       klo_F_gb:khi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_z(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real dx,dy,dz
      real epsilon
      real one_over_epsilon
      integer i,j,k
      real phi_cur
      real delta
      real norm_grad_phi
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     compute one_over_epsilon
      one_over_epsilon = 1.d0/epsilon

c     initialize int_F to zero
      int_F = 0.0d0
 
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                phi_cur = phi(i,j,k)
    
                if (abs(phi_cur) .lt. epsilon) then
                  delta = 0.5d0*one_over_epsilon
     &                  * ( 1.d0+cos(pi*phi_cur*one_over_epsilon) ) 

                  norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

                  int_F = int_F + delta*norm_grad_phi*F(i,j,k)*dV
                endif
       
            enddo
          enddo
        enddo
c       } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm3dMaxNormDiffControlVolume(
     &  max_norm_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  klo_field1_gb, khi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  klo_field2_gb, khi_field2_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer klo_field1_gb, khi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer klo_field2_gb, khi_field2_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb,
     &            klo_field1_gb:khi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb,
     &            klo_field2_gb:khi_field2_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      real max_norm_diff
      real next_diff
      integer i,j,k


c     initialize max_norm_diff
      max_norm_diff = abs( field1(ilo_ib,jlo_ib,klo_ib) 
     &                   - field2(ilo_ib,jlo_ib,klo_ib))

      if (control_vol_sgn .gt. 0) then   
c       loop over included cells { 
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in max norm calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then

                next_diff = abs(field1(i,j,k) - field2(i,j,k))
                if (next_diff .gt. max_norm_diff) then
                  max_norm_diff = next_diff
                endif

              endif
  
            enddo
          enddo
        enddo
c       } end loop over grid 

      else
c       loop over included cells { 
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in max norm calculation if it has a 
c             negative control volume
              if (control_vol(i,j,k) .lt. 0.d0) then

                next_diff = abs(field1(i,j,k) - field2(i,j,k))
                if (next_diff .gt. max_norm_diff) then
                  max_norm_diff = next_diff
                endif

              endif
  
            enddo
          enddo
        enddo
c       } end loop over grid 

      endif      
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dComputeStableAdvectionDtControlVolume(
     &  dt,
     &  vel_x, vel_y, vel_z,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  klo_vel_gb, khi_vel_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer klo_vel_gb, khi_vel_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb,
     &           klo_vel_gb:khi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb,
     &           klo_vel_gb:khi_vel_gb)
      real vel_z(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb,
     &           klo_vel_gb:khi_vel_gb)
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      real dx, dy, dz
      real inv_dx, inv_dy, inv_dz
      real cfl_number
      integer i, j, k
      real max_U_over_dX
      real U_over_dX_cur
      real small_number
      parameter (small_number = 1.d-99)

c     initialize max_U_over_dX to -1
      max_U_over_dX = -1.0d0

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy
      inv_dz = 1.d0/dz

      if (control_vol_sgn .gt. 0) then    
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
c             only include cell in dt calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then 

                U_over_dX_cur = abs(vel_x(i,j,k))*inv_dx
     &                        + abs(vel_y(i,j,k))*inv_dy
     &                        + abs(vel_z(i,j,k))*inv_dz

                if (U_over_dX_cur .gt. max_U_over_dX) then
                  max_U_over_dX = U_over_dX_cur  
                endif
		
              endif
            enddo
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
c             only include cell in dt calculation if it has a 
c             negative control volume
              if (control_vol(i,j,k) .lt. 0.d0) then 

                U_over_dX_cur = abs(vel_x(i,j,k))*inv_dx
     &                        + abs(vel_y(i,j,k))*inv_dy
     &                        + abs(vel_z(i,j,k))*inv_dz

                if (U_over_dX_cur .gt. max_U_over_dX) then
                  max_U_over_dX = U_over_dX_cur  
                endif
		
              endif
            enddo
          enddo
        enddo
c       } end loop over grid

      endif
c     set dt
      dt = cfl_number / (max_U_over_dX + small_number);

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dComputeStableNormalVelDtControlVolume(
     &  dt,
     &  vel_n,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  klo_vel_gb, khi_vel_gb,
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz, 
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer klo_vel_gb, khi_vel_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real vel_n(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb,
     &           klo_vel_gb:khi_vel_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_z_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_z_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy,dz
      real inv_dx, inv_dy, inv_dz
      real max_dx_sq
      real cfl_number
      integer i,j,k
      real max_H_over_dX
      real H_over_dX_cur
      real phi_x_cur, phi_y_cur, phi_z_cur
      real norm_grad_phi
      real small_number
      parameter (small_number = 1.d-99)

c     compute max_dx_sq
      max_dx_sq = max(dx,dy,dz)
      max_dx_sq = max(dx,dy,dz) * max(dx,dy,dz)

c     initialize max_H_over_dX to -1
      max_H_over_dX = -1.0d0

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy
      inv_dz = 1.d0/dz
      
      if (control_vol_sgn .gt. 0) then    
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in dt calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then 
                  phi_x_cur = max(abs(phi_x_plus(i,j,k)),
     &                          abs(phi_x_minus(i,j,k)))
                  phi_y_cur = max(abs(phi_y_plus(i,j,k)),
     &                          abs(phi_y_minus(i,j,k)))
                  phi_z_cur = max(abs(phi_z_plus(i,j,k)),
     &                          abs(phi_z_minus(i,j,k)))
                  norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                              + phi_y_cur*phi_y_cur 
     &                              + phi_z_cur*phi_z_cur + max_dx_sq )

                  H_over_dX_cur = abs(vel_n(i,j,k)) / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy 
     &                        + phi_z_cur*inv_dz )

                if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
                endif
              endif
	      
            enddo
          enddo
        enddo
c       } end loop over grid
      
      else
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in dt calculation if it has a 
c             negative control volume
              if (control_vol(i,j,k) .lt. 0.d0) then 
                  phi_x_cur = max(abs(phi_x_plus(i,j,k)),
     &                          abs(phi_x_minus(i,j,k)))
                  phi_y_cur = max(abs(phi_y_plus(i,j,k)),
     &                          abs(phi_y_minus(i,j,k)))
                  phi_z_cur = max(abs(phi_z_plus(i,j,k)),
     &                          abs(phi_z_minus(i,j,k)))
                  norm_grad_phi = sqrt( phi_x_cur*phi_x_cur 
     &                              + phi_y_cur*phi_y_cur 
     &                              + phi_z_cur*phi_z_cur + max_dx_sq )

                  H_over_dX_cur = abs(vel_n(i,j,k)) / norm_grad_phi
     &                        * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy 
     &                        + phi_z_cur*inv_dz )

                if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
                endif
              endif
	      
            enddo
          enddo
        enddo
c       } end loop over grid
      endif
      
c     set dt
      dt = cfl_number / (max_H_over_dX + small_number);

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dComputeStableConstNormalVelDtControlVolume(
     &  dt,
     &  vel_n,
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real vel_n
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_z_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb,
     &                klo_grad_phi_plus_gb:khi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      real phi_z_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb,
     &                 klo_grad_phi_minus_gb:khi_grad_phi_minus_gb)
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy,dz
      real inv_dx,inv_dy,inv_dz
      real max_dx_sq
      real cfl_number
      integer i,j,k
      real max_H_over_dX
      real H_over_dX_cur
      real phi_x_cur, phi_y_cur, phi_z_cur
      real small_number
      parameter (small_number = 1.d-99)

c     compute max_dx_sq
      max_dx_sq = max(dx,dy,dz)
      max_dx_sq = max_dx_sq * max_dx_sq

c     initialize max_H_over_dX to -1
      max_H_over_dX = -1

c     compute inv_dx, inv_dy, inv_dz
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy
      inv_dz = 1.d0/dz

      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in dt calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then  
                phi_x_cur = max(abs(phi_x_plus(i,j,k)),
     &         	              abs(phi_x_minus(i,j,k)))
                phi_y_cur = max(abs(phi_y_plus(i,j,k)),
     &	                        abs(phi_y_minus(i,j,k)))
                phi_z_cur = max(abs(phi_z_plus(i,j,k)),
     &	                        abs(phi_z_minus(i,j,k)))

                H_over_dX_cur = abs(vel_n) 
     &                    / sqrt( phi_x_cur*phi_x_cur 
     &                          + phi_y_cur*phi_y_cur 
     &                          + phi_z_cur*phi_z_cur + max_dx_sq )
     &                    * ( phi_x_cur*inv_dx 
     &                      + phi_y_cur*inv_dy 
     &                      + phi_z_cur*inv_dz  )
  
                if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
                endif

              endif
            enddo
          enddo
        enddo
c       } end loop over grid
      else
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in dt calculation if it has a 
c             negative control volume
              if (control_vol(i,j,k) .lt. 0.d0) then  
                phi_x_cur = max(abs(phi_x_plus(i,j,k)),
     &         	              abs(phi_x_minus(i,j,k)))
                phi_y_cur = max(abs(phi_y_plus(i,j,k)),
     &	                        abs(phi_y_minus(i,j,k)))
                phi_z_cur = max(abs(phi_z_plus(i,j,k)),
     &	                        abs(phi_z_minus(i,j,k)))

                H_over_dX_cur = abs(vel_n) 
     &                    / sqrt( phi_x_cur*phi_x_cur 
     &                          + phi_y_cur*phi_y_cur 
     &                          + phi_z_cur*phi_z_cur + max_dx_sq )
     &                    * ( phi_x_cur*inv_dx 
     &                      + phi_y_cur*inv_dy 
     &                      + phi_z_cur*inv_dz )

                if (H_over_dX_cur .gt. max_H_over_dX) then
                  max_H_over_dX = H_over_dX_cur  
                endif

              endif
            enddo
          enddo
        enddo
c       } end loop over grid
      endif
      
c     set dt
      dt = cfl_number / (max_H_over_dX + small_number);

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dVolumeIntegralPhiLessThanZeroControlVolume(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  klo_F_gb, khi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer klo_F_gb, khi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb,
     &       klo_F_gb:khi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy,dz
      real epsilon
      integer i,j,k
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     initialize int_F to zero
      int_F = 0.0d0

      if (control_vol_sgn .gt. 0) then    
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
  
c             only include cell in integral calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then

                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .lt. -epsilon) then
                  int_F = int_F + F(i,j,k)*dV
                elseif (phi_cur .lt. epsilon) then
                  one_minus_H = 
     &		    0.5d0*(1.d0-phi_cur_over_epsilon
     &                   -one_over_pi*sin(pi*phi_cur_over_epsilon))
                  int_F = int_F + one_minus_H*F(i,j,k)*dV
                endif

              endif
      
            enddo
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
  
c             only include cell in integral calculation if it has a 
c             negative control volume
              if (control_vol(i,j,k) .lt. 0.d0) then

                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .lt. -epsilon) then
                  int_F = int_F + F(i,j,k)*dV
                elseif (phi_cur .lt. epsilon) then
                  one_minus_H = 
     &		    0.5d0*(1.d0-phi_cur_over_epsilon
     &                   -one_over_pi*sin(pi*phi_cur_over_epsilon))
                  int_F = int_F + one_minus_H*F(i,j,k)*dV
                endif

              endif
      
            enddo
          enddo
        enddo
c       } end loop over grid

      endif      
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dVolumeIntegralPhiGreaterThanZeroControlVolume(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  klo_F_gb, khi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer klo_F_gb, khi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb,
     &       klo_F_gb:khi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy,dz
      real epsilon
      integer i,j,k
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     initialize int_F to zero
      int_F = 0.0d0

      if (control_vol_sgn .gt. 0) then
      
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
  
c             only include cell in integral calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then

                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .gt. epsilon) then
                  int_F = int_F + F(i,j,k)*dV
                elseif (phi_cur .gt. -epsilon) then
                  H = 0.5d0*( 1.d0 + phi_cur_over_epsilon 
     &                             + one_over_pi
     &                             * sin(pi*phi_cur_over_epsilon) )
                  int_F = int_F + H*F(i,j,k)*dV
                endif

              endif
          
            enddo
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
c             only include cell in integral calculation if it has a 
c             negative control volume
              if (control_vol(i,j,k) .lt. 0.d0) then

                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .gt. epsilon) then
                  int_F = int_F + F(i,j,k)*dV
                elseif (phi_cur .gt. -epsilon) then
                  H = 0.5d0*( 1.d0 + phi_cur_over_epsilon 
     &                             + one_over_pi
     &                             * sin(pi*phi_cur_over_epsilon) )
                  int_F = int_F + H*F(i,j,k)*dV
                endif

              endif
          
            enddo
          enddo
        enddo
c       } end loop over grid
      endif

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dSurfaceIntegralControlVolume(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  klo_F_gb, khi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  phi_x, phi_y, phi_z,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  klo_grad_phi_gb, khi_grad_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib,
     &  dx, dy, dz,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer klo_F_gb, khi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb,
     &       klo_F_gb:khi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_z(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy,dz
      real epsilon
      real one_over_epsilon
      integer i,j,k
      real phi_cur
      real delta
      real norm_grad_phi
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     compute one_over_epsilon
      one_over_epsilon = 1.d0/epsilon

c     initialize int_F to zero
      int_F = 0.0d0

      if (control_vol_sgn .gt. 0) then
   
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in integral calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then

                phi_cur = phi(i,j,k)
    
                if (abs(phi_cur) .lt. epsilon) then
                  delta = 0.5d0*one_over_epsilon
     &                  * ( 1.d0+cos(pi*phi_cur*one_over_epsilon) ) 

                  norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

                  int_F = int_F + delta*norm_grad_phi*F(i,j,k)*dV
                endif

              endif
        
            enddo
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in integral calculation if it has a 
c             negative control volume
              if (control_vol(i,j,k) .lt. 0.d0) then

                phi_cur = phi(i,j,k)
    
                if (abs(phi_cur) .lt. epsilon) then
                  delta = 0.5d0*one_over_epsilon
     &                  * ( 1.d0+cos(pi*phi_cur*one_over_epsilon) ) 

                  norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

                  int_F = int_F + delta*norm_grad_phi*F(i,j,k)*dV
                endif

              endif
        
            enddo
          enddo
        enddo
c       } end loop over grid
      endif

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm3dVoxelCountGreaterThanZeroControlVolume(
     &  count,
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      integer count
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      
c     local variables      
      integer i,j,k
      
      count = 0
      if (control_vol_sgn .gt. 0) then      
c       loop over grid { 
        do k=klo_ib,khi_ib
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j,k) .gt. 0.d0 .and.
     &	        control_vol(i,j,k) .ge. 0.d0) then
	           count = count + 1
	     endif
	        
          enddo
	 enddo 
        enddo	
c       } end loop over grid     
      else
        
c       loop over grid { 
        do k=klo_ib,khi_ib
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j,k) .gt. 0.d0 .and.
     &	        control_vol(i,j,k) .le. 0.d0) then
	           count = count + 1
	     endif
	        
          enddo
	 enddo 
        enddo	
c       } end loop over grid     
      
      endif
      	 
      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm3dVoxelCountLessThanZeroControlVolume(
     &  count,
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  klo_gb, khi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  klo_control_vol_gb, khi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  klo_ib, khi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      integer count
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer klo_gb, khi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb,klo_gb:khi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb,
     &                 klo_control_vol_gb:khi_control_vol_gb)
      integer control_vol_sgn
      
c     local variables      
      integer i,j,k
      
      count = 0
      if (control_vol_sgn .gt. 0) then      
c       loop over grid { 
        do k=klo_ib,khi_ib
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j,k) .lt. 0.d0 .and.
     &	        control_vol(i,j,k) .ge. 0.d0) then
	           count = count + 1
	     endif
	        
          enddo
	 enddo 
        enddo	
c       } end loop over grid     
      else
        
c       loop over grid { 
        do k=klo_ib,khi_ib
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j,k) .lt. 0.d0 .and.
     &	        control_vol(i,j,k) .le. 0.d0) then
	           count = count + 1
	     endif
	        
          enddo
	 enddo 
        enddo	
c       } end loop over grid     
      
      endif
      	 
      return
      end
c } end subroutine
c***********************************************************************