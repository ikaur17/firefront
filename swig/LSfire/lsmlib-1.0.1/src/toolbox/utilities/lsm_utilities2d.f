c***********************************************************************
c
c  File:        lsm_utilities2d.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 routines for 2D level set method utility subroutines
c
c***********************************************************************


c***********************************************************************
      subroutine lsm2dMaxNormDiff(
     &  max_norm_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb)
      real max_norm_diff
      real next_diff
      integer i,j

c     initialize max_norm_diff
      max_norm_diff = abs( field1(ilo_ib,jlo_ib) 
     &                   - field2(ilo_ib,jlo_ib))

c       loop over included cells { 
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

              next_diff = abs(field1(i,j) - field2(i,j))
              if (next_diff .gt. max_norm_diff) then
                max_norm_diff = next_diff
              endif

          enddo
        enddo
c       } end loop over grid 
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dAveAbsDiff(
     &  ave_abs_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      real ave_abs_diff
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &                        jlo_field1_gb:jhi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &                        jlo_field2_gb:jhi_field2_gb)
      
c     local variables      
      real sum_abs_diff, num_pts, next_diff
      real zero, one
      parameter (zero=0.d0, one=1.d0)
      integer i,j

c     initialize max_norm_diff
      sum_abs_diff = zero
      num_pts = zero

c       loop over grid { 
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	      next_diff = abs(field1(i,j) - field2(i,j))
              sum_abs_diff = sum_abs_diff + next_diff
              num_pts = num_pts + one
	      
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
      subroutine lsm2dVoxelCountGreaterThanZero(
     &  count,
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      integer count
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
c     local variables      
      integer i,j

      count = 0
c       loop over grid { 
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j) .gt. 0) then
	       count = count + 1
	     endif
	        
          enddo
        enddo
c       } end loop over grid     
      	 
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dVoxelCountLessThanZero(
     &  count,
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      integer count
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      
c     local variables      
      integer i,j

      count = 0
c       loop over grid { 
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j) .lt. 0.d0) then
	       count = count + 1
	     endif
	        
          enddo
        enddo
c       } end loop over grid     
      	 
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dComputeStableAdvectionDt(
     &  dt,
     &  vel_x, vel_y,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real dx, dy
      real inv_dx, inv_dy
      real cfl_number
      integer i, j
      real max_U_over_dX
      real U_over_dX_cur
      real small_number
      parameter (small_number = 1.d-99)

c     initialize max_U_over_dX to -1
      max_U_over_dX = -1.0d0

c     compute inv_dx and inv_dx
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy

c       loop over included cells {
    
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
  
              U_over_dX_cur = abs(vel_x(i,j))*inv_dx
     &                    + abs(vel_y(i,j))*inv_dy

              if (U_over_dX_cur .gt. max_U_over_dX) then
                max_U_over_dX = U_over_dX_cur  
              endif
  
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
      subroutine lsm2dComputeStableNormalVelDt(
     &  dt,
     &  vel_n,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real vel_n(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real dx,dy
      real inv_dx,inv_dy
      real max_dx_sq
      real cfl_number
      integer i,j
      real max_H_over_dX
      real H_over_dX_cur
      real phi_x_cur, phi_y_cur
      real small_number
      parameter (small_number = 1.d-99)

c     compute max_dx_sq
      max_dx_sq = max(dx,dy)
      max_dx_sq = max_dx_sq * max_dx_sq

c     initialize max_H_over_dX to -1
      max_H_over_dX = -1.0d0

c     compute inv_dx and inv_dy
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy

c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

              phi_x_cur = max(abs(phi_x_plus(i,j)),
     & 	                      abs(phi_x_minus(i,j)))
              phi_y_cur = max(abs(phi_y_plus(i,j)),
     &	                      abs(phi_y_minus(i,j)))

              H_over_dX_cur = abs(vel_n(i,j)) 
     &                      / sqrt( phi_x_cur*phi_x_cur 
     &                        + phi_y_cur*phi_y_cur + max_dx_sq )
     &                      * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy )

              if (H_over_dX_cur .gt. max_H_over_dX) then
                max_H_over_dX = H_over_dX_cur  
              endif

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
      subroutine lsm2dComputeStableConstNormalVelDt(
     &  dt,
     &  vel_n,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real vel_n
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real dx,dy
      real inv_dx,inv_dy
      real max_dx_sq
      real cfl_number
      integer i,j
      real max_H_over_dX
      real H_over_dX_cur
      real phi_x_cur, phi_y_cur
      real small_number
      parameter (small_number = 1.d-99)

c     compute max_dx_sq
      max_dx_sq = max(dx,dy)
      max_dx_sq = max_dx_sq * max_dx_sq

c     initialize max_H_over_dX to -1
      max_H_over_dX = -1

c     compute inv_dx and inv_dy
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy

c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

              phi_x_cur = max(abs(phi_x_plus(i,j)),
     &       	              abs(phi_x_minus(i,j)))
              phi_y_cur = max(abs(phi_y_plus(i,j)),
     &	                      abs(phi_y_minus(i,j)))

              H_over_dX_cur = abs(vel_n) 
     &                  / sqrt( phi_x_cur*phi_x_cur 
     &                        + phi_y_cur*phi_y_cur + max_dx_sq )
     &                  * ( phi_x_cur*inv_dx 
     &                    + phi_y_cur*inv_dy  )

              if (H_over_dX_cur .gt. max_H_over_dX) then
                max_H_over_dX = H_over_dX_cur  
              endif

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
      subroutine lsm2dVolumeIntegralPhiLessThanZero(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy
      dV = dx * dy

c     initialize int_F to zero
      int_F = 0.0d0

c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon

              if (phi_cur .lt. -epsilon) then
                int_F = int_F + F(i,j)*dV
              elseif (phi_cur .lt. epsilon) then
                one_minus_H = 0.5d0*(1.d0 - phi_cur_over_epsilon
     &                                  - one_over_pi
     &                                  * sin(pi*phi_cur_over_epsilon))
                int_F = int_F + one_minus_H*F(i,j)*dV
              endif
     
          enddo
        enddo
c       } end loop over grid
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dVolumeIntegralPhiGreaterThanZero(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy
      dV = dx * dy

c     initialize int_F to zero
      int_F = 0.0d0

c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon
  
              if (phi_cur .gt. epsilon) then
                int_F = int_F + F(i,j)*dV
              elseif (phi_cur .gt. -epsilon) then
                H = 0.5d0*( 1.d0 + phi_cur_over_epsilon 
     &                           + one_over_pi
     &                           * sin(pi*phi_cur_over_epsilon) )
                int_F = int_F + H*F(i,j)*dV
              endif
       
          enddo
        enddo
c       } end loop over grid
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dSurfaceIntegral(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real dx,dy
      real epsilon
      real one_over_epsilon
      integer i,j
      real phi_cur
      real delta
      real norm_grad_phi
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      

c     compute dV = dx * dy
      dV = dx * dy

c     compute one_over_epsilon
      one_over_epsilon = 1.d0/epsilon

c     initialize int_F to zero
      int_F = 0.0d0
      
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

              phi_cur = phi(i,j)
  
              if (abs(phi_cur) .lt. epsilon) then
                delta = 0.5d0*one_over_epsilon
     &              * ( 1.d0+cos(pi*phi_cur*one_over_epsilon) ) 

                norm_grad_phi = sqrt(
     &              phi_x(i,j)*phi_x(i,j)
     &            + phi_y(i,j)*phi_y(i,j) )

                int_F = int_F + delta*norm_grad_phi*F(i,j)*dV
              endif
   
          enddo
        enddo
c       } end loop over grid
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dSurfaceIntegralPrecomputedDelta(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  delta_phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  grad_phi_mag,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb)
      real delta_phi(ilo_phi_gb:ihi_phi_gb,
     &               jlo_phi_gb:jhi_phi_gb)
      real grad_phi_mag(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                  jlo_grad_phi_gb:jhi_grad_phi_gb)     
       
      real dx,dy
      integer i,j
      real dV
      

c     compute dV = dx * dy
      dV = dx * dy

c     initialize int_F to zero
      int_F = 0.0d0
      
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
  
              int_F = int_F + delta_phi(i,j)*grad_phi_mag(i,j)*F(i,j)*dV
          enddo
        enddo
c       } end loop over grid
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dMaxNormDiffControlVolume(
     &  max_norm_diff,
     &  field1,
     &  ilo_field1_gb, ihi_field1_gb,
     &  jlo_field1_gb, jhi_field1_gb,
     &  field2,
     &  ilo_field2_gb, ihi_field2_gb,
     &  jlo_field2_gb, jhi_field2_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,     
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_field1_gb, ihi_field1_gb
      integer jlo_field1_gb, jhi_field1_gb
      integer ilo_field2_gb, ihi_field2_gb
      integer jlo_field2_gb, jhi_field2_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real field1(ilo_field1_gb:ihi_field1_gb,
     &            jlo_field1_gb:jhi_field1_gb)
      real field2(ilo_field2_gb:ihi_field2_gb,
     &            jlo_field2_gb:jhi_field2_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real max_norm_diff
      real next_diff
      integer i,j

c     initialize max_norm_diff
      max_norm_diff = abs( field1(ilo_ib,jlo_ib) 
     &                   - field2(ilo_ib,jlo_ib))

      if (control_vol_sgn .gt. 0) then   
c       loop over included cells { 
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in max norm calculation if it has a 
c           positive control volume
            if (control_vol(i,j) .gt. 0.d0) then

              next_diff = abs(field1(i,j) - field2(i,j))
              if (next_diff .gt. max_norm_diff) then
                max_norm_diff = next_diff
              endif

            endif

          enddo
        enddo
c       } end loop over grid 

      else
c       loop over included cells { 
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in max norm calculation if it has a 
c           negative control volume
            if (control_vol(i,j) .lt. 0.d0) then

              next_diff = abs(field1(i,j) - field2(i,j))
              if (next_diff .gt. max_norm_diff) then
                max_norm_diff = next_diff
              endif

            endif

          enddo
        enddo
c       } end loop over grid 
      endif
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dComputeStableAdvectionDtControlVolume(
     &  dt,
     &  vel_x, vel_y,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx, dy
      real inv_dx, inv_dy
      real cfl_number
      integer i, j
      real max_U_over_dX
      real U_over_dX_cur
      real small_number
      parameter (small_number = 1.d-99)

c     initialize max_U_over_dX to -1
      max_U_over_dX = -1.0d0

c     compute inv_dx and inv_dx
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy
      
      if (control_vol_sgn .gt. 0) then   
c       loop over included cells {
    
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
c           only include cell in dt calculation if it has a 
c           positive control volume
            if (control_vol(i,j) .gt. 0.d0) then  
              U_over_dX_cur = abs(vel_x(i,j))*inv_dx
     &                    + abs(vel_y(i,j))*inv_dy

              if (U_over_dX_cur .gt. max_U_over_dX) then
                max_U_over_dX = U_over_dX_cur  
              endif
            endif	      
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {   
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
c           only include cell in dt calculation if it has a 
c           negative control volume
            if (control_vol(i,j) .lt. 0.d0) then  
              U_over_dX_cur = abs(vel_x(i,j))*inv_dx
     &                    + abs(vel_y(i,j))*inv_dy

              if (U_over_dX_cur .gt. max_U_over_dX) then
                max_U_over_dX = U_over_dX_cur  
              endif
            endif	      
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
      subroutine lsm2dComputeStableNormalVelDtControlVolume(
     &  dt,
     &  vel_n,
     &  ilo_vel_gb, ihi_vel_gb,
     &  jlo_vel_gb, jhi_vel_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_vel_gb, ihi_vel_gb
      integer jlo_vel_gb, jhi_vel_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real vel_n(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy
      real inv_dx,inv_dy
      real max_dx_sq
      real cfl_number
      integer i,j
      real max_H_over_dX
      real H_over_dX_cur
      real phi_x_cur, phi_y_cur
      real small_number
      parameter (small_number = 1.d-99)

c     compute max_dx_sq
      max_dx_sq = max(dx,dy)
      max_dx_sq = max_dx_sq * max_dx_sq

c     initialize max_H_over_dX to -1
      max_H_over_dX = -1.0d0

c     compute inv_dx and inv_dy
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy

      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
c           only include cell in dt calculation if it has a 
c           positive control volume
            if (control_vol(i,j) .gt. 0.d0) then  
              phi_x_cur = max(abs(phi_x_plus(i,j)),
     & 	                      abs(phi_x_minus(i,j)))
              phi_y_cur = max(abs(phi_y_plus(i,j)),
     &	                      abs(phi_y_minus(i,j)))

              H_over_dX_cur = abs(vel_n(i,j)) 
     &                      / sqrt( phi_x_cur*phi_x_cur 
     &                        + phi_y_cur*phi_y_cur + max_dx_sq )
     &                      * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy  )

              if (H_over_dX_cur .gt. max_H_over_dX) then
                max_H_over_dX = H_over_dX_cur  
              endif
            endif
          enddo
        enddo
c       } end loop over grid

      else      
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib
c           only include cell in dt calculation if it has a 
c           negative control volume
            if (control_vol(i,j) .lt. 0.d0) then  
              phi_x_cur = max(abs(phi_x_plus(i,j)),
     & 	                      abs(phi_x_minus(i,j)))
              phi_y_cur = max(abs(phi_y_plus(i,j)),
     &	                      abs(phi_y_minus(i,j)))
              H_over_dX_cur = abs(vel_n(i,j)) 
     &                      / sqrt( phi_x_cur*phi_x_cur 
     &                        + phi_y_cur*phi_y_cur + max_dx_sq )
     &                      * ( phi_x_cur*inv_dx 
     &                        + phi_y_cur*inv_dy )

              if (H_over_dX_cur .gt. max_H_over_dX) then
                max_H_over_dX = H_over_dX_cur  
              endif
            endif
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
      subroutine lsm2dComputeStableConstNormalVelDtControlVolume(
     &  dt,
     &  vel_n,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  cfl_number)
c***********************************************************************
c { begin subroutine
      implicit none

      real dt

c     _gb refers to ghostbox 
c     _ib refers to box to include in dt calculation
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real vel_n
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy
      real inv_dx,inv_dy
      real max_dx_sq
      real cfl_number
      integer i,j
      real max_H_over_dX
      real H_over_dX_cur
      real phi_x_cur, phi_y_cur
      real small_number
      parameter (small_number = 1.d-99)

c     compute max_dx_sq
      max_dx_sq = max(dx,dy)
      max_dx_sq = max_dx_sq * max_dx_sq

c     initialize max_H_over_dX to -1
      max_H_over_dX = -1

c     compute inv_dx and inv_dy
      inv_dx = 1.d0/dx
      inv_dy = 1.d0/dy

      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in dt calculation if it has a 
c           positive control volume
            if (control_vol(i,j) .gt. 0.d0) then  
              phi_x_cur = max(abs(phi_x_plus(i,j)),
     &       	              abs(phi_x_minus(i,j)))
              phi_y_cur = max(abs(phi_y_plus(i,j)),
     &	                      abs(phi_y_minus(i,j)))

              H_over_dX_cur = abs(vel_n) 
     &                  / sqrt( phi_x_cur*phi_x_cur 
     &                        + phi_y_cur*phi_y_cur + max_dx_sq )
     &                  * ( phi_x_cur*inv_dx 
     &                    + phi_y_cur*inv_dy  )

              if (H_over_dX_cur .gt. max_H_over_dX) then
                max_H_over_dX = H_over_dX_cur  
              endif

            endif
          enddo
        enddo
c       } end loop over grid
      else
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in dt calculation if it has a 
c           negative control volume
            if (control_vol(i,j) .lt. 0.d0) then  
              phi_x_cur = max(abs(phi_x_plus(i,j)),
     &       	              abs(phi_x_minus(i,j)))
              phi_y_cur = max(abs(phi_y_plus(i,j)),
     &	                      abs(phi_y_minus(i,j)))
              H_over_dX_cur = abs(vel_n) 
     &                  / sqrt( phi_x_cur*phi_x_cur 
     &                        + phi_y_cur*phi_y_cur + max_dx_sq )
     &                  * ( phi_x_cur*inv_dx 
     &                    + phi_y_cur*inv_dy  )

              if (H_over_dX_cur .gt. max_H_over_dX) then
                max_H_over_dX = H_over_dX_cur  
              endif

            endif
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
      subroutine lsm2dVolumeIntegralPhiLessThanZeroControlVolume(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real one_minus_H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy
      dV = dx * dy

c     initialize int_F to zero
      int_F = 0.0d0

      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in integral if it has a positive control volume
            if (control_vol(i,j) .gt. 0.d0) then

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon

              if (phi_cur .lt. -epsilon) then
                int_F = int_F + F(i,j)*dV
              elseif (phi_cur .lt. epsilon) then
                one_minus_H = 0.5d0*(1.d0 - phi_cur_over_epsilon
     &                                  - one_over_pi
     &                                  * sin(pi*phi_cur_over_epsilon))
                int_F = int_F + one_minus_H*F(i,j)*dV
              endif

            endif
      
          enddo
        enddo
c       } end loop over grid
    
      else
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in integral if it has a negative control volume
            if (control_vol(i,j) .lt. 0.d0) then

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon

              if (phi_cur .lt. -epsilon) then
                int_F = int_F + F(i,j)*dV
              elseif (phi_cur .lt. epsilon) then
                one_minus_H = 0.5d0*(1.d0 - phi_cur_over_epsilon
     &                                  - one_over_pi
     &                                  * sin(pi*phi_cur_over_epsilon))
                int_F = int_F + one_minus_H*F(i,j)*dV
              endif

            endif
      
          enddo
        enddo
c       } end loop over grid      

      endif
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dVolumeIntegralPhiGreaterThanZeroControlVolume(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy
      real epsilon
      integer i,j
      real phi_cur
      real phi_cur_over_epsilon
      real H
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      real one_over_pi
      parameter (one_over_pi=0.31830988618379d0)
      

c     compute dV = dx * dy
      dV = dx * dy

c     initialize int_F to zero
      int_F = 0.0d0

      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in integral if it has a positive control volume
            if (control_vol(i,j) .gt. 0.d0) then

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon
  
              if (phi_cur .gt. epsilon) then
                int_F = int_F + F(i,j)*dV
              elseif (phi_cur .gt. -epsilon) then
                H = 0.5d0*( 1.d0 + phi_cur_over_epsilon 
     &                           + one_over_pi
     &                           * sin(pi*phi_cur_over_epsilon) )
                int_F = int_F + H*F(i,j)*dV
              endif

            endif
        
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in integral if it has a negative control volume
            if (control_vol(i,j) .lt. 0.d0) then

              phi_cur = phi(i,j)
              phi_cur_over_epsilon = phi_cur/epsilon
  
              if (phi_cur .gt. epsilon) then
                int_F = int_F + F(i,j)*dV
              elseif (phi_cur .gt. -epsilon) then
                H = 0.5d0*( 1.d0 + phi_cur_over_epsilon 
     &                           + one_over_pi
     &                           * sin(pi*phi_cur_over_epsilon) )
                int_F = int_F + H*F(i,j)*dV
              endif

            endif
        
          enddo
        enddo
c       } end loop over grid

      endif
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dSurfaceIntegralControlVolume(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy,
     &  epsilon)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      real dx,dy
      real epsilon
      real one_over_epsilon
      integer i,j
      real phi_cur
      real delta
      real norm_grad_phi
      real dV
      real pi
      parameter (pi=3.14159265358979323846d0)
      

c     compute dV = dx * dy
      dV = dx * dy

c     compute one_over_epsilon
      one_over_epsilon = 1.d0/epsilon

c     initialize int_F to zero
      int_F = 0.0d0
      
      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in integral if it has a positive control volume
            if (control_vol(i,j) .gt. 0.d0) then

              phi_cur = phi(i,j)
  
              if (abs(phi_cur) .lt. epsilon) then
                delta = 0.5d0*one_over_epsilon
     &              * ( 1.d0+cos(pi*phi_cur*one_over_epsilon) ) 

                norm_grad_phi = sqrt(
     &              phi_x(i,j)*phi_x(i,j)
     &            + phi_y(i,j)*phi_y(i,j) )

                int_F = int_F + delta*norm_grad_phi*F(i,j)*dV
              endif

            endif
      
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in integral if it has a positive control volume
            if (control_vol(i,j) .lt. 0.d0) then

              phi_cur = phi(i,j)
  
              if (abs(phi_cur) .lt. epsilon) then
                delta = 0.5d0*one_over_epsilon
     &              * ( 1.d0+cos(pi*phi_cur*one_over_epsilon) ) 

                norm_grad_phi = sqrt(
     &              phi_x(i,j)*phi_x(i,j)
     &            + phi_y(i,j)*phi_y(i,j) )

                int_F = int_F + delta*norm_grad_phi*F(i,j)*dV
              endif

            endif
      
          enddo
        enddo
c       } end loop over grid      

      endif
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dVoxelCountGreaterThanZeroControlVolume(
     &  count,
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      integer count
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      
c     local variables      
      integer i,j
      
      count = 0
      if (control_vol_sgn .gt. 0) then      
c       loop over grid { 
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j) .gt. 0.d0 .and.
     &	        control_vol(i,j) .ge. 0.d0) then
	           count = count + 1
	     endif
	        
	 enddo 
        enddo	
c       } end loop over grid     
      else
        
c       loop over grid { 
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j) .gt. 0.d0 .and.
     &	        control_vol(i,j) .le. 0.d0) then
	           count = count + 1
	     endif
	        
          enddo
	 enddo
c       } end loop over grid     
      
      endif
      	 
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dVoxelCountLessThanZeroControlVolume(
     &  count,
     &  phi,
     &  ilo_gb, ihi_gb,
     &  jlo_gb, jhi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib)
c***********************************************************************
c { begin subroutine
      implicit none

      integer count
c     _gb refers to ghostbox 
c     _ib refers to box to include in norm calculation
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real phi(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      
c     local variables      
      integer i,j
      
      count = 0
      if (control_vol_sgn .gt. 0) then      
c       loop over grid { 
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j) .lt. 0.d0 .and.
     &	        control_vol(i,j) .ge. 0.d0) then
	           count = count + 1
	     endif
	        
	 enddo 
        enddo	
c       } end loop over grid     
      else
        
c       loop over grid { 
	 do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

	     if(phi(i,j) .lt. 0.d0 .and.
     &	        control_vol(i,j) .le. 0.d0) then
	           count = count + 1
	     endif
	        
          enddo
	 enddo
c       } end loop over grid     
      
      endif
      	 
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dSurfaceIntegralPrecomputedDeltaControlVolume(
     &  int_F,
     &  F,
     &  ilo_F_gb, ihi_F_gb,
     &  jlo_F_gb, jhi_F_gb,
     &  delta_phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  grad_phi_mag,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  control_vol,
     &  ilo_control_vol_gb, ihi_control_vol_gb,
     &  jlo_control_vol_gb, jhi_control_vol_gb,
     &  control_vol_sgn,
     &  ilo_ib, ihi_ib,
     &  jlo_ib, jhi_ib,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

      real int_F

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_F_gb, ihi_F_gb
      integer jlo_F_gb, jhi_F_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      real F(ilo_F_gb:ihi_F_gb,
     &       jlo_F_gb:jhi_F_gb)
      real delta_phi(ilo_phi_gb:ihi_phi_gb,
     &               jlo_phi_gb:jhi_phi_gb)
      real grad_phi_mag(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                  jlo_grad_phi_gb:jhi_grad_phi_gb)     
      real control_vol(ilo_control_vol_gb:ihi_control_vol_gb,
     &                 jlo_control_vol_gb:jhi_control_vol_gb)
      integer control_vol_sgn
      
      real dx,dy
      integer i,j
      real dV
      

c     compute dV = dx * dy
      dV = dx * dy

c     initialize int_F to zero
      int_F = 0.0d0
      
      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in integral if it has a positive control volume
            if (control_vol(i,j) .gt. 0.d0) then
  
              int_F = int_F + delta_phi(i,j)*grad_phi_mag(i,j)*F(i,j)*dV

            endif
      
          enddo
        enddo
c       } end loop over grid

      else
c       loop over included cells {
        do j=jlo_ib,jhi_ib
          do i=ilo_ib,ihi_ib

c           only include cell in integral if it has a negative control volume
            if (control_vol(i,j) .lt. 0.d0) then

              int_F = int_F + delta_phi(i,j)*grad_phi_mag(i,j)*F(i,j)*dV
      
            endif
      
          enddo
        enddo
c       } end loop over grid      

      endif
      
      return
      end
c } end subroutine
c***********************************************************************