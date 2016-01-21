c***********************************************************************
c
c  File:        lsm_geometry3d_fort.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 routines for several common level set method
c               geometry calculations
c
c***********************************************************************

c***********************************************************************
c
c  lsm3dComputeUnitNormal() computes the unit normal vector to the 
c  interface from grad(phi). 
c
c  Arguments:
c    normal_* (out):   components of unit normal vector
c    phi_* (in):       components of grad(phi) 
c    dx, dy, dz (in):  grid spacing
c    *_gb (in):        index range for ghostbox
c    *_fb (in):        index range for fillbox
c
c***********************************************************************
      subroutine lsm3dComputeUnitNormal(
     &  normal_x, normal_y, normal_z,
     &  ilo_normal_gb, ihi_normal_gb,
     &  jlo_normal_gb, jhi_normal_gb,
     &  klo_normal_gb, khi_normal_gb,
     &  phi_x, phi_y, phi_z,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  klo_grad_phi_gb, khi_grad_phi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  klo_fb, khi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostboxes 
c     _fb refers to fill-box for normal data

      integer ilo_normal_gb, ihi_normal_gb
      integer jlo_normal_gb, jhi_normal_gb
      integer klo_normal_gb, khi_normal_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      integer klo_fb, khi_fb
      real normal_x(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb,
     &              klo_normal_gb:khi_normal_gb)
      real normal_y(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb,
     &              klo_normal_gb:khi_normal_gb)
      real normal_z(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb,
     &              klo_normal_gb:khi_normal_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_z(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real norm_grad_phi, inv_norm_grad_phi
      integer i,j,k
      real half
      parameter (half=0.5d0)
      real zero_tol
      parameter (zero_tol=1.d-11)

c     { begin loop over grid
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

c           compute unit normal 

            norm_grad_phi = sqrt( phi_x(i,j,k)*phi_x(i,j,k)
     &                          + phi_y(i,j,k)*phi_y(i,j,k)
     &                          + phi_z(i,j,k)*phi_z(i,j,k) )

            if (norm_grad_phi .ge. zero_tol) then
              inv_norm_grad_phi = 1.0d0/norm_grad_phi
              normal_x(i,j,k) = phi_x(i,j,k)*inv_norm_grad_phi
              normal_y(i,j,k) = phi_y(i,j,k)*inv_norm_grad_phi
              normal_z(i,j,k) = phi_z(i,j,k)*inv_norm_grad_phi
            else
              normal_x(i,j,k) = 1.0d0
              normal_y(i,j,k) = 0.0d0
              normal_z(i,j,k) = 0.0d0
            endif

          enddo
        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dComputeSignedUnitNormal() computes the signed unit normal 
c  vector (sgn(phi)*normal) to the interface from grad(phi) using 
c  the following smoothed sgn function 
c
c    sgn(phi) = phi/sqrt( phi^2 + |grad(phi)|^2 * dx^2 )
c
c  Arguments:
c    normal_* (out):     components of unit normal vector
c    phi_* (in):         components of grad(phi) 
c    phi (in):           level set function
c    dx, dy, dz (in):    grid spacing
c    *_gb (in):          index range for ghostbox
c    *_fb (in):          index range for fillbox
c
c***********************************************************************
      subroutine lsm3dComputeSignedUnitNormal(
     &  normal_x, normal_y, normal_z,
     &  ilo_normal_gb, ihi_normal_gb,
     &  jlo_normal_gb, jhi_normal_gb,
     &  klo_normal_gb, khi_normal_gb,
     &  phi_x, phi_y, phi_z,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  klo_grad_phi_gb, khi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  klo_fb, khi_fb,
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostboxes 
c     _fb refers to fill-box for normal data

      integer ilo_normal_gb, ihi_normal_gb
      integer jlo_normal_gb, jhi_normal_gb
      integer klo_normal_gb, khi_normal_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      integer klo_fb, khi_fb
      real normal_x(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb,
     &              klo_normal_gb:khi_normal_gb)
      real normal_y(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb,
     &              klo_normal_gb:khi_normal_gb)
      real normal_z(ilo_normal_gb:ihi_normal_gb,
     &              jlo_normal_gb:jhi_normal_gb,
     &              klo_normal_gb:khi_normal_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_z(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb,
     &           klo_grad_phi_gb:khi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real dx, dy, dz
      real phi_cur
      real sgn_phi
      real norm_grad_phi_sq, inv_norm_grad_phi
      real dx_sq
      integer i,j,k
      real zero_tol
      parameter (zero_tol=1.d-11)

c     set value of dx_sq to be square of max{dx,dy,dz}
      dx_sq = max(dx,dy,dz)
      dx_sq = dx_sq*dx_sq

c     { begin loop over grid
      do k=klo_fb,khi_fb
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb

c           cache phi_cur
            phi_cur = phi(i,j,k)

c           compute sgn(phi)*normal
            if (abs(phi_cur) .gt. zero_tol) then
              norm_grad_phi_sq = phi_x(i,j,k)*phi_x(i,j,k)
     &                         + phi_y(i,j,k)*phi_y(i,j,k)
     &                         + phi_z(i,j,k)*phi_z(i,j,k)

              if (norm_grad_phi_sq .ge. zero_tol) then

                sgn_phi = phi_cur
     &                  / sqrt(phi_cur*phi_cur + norm_grad_phi_sq*dx_sq)

                inv_norm_grad_phi = 1.d0/sqrt(norm_grad_phi_sq)

                normal_x(i,j,k) = sgn_phi*phi_x(i,j,k)*inv_norm_grad_phi
                normal_y(i,j,k) = sgn_phi*phi_y(i,j,k)*inv_norm_grad_phi
                normal_z(i,j,k) = sgn_phi*phi_z(i,j,k)*inv_norm_grad_phi

              else
                normal_x(i,j,k) = 1.0d0
                normal_y(i,j,k) = 0.0d0
                normal_z(i,j,k) = 0.0d0
              endif

            else

              normal_x(i,j,k) = 0.0d0
              normal_y(i,j,k) = 0.0d0
              normal_z(i,j,k) = 0.0d0

            endif

          enddo
        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dVolumeRegionPhiLessThanZero() computes the volume of the 
c  region where the level set function is less than 0.  
c
c  Arguments:
c    volume (out):          volume of the region where phi < 0
c    phi (in):              level set function
c    dx, dy, dz (in):       grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm3dVolumeRegionPhiLessThanZero(
     &  volume,
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

      real volume

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
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

c     initialize volume to zero
      volume = 0.0d0
           
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
  
                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
    
                if (phi_cur .lt. -epsilon) then
                  volume = volume + dV
                elseif (phi_cur .lt. epsilon) then
                  one_minus_H = 0.5d0*(1 - phi_cur_over_epsilon
     &                                   - one_over_pi
     &                                   * sin(pi*phi_cur_over_epsilon))
                  volume = volume + one_minus_H*dV
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
c
c  lsm3dVolumeRegionPhiGreaterThanZero() computes the volume of the 
c  region where the level set function is greater than 0.  
c
c  Arguments:
c    volume (out):          volume of the region where phi > 0
c    phi (in):              level set function
c    dx, dy, dz (in):       grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm3dVolumeRegionPhiGreaterThanZero(
     &  volume,
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

      real volume

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
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

c     initialize volume to zero
      volume = 0.0d0

c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .gt. epsilon) then
                  volume = volume + dV
                elseif (phi_cur .gt. -epsilon) then
                  H = 0.5*(1 + phi_cur_over_epsilon 
     &                       + one_over_pi*sin(pi*phi_cur_over_epsilon))
                  volume = volume + H*dV
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
c
c  lsm3dSurfaceAreaZeroLevelSet() computes the surface area of the 
c  surface defined by the zero level set. 
c
c  Arguments:
c    area (out):            area of the surface defined by the zero level 
c                           set
c    phi (in):              level set function
c    phi_* (in):            components of grad(phi)
c    dx, dy, dz (in):       grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm3dSurfaceAreaZeroLevelSet(
     &  area,
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

      real area

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
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

c     initialize area to zero
      area = 0.0d0

c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

                phi_cur = phi(i,j,k)
    
                if (abs(phi_cur) .lt. epsilon) then
                  delta = 0.5d0*one_over_epsilon
     &                         *( 1+cos(pi*phi_cur*one_over_epsilon) ) 

                  norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

                  area = area + delta*norm_grad_phi*dV
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
      subroutine lsm3dSurfaceAreaZeroLevelSetDelta(
     &  area,
     &  delta_phi,
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
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none

      real area

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
      real delta_phi(ilo_phi_gb:ihi_phi_gb,
     &               jlo_phi_gb:jhi_phi_gb,
     &               klo_phi_gb:khi_phi_gb)
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
      integer i,j,k
      real norm_grad_phi
      real dV
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     initialize area to zero
      area = 0.0d0

c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

             if( delta_phi(i,j,k) .gt. 0.d0) then
              norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

              area = area + delta_phi(i,j,k)*norm_grad_phi*dV
      
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
c
c  lsm3dVolumeRegionPhiLessThanZeroControlVolume() computes the volume 
c  of the region of the computational domain where the level set
c  function is less than 0.  The computational domain contains only
c  those cells that are included by the control volume data.
c
c  Arguments:
c    volume (out):          volume of the region where phi < 0
c    phi (in):              level set function
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx, dy, dz (in):       grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm3dVolumeRegionPhiLessThanZeroControlVolume(
     &  volume,
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

      real volume

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
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

c     initialize volume to zero
      volume = 0.0d0

      if (control_vol_sgn .gt. 0) then
      
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
   
c             only include cell in max norm calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then
  
                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
    
                if (phi_cur .lt. -epsilon) then
                  volume = volume + dV
                elseif (phi_cur .lt. epsilon) then
                  one_minus_H = 0.5d0*(1 - phi_cur_over_epsilon
     &                                   - one_over_pi
     &                                   * sin(pi*phi_cur_over_epsilon))
                  volume = volume + one_minus_H*dV
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
  
                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
    
                if (phi_cur .lt. -epsilon) then
                  volume = volume + dV
                elseif (phi_cur .lt. epsilon) then
                  one_minus_H = 0.5d0*(1 - phi_cur_over_epsilon
     &                                   - one_over_pi
     &                                   * sin(pi*phi_cur_over_epsilon))
                  volume = volume + one_minus_H*dV
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
c
c  lsm3dVolumeRegionPhiGreaterThanZeroControlVolume() computes the 
c  volume of the region of the computational domain where the level 
c  set function is greater than 0.  The computational domain contains 
c  only those cells that are included by the control volume data.
c
c  Arguments:
c    volume (out):          volume of the region where phi > 0
c    phi (in):              level set function
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx, dy, dz (in):       grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm3dVolumeRegionPhiGreaterThanZeroControlVolume(
     &  volume,
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

      real volume

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_control_vol_gb, ihi_control_vol_gb
      integer jlo_control_vol_gb, jhi_control_vol_gb
      integer klo_control_vol_gb, khi_control_vol_gb
      integer ilo_ib, ihi_ib
      integer jlo_ib, jhi_ib
      integer klo_ib, khi_ib
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

c     initialize volume to zero
      volume = 0.0d0

      if (control_vol_sgn .gt. 0) then
      
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib
  
c             only include cell in max norm calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then

                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .gt. epsilon) then
                  volume = volume + dV
                elseif (phi_cur .gt. -epsilon) then
                  H = 0.5*(1 + phi_cur_over_epsilon 
     &                       + one_over_pi*sin(pi*phi_cur_over_epsilon))
                  volume = volume + H*dV
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

                phi_cur = phi(i,j,k)
                phi_cur_over_epsilon = phi_cur/epsilon
  
                if (phi_cur .gt. epsilon) then
                  volume = volume + dV
                elseif (phi_cur .gt. -epsilon) then
                  H = 0.5*(1 + phi_cur_over_epsilon 
     &                       + one_over_pi*sin(pi*phi_cur_over_epsilon))
                  volume = volume + H*dV
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
c
c  lsm3dSurfaceAreaZeroLevelSetControlVolume() computes the surface 
c  area of the surface defined by the zero level set within the 
c  computational domain.  The computational domain contains only those
c  cells that are included by the control volume data.
c
c  Arguments:
c    area (out):            area of the surface defined by the zero level 
c                           set
c    phi (in):              level set function
c    phi_* (in):            components of grad(phi)
c    control_vol (in):      control volume data (used to exclude cells
c                           from the integral calculation)
c    control_vol_sgn (in):  1 (-1) if positive (negative) control volume
c                           points should be used
c    dx, dy, dz (in):       grid spacing
c    epsilon (in):          width of numerical smoothing to use for 
c                           Heaviside function
c    *_gb (in):             index range for ghostbox
c    *_ib (in):             index range for interior box
c
c***********************************************************************
      subroutine lsm3dSurfaceAreaZeroLevelSetControlVolume(
     &  area,
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

      real area

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
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

c     initialize area to zero
      area = 0.0d0

      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in max norm calculation if it has a 
c             positive control volume
              if (control_vol(i,j,k) .gt. 0.d0) then

                phi_cur = phi(i,j,k)
    
                if (abs(phi_cur) .lt. epsilon) then
                  delta = 0.5d0*one_over_epsilon
     &                         *( 1+cos(pi*phi_cur*one_over_epsilon) ) 

                  norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

                  area = area + delta*norm_grad_phi*dV
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

                phi_cur = phi(i,j,k)
    
                if (abs(phi_cur) .lt. epsilon) then
                  delta = 0.5d0*one_over_epsilon
     &                         *( 1+cos(pi*phi_cur*one_over_epsilon) ) 

                  norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

                  area = area + delta*norm_grad_phi*dV
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
      subroutine lsm3dSurfaceAreaZeroLevelSetDeltaControlVolume(
     &  area,
     &  delta_phi,
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
     &  dx, dy, dz)
c***********************************************************************
c { begin subroutine
      implicit none

      real area

c     _gb refers to ghostbox 
c     _ib refers to box to include in integral calculation
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
      real delta_phi(ilo_phi_gb:ihi_phi_gb,
     &               jlo_phi_gb:jhi_phi_gb,
     &               klo_phi_gb:khi_phi_gb)
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
      integer i,j,k
      real norm_grad_phi
      real dV
      

c     compute dV = dx * dy * dz
      dV = dx * dy * dz

c     initialize area to zero
      area = 0.0d0

      if (control_vol_sgn .gt. 0) then
c       loop over included cells {
        do k=klo_ib,khi_ib
          do j=jlo_ib,jhi_ib
            do i=ilo_ib,ihi_ib

c             only include cell in max norm calculation if it has a 
c             positive control volume
              if ((control_vol(i,j,k) .gt. 0.d0) .and.
     &            (delta_phi(i,j,k) .gt. 0.d0)) then
 
                  norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

                  area = area + delta_phi(i,j,k)*norm_grad_phi*dV
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
              if ((control_vol(i,j,k) .lt. 0.d0) .and.
     &	            (delta_phi(i,j,k) .gt. 0.d0)) then

                  norm_grad_phi = sqrt(
     &                phi_x(i,j,k)*phi_x(i,j,k)
     &              + phi_y(i,j,k)*phi_y(i,j,k)
     &              + phi_z(i,j,k)*phi_z(i,j,k) )

                  area = area + delta_phi(i,j,k)*norm_grad_phi*dV

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
