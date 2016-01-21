c***********************************************************************
c
c  File:        lsm_curvature2d.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 routines for computing 2D curvature
c
c***********************************************************************


c***********************************************************************
      subroutine lsm2dComputeMeanCurvatureOrder2(
     &  kappa,
     &  ilo_kappa_gb, ihi_kappa_gb,
     &  jlo_kappa_gb, jhi_kappa_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  phi_x, phi_y, grad_phi_mag,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  dx,dy)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_kappa_gb, ihi_kappa_gb
      integer jlo_kappa_gb, jhi_kappa_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      real kappa(ilo_kappa_gb:ihi_kappa_gb,
     &                       jlo_kappa_gb:jhi_kappa_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                       jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                       jlo_grad_phi_gb:jhi_grad_phi_gb)
      real grad_phi_mag(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                              jlo_grad_phi_gb:jhi_grad_phi_gb)
      real dx, dy
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb    
         
c     local variables      
      integer i,j,l
      real phi_xx, phi_xy, phi_yy
      real dxsq_factor, dysq_factor, dxdy_factor
      real zero_tol, denominator
      parameter (zero_tol=1.d-11)
 
      dxsq_factor = 1.d0/dx/dx
      dysq_factor = 1.d0/dy/dy
      dxdy_factor = 0.25d0/dx/dy
      
c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          phi_xx = ( phi(i+1,j) - 2.0d0*phi(i,j) 
     &             + phi(i-1,j) )*dxsq_factor
          phi_yy = ( phi(i,j+1) - 2.0d0*phi(i,j)
     &             + phi(i,j-1) )*dysq_factor
          phi_xy = ( phi(i+1,j+1) - phi(i+1,j-1)
     &              -phi(i-1,j+1) + phi(i-1,j-1))*dxdy_factor

          denominator = grad_phi_mag(i,j)*grad_phi_mag(i,j)*
     &          grad_phi_mag(i,j)       

          if ( denominator .lt. zero_tol ) then
             kappa(i,j) = 0.d0
          else 
             kappa(i,j) =  phi_xx*phi_y(i,j)*phi_y(i,j)  
     &                  +  phi_yy*phi_x(i,j)*phi_x(i,j)  
     &                  -2*phi_xy*phi_x(i,j)*phi_y(i,j)
             kappa(i,j) = kappa(i,j)/denominator
          endif

        enddo
      enddo 
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************



c**********************************************************************
      subroutine lsm2dComputeMeanCurvature9Stencil(
     &  kappa,
     &  ilo_kappa_gb, ihi_kappa_gb,
     &  jlo_kappa_gb, jhi_kappa_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  dx,dy)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_kappa_gb, ihi_kappa_gb
      integer jlo_kappa_gb, jhi_kappa_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real kappa(ilo_kappa_gb:ihi_kappa_gb,
     &                         jlo_kappa_gb:jhi_kappa_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &                         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      
c     local variables      
      integer i,j
      real phi_x_plus, phi_x_minus
      real norm_x_plus, norm_x_minus
      real phi_y_plus, phi_y_minus
      real norm_y_plus, norm_y_minus
      real zero_tol, tmp, c
      parameter (zero_tol=1.d-11, c=1.d0/16.d0)

          
c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

        phi_x_plus = (phi(i+1,j) - phi(i,j))/dx
        tmp=(phi(i,j+1) - phi(i,j-1) +  phi(i+1,j+1) - phi(i+1,j-1))/dy
        norm_x_plus = sqrt( phi_x_plus * phi_x_plus + c*tmp*tmp)

        phi_x_minus = (phi(i,j) - phi(i-1,j))/dx
        tmp=(phi(i,j+1) - phi(i,j-1) +  phi(i-1,j+1) - phi(i-1,j-1))/dy
        norm_x_minus = sqrt( phi_x_minus * phi_x_minus + c*tmp*tmp)

        phi_y_plus = (phi(i,j+1) - phi(i,j))/dy
        tmp=(phi(i+1,j) - phi(i-1,j) + phi(i+1,j+1) - phi(i-1,j+1))/dx
        norm_y_plus = sqrt( phi_y_plus*phi_y_plus + c*tmp*tmp)

        phi_y_minus = (phi(i,j) - phi(i,j-1))/dy
        tmp=(phi(i+1,j) - phi(i-1,j) + phi(i+1,j-1) - phi(i-1,j-1))/dx
        norm_y_minus = sqrt( phi_y_minus*phi_y_minus + c*tmp*tmp)

        if((norm_x_plus .lt. zero_tol) .or.
     &     (norm_x_minus .lt. zero_tol) .or.
     &     (norm_y_plus .lt. zero_tol) .or. 
     &     (norm_y_minus .lt. zero_tol) ) then
          kappa(i,j) = 0.d0
        else 
          kappa(i,j)=(  phi_x_plus/norm_x_plus 
     &                - phi_x_minus/norm_x_minus )/dx
     &              +(  phi_y_plus/norm_y_plus
     &                - phi_y_minus/norm_y_minus )/dy
        endif

        enddo  
      enddo 
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dComputeMeanCurvatureSgnDist(
     &  kappa,
     &  ilo_kappa_gb, ihi_kappa_gb,
     &  jlo_kappa_gb, jhi_kappa_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,     
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb,
     &  dx,dy)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_kappa_gb, ihi_kappa_gb
      integer jlo_kappa_gb, jhi_kappa_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb      
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real kappa(ilo_kappa_gb:ihi_kappa_gb,
     &                       jlo_kappa_gb:jhi_kappa_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      
c     local variables      
      integer i,j,l
      real laplacian
      real zero_tol, denominator, one, zero      
      parameter (zero_tol=1.d-11, one=1.d0, zero = 0.d0)
 
      real inv_dx_sq
      real inv_dy_sq

c     compute denominator values
      inv_dx_sq = 1.0d0/dx/dx
      inv_dy_sq = 1.0d0/dy/dy
      
c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          laplacian = 
     &        inv_dx_sq *
     &        ( phi(i+1,j) - 2.0d0*phi(i,j) + phi(i-1,j) ) 
     &      + inv_dy_sq *
     &        ( phi(i,j+1) - 2.0d0*phi(i,j) + phi(i,j-1) )
      
          if ( abs(laplacian) .lt. zero_tol) then
            laplacian = zero
          endif 

          denominator = one - phi(i,j) * laplacian       

          if ( abs(denominator) .lt. zero_tol ) then
             denominator  = 1.d0
          endif

          kappa(i,j)= laplacian/denominator

        enddo  
      enddo 
c     } end loop over grid


      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dComputeMeanCurvatureOrder4(
     &  kappa,
     &  ilo_kappa_gb, ihi_kappa_gb,
     &  jlo_kappa_gb, jhi_kappa_gb,     
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  phi_x, phi_y, grad_phi_mag,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  ilo_fb, ihi_fb, 
     &  jlo_fb, jhi_fb,     
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

      integer ilo_kappa_gb, ihi_kappa_gb
      integer jlo_kappa_gb, jhi_kappa_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb      
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      real kappa(ilo_kappa_gb:ihi_kappa_gb,
     &                         jlo_kappa_gb:jhi_kappa_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb)
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                       jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                       jlo_grad_phi_gb:jhi_grad_phi_gb)
      real grad_phi_mag(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                              jlo_grad_phi_gb:jhi_grad_phi_gb)
      real dx, dy
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      
c     local variables      
      integer i,j
      real phi_xx, phi_xy, phi_yy
      real dxsq_factor, dy_factor, dysq_factor
      real zero_tol, tmp
      parameter (zero_tol=1.d-11)
      real sixteen, thirty, eight
      parameter (sixteen = 16.0d0, thirty=30.0d0)
      parameter (eight = 8.0d0)
 
c     compute denominator values
      dxsq_factor = 0.0833333333333333333333d0/dx/dx
      dysq_factor = 0.0833333333333333333333d0/dy/dy
      dy_factor   = 0.0833333333333333333333d0/dy
      
c       { begin loop over grid
        do j=jlo_fb,jhi_fb
          do i=ilo_fb,ihi_fb  
       
            phi_xx = ( -phi(i+2,j) + sixteen*phi(i+1,j)
     &                 -thirty*phi(i,j)
     &                 -phi(i-2,j) + sixteen*phi(i-1,j) )
     &               * dxsq_factor
            phi_yy = ( -phi(i,j+2) + sixteen*phi(i,j+1) 
     &                 -thirty*phi(i,j) 
     &                 -phi(i,j-2) + sixteen*phi(i,j-1) )
     &               * dysq_factor
            phi_xy = ( -phi_x(i,j+2) + eight*phi_x(i,j+1) 
     &                 +phi_x(i,j-2) - eight*phi_x(i,j-1) )
     &               * dy_factor

            tmp = grad_phi_mag(i,j)*grad_phi_mag(i,j)*
     &            grad_phi_mag(i,j)

            if ( tmp .lt. zero_tol ) then
              kappa(i,j) = 0.d0
            else
              kappa(i,j) =  phi_xx*phi_y(i,j)*phi_y(i,j) 
     &                   +  phi_yy*phi_x(i,j)*phi_x(i,j)
     &                   -2*phi_xy*phi_x(i,j)*phi_y(i,j)
              kappa(i,j) = kappa(i,j)/tmp
            endif

       enddo
      enddo 
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
