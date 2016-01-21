c***********************************************************************
c
c  File:        lsm_spatial_derivatives2d.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 routines for computing 2D ENO/WENO spatial derivatives
c
c***********************************************************************

c***********************************************************************
c The algorithms and notation in these subroutines closely follows
c the discussion in Osher & Fedkiw (2003).
c***********************************************************************

c***********************************************************************
c
c  lsm2dComputeDn() computes the n-th undivided differences in the
c  specified direction given the (n-1)-th undivided differences.  The 
c  undivided differences in cells with insufficient data is set to a 
c  large number.
c
c  Arguments:
c    Dn (out):           n-th undivided differences 
c    Dn_minus_one (in):  (n-1)-th undivided differences 
c    n (in):             order of undivided differences to compute
c    *_gb (in):          index range for ghostbox
c    *_fb (in):          index range for fillbox
c
c  NOTES:
c   - The index ranges for all ghostboxes and the fillbox should 
c     correspond to the index range for cell-centered data.
c   - The undivided differences for odd n are face-centered (i.e.
c     indices are of the form (i+1/2)).  In this situation, the array
c     index corresponding to the (i+1/2)-th undivided difference is
c     i (i.e. the index shifted down to the nearest integer index). 
c   - When n is odd, Dn is computed on the faces of the grid cells
c     specified by the fillbox indices.  The index range for the 
c     undivided differences to be computed is ilo_fb to (ihi_fb+1); 
c     that is, the number of undivided difference computed is equal
c     to the number of faces associated with the fillbox grid cells
c     (ihi_fb - ilo_fb + 2).
c   - The ghostbox for Dn_minus_one MUST be at least one ghostcell width
c     larger than the fillbox.
c
c***********************************************************************
      subroutine lsm2dComputeDn(
     &  Dn,
     &  ilo_Dn_gb, ihi_Dn_gb, 
     &  jlo_Dn_gb, jhi_Dn_gb,
     &  Dn_minus_one,
     &  ilo_Dn_minus_one_gb, ihi_Dn_minus_one_gb, 
     &  jlo_Dn_minus_one_gb, jhi_Dn_minus_one_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  n,
     &  dir)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_Dn_gb, ihi_Dn_gb, jlo_Dn_gb, jhi_Dn_gb
      integer ilo_Dn_minus_one_gb, ihi_Dn_minus_one_gb
      integer jlo_Dn_minus_one_gb, jhi_Dn_minus_one_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real Dn(ilo_Dn_gb:ihi_Dn_gb,jlo_Dn_gb:jhi_Dn_gb)
      real Dn_minus_one(ilo_Dn_minus_one_gb:ihi_Dn_minus_one_gb,
     &                  jlo_Dn_minus_one_gb:jhi_Dn_minus_one_gb)
      integer n
      integer dir 
      integer i,j
      integer offset(1:2)
      integer fillbox_shift(1:2)
      real sign_multiplier
      real big
      parameter (big=1.d10)

c     calculate offsets, fillbox shifts, and sign_multiplier used 
c     when computing undivided differences.
c     NOTE:  even and odd undivided differences are taken in
c            opposite order because of the discrepancy between
c            face- and cell-centered data.  the sign discrepancy 
c            is taken into account by sign_multiplier
      do i=1,2
        offset(i) = 0
        fillbox_shift(i) = 0
      enddo
      if (mod(n,2).eq.1) then
        offset(dir) = 1
        sign_multiplier = 1.0
        fillbox_shift(dir) = 1
      else
        offset(dir) = -1
        sign_multiplier = -1.0
        fillbox_shift(dir) = 0
      endif

c     loop over cells with sufficient data {
      do j=jlo_fb,jhi_fb+fillbox_shift(2)
        do i=ilo_fb,ihi_fb+fillbox_shift(1)

          Dn(i,j) = sign_multiplier
     &            * ( Dn_minus_one(i,j)
     &              - Dn_minus_one(i-offset(1),j-offset(2)))

        enddo
      enddo
c     } end loop over grid 

c     set undivided differences for cells with insufficient data to big {
      do j=jlo_Dn_gb,jhi_Dn_gb
        do i=ilo_Dn_gb,ilo_fb-1
          Dn(i,j) = big
        enddo
      enddo

      do j=jlo_Dn_gb,jhi_Dn_gb
        do i=ihi_fb+fillbox_shift(1)+1,ihi_Dn_gb
          Dn(i,j) = big
        enddo
      enddo

      do j=jlo_Dn_gb,jlo_fb-1
        do i=ilo_Dn_gb,ihi_Dn_gb
          Dn(i,j) = big
        enddo
      enddo

      do j=jhi_fb+fillbox_shift(2)+1,jhi_Dn_gb
        do i=ilo_Dn_gb,ihi_Dn_gb
          Dn(i,j) = big
        enddo
      enddo

c     } end setting big value for cells near boundary of ghostcell box

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dHJENO1() computes the forward (plus) and backward (minus)
c  first-order Hamilton-Jacobi ENO approximations to the gradient of 
c  phi.
c
c  Arguments:
c    phi_*_plus (out):   components of grad(phi) in plus direction
c    phi_*_minus (out):  components of grad(phi) in minus direction
c    phi (in):           phi
c    D1 (in):            scratch space for holding undivided first-differences
c    dx, dy (in):        grid spacing
c    *_gb (in):          index range for ghostbox
c    *_fb (in):          index range for fillbox
c
c  NOTES:
c   - it is assumed that BOTH the plus AND minus derivatives have
c     the same fillbox
c
c***********************************************************************
      subroutine lsm2dHJENO1(
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb, 
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_plus_gb refers to ghostbox for grad_phi plus data
c     _grad_phi_minus_gb refers to ghostbox for grad_phi minus data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb)
      real dx, dy
      real inv_dx, inv_dy
      integer i,j
      integer order
      parameter (order=1)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)


c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------

c     compute first undivided differences in x-direction
      call lsm2dComputeDn(D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb,
     &                    jlo_fb, jhi_fb, 
     &                    order, x_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          phi_x_plus(i,j) = D1(i+1,j)*inv_dx
          phi_x_minus(i,j) = D1(i,j)*inv_dx
   
        enddo
      enddo
c     } end loop over grid 

c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------

c     compute first undivided differences in y-direction
      call lsm2dComputeDn(D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb, 
     &                    jlo_fb, jhi_fb, 
     &                    order, y_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          phi_y_plus(i,j) = D1(i,j+1)*inv_dy
          phi_y_minus(i,j) = D1(i,j)*inv_dy
   
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dHJENO2() computes the forward (plus) and backward (minus)
c  second-order Hamilton-Jacobi ENO approximations to the gradient of 
c  phi.
c
c  Arguments:
c    phi_*_plus (out):   components of grad(phi) in plus direction
c    phi_*_minus (out):  components of grad(phi) in minus direction
c    phi (in):           phi
c    D1 (in):            scratch space for holding undivided first-differences
c    D2 (in):            scratch space for holding undivided second-differences
c    dx, dy (in):        grid spacing
c    *_gb (in):          index range for ghostbox
c    *_fb (in):          index range for fillbox
c
c  NOTES:
c   - it is assumed that BOTH the plus AND minus derivatives have
c     the same fillbox
c
c***********************************************************************
      subroutine lsm2dHJENO2(
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb, 
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb,
     &  D2,
     &  ilo_D2_gb, ihi_D2_gb, jlo_D2_gb, jhi_D2_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_plus_gb refers to ghostbox for grad_phi plus data
c     _grad_phi_minus_gb refers to ghostbox for grad_phi minus data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
      integer ilo_D2_gb, ihi_D2_gb, jlo_D2_gb, jhi_D2_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb)
      real D2(ilo_D2_gb:ihi_D2_gb,
     &        jlo_D2_gb:jhi_D2_gb)
      real dx, dy
      real inv_dx, inv_dy
      integer i, j
      real half
      parameter (half=0.5d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)


c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------

c     compute first undivided differences in x-direction
      call lsm2dComputeDn(D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb-1, ihi_fb+1, 
     &                    jlo_fb, jhi_fb,
     &                    order_1, x_dir)

c     compute second undivided differences in x-direction
      call lsm2dComputeDn(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    ilo_fb-1, ihi_fb+1, 
     &                    jlo_fb, jhi_fb, 
     &                    order_2, x_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         phi_x_plus
          if (abs(D2(i,j)).lt.abs(D2(i+1,j))) then
            phi_x_plus(i,j) = (D1(i+1,j) - half*D2(i,j))*inv_dx
          else
            phi_x_plus(i,j) = (D1(i+1,j) - half*D2(i+1,j))*inv_dx
          endif

c         phi_x_minus 
          if (abs(D2(i-1,j)).lt.abs(D2(i,j))) then
            phi_x_minus(i,j) = (D1(i,j) + half*D2(i-1,j))*inv_dx
          else
            phi_x_minus(i,j) = (D1(i,j) + half*D2(i,j))*inv_dx
          endif

        enddo
      enddo
c     } end loop over grid 

c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------

c     compute first undivided differences in y-direction
      call lsm2dComputeDn(D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb, 
     &                    jlo_fb-1, jhi_fb+1, 
     &                    order_1, y_dir)

c     compute second undivided differences in y-direction
      call lsm2dComputeDn(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    ilo_fb, ihi_fb, 
     &                    jlo_fb-1, jhi_fb+1, 
     &                    order_2, y_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         phi_y_plus
          if (abs(D2(i,j)).lt.abs(D2(i,j+1))) then
            phi_y_plus(i,j) = (D1(i,j+1) - half*D2(i,j))*inv_dy
          else
            phi_y_plus(i,j) = (D1(i,j+1) - half*D2(i,j+1))*inv_dy
          endif

c         phi_y_minus
          if (abs(D2(i,j-1)).lt.abs(D2(i,j))) then
            phi_y_minus(i,j) = (D1(i,j) + half*D2(i,j-1))*inv_dy
          else
            phi_y_minus(i,j) = (D1(i,j) + half*D2(i,j))*inv_dy
          endif

        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dHJENO3() computes the forward (plus) and backward (minus)
c  third-order Hamilton-Jacobi ENO approximations to the gradient of 
c  phi.
c
c  Arguments:
c    phi_*_plus (out):   components of grad(phi) in plus direction
c    phi_*_minus (out):  components of grad(phi) in minus direction
c    phi (in):           phi
c    D1 (in):            scratch space for holding undivided first-differences
c    D2 (in):            scratch space for holding undivided second-differences
c    D3 (in):            scratch space for holding undivided third-differences
c    dx, dy (in):        grid spacing
c    *_gb (in):          index range for ghostbox
c    *_fb (in):          index range for fillbox
c
c  NOTES:
c   - it is assumed that BOTH the plus AND minus derivatives have
c     the same fillbox
c
c***********************************************************************
      subroutine lsm2dHJENO3(
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb, 
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb,
     &  D2,
     &  ilo_D2_gb, ihi_D2_gb, jlo_D2_gb, jhi_D2_gb,
     &  D3,
     &  ilo_D3_gb, ihi_D3_gb, jlo_D3_gb, jhi_D3_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_plus_gb refers to ghostbox for grad_phi plus data
c     _grad_phi_minus_gb refers to ghostbox for grad_phi minus data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
      integer ilo_D2_gb, ihi_D2_gb, jlo_D2_gb, jhi_D2_gb
      integer ilo_D3_gb, ihi_D3_gb, jlo_D3_gb, jhi_D3_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb)
      real D2(ilo_D2_gb:ihi_D2_gb,
     &        jlo_D2_gb:jhi_D2_gb)
      real D3(ilo_D3_gb:ihi_D3_gb,
     &        jlo_D3_gb:jhi_D3_gb)
      real dx, dy
      real inv_dx, inv_dy
      integer i,j
      real zero, half, third, sixth
      parameter (zero=0.0d0, half=0.5d0, third=1.d0/3.d0)
      parameter (sixth=1.d0/6.d0)
      integer order_1, order_2, order_3
      parameter (order_1=1,order_2=2,order_3=3)
      integer x_dir, y_dir
      parameter (x_dir=1, y_dir=2)


c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------

c     compute first undivided differences in x-direction
      call lsm2dComputeDn(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb-2, ihi_fb+2,
     &                    jlo_fb, jhi_fb,
     &                    order_1, x_dir)

c     compute second undivided differences in x-direction
      call lsm2dComputeDn(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    ilo_fb-2, ihi_fb+2,
     &                    jlo_fb, jhi_fb,
     &                    order_2, x_dir)

c     compute third undivided differences in x-direction
      call lsm2dComputeDn(D3, 
     &                    ilo_D3_gb, ihi_D3_gb, 
     &                    jlo_D3_gb, jhi_D3_gb, 
     &                    D2,
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    ilo_fb-1, ihi_fb+1, 
     &                    jlo_fb, jhi_fb,
     &                    order_3, x_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         { begin calculation of phi_x_plus
          phi_x_plus(i,j) = D1(i+1,j)

          if (abs(D2(i,j)).lt.abs(D2(i+1,j))) then
            phi_x_plus(i,j) = phi_x_plus(i,j) - half*D2(i,j) 
            if (abs(D3(i,j)).lt.abs(D3(i+1,j))) then
              phi_x_plus(i,j) = phi_x_plus(i,j) - sixth*D3(i,j)
            else
              phi_x_plus(i,j) = phi_x_plus(i,j) - sixth*D3(i+1,j)
            endif
          else
            phi_x_plus(i,j) = phi_x_plus(i,j) - half*D2(i+1,j) 
            if (abs(D3(i+1,j)).lt.abs(D3(i+2,j))) then
              phi_x_plus(i,j) = phi_x_plus(i,j) + third*D3(i+1,j)
            else
              phi_x_plus(i,j) = phi_x_plus(i,j) + third*D3(i+2,j)
            endif
          endif
  
c         divide phi_x_plus by dx
          phi_x_plus(i,j) = phi_x_plus(i,j)*inv_dx

c         } end calculation of phi_x_plus

c         { begin calculation of phi_x_minus
          phi_x_minus(i,j) = D1(i,j)

          if (abs(D2(i-1,j)).lt.abs(D2(i,j))) then
            phi_x_minus(i,j) = phi_x_minus(i,j) + half*D2(i-1,j) 
            if (abs(D3(i-1,j)).lt.abs(D3(i,j))) then
              phi_x_minus(i,j) = phi_x_minus(i,j) + third*D3(i-1,j)
            else
              phi_x_minus(i,j) = phi_x_minus(i,j) + third*D3(i,j)
            endif
          else
            phi_x_minus(i,j) = phi_x_minus(i,j) + half*D2(i,j) 
            if (abs(D3(i,j)).lt.abs(D3(i+1,j))) then
              phi_x_minus(i,j) = phi_x_minus(i,j) - sixth*D3(i,j)
            else
              phi_x_minus(i,j) = phi_x_minus(i,j) - sixth*D3(i+1,j)
            endif
          endif

c         divide phi_x_minus by dx
          phi_x_minus(i,j) = phi_x_minus(i,j)*inv_dx

c         } end calculation of phi_x_minus

        enddo
      enddo
c     } end loop over grid 


c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDn(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb,
     &                    jlo_fb-2, jhi_fb+2,
     &                    order_1, y_dir)

c     compute second undivided differences in y-direction
      call lsm2dComputeDn(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    ilo_fb, ihi_fb,
     &                    jlo_fb-2, jhi_fb+2,
     &                    order_2, y_dir)

c     compute third undivided differences in y-direction
      call lsm2dComputeDn(D3, 
     &                    ilo_D3_gb, ihi_D3_gb, 
     &                    jlo_D3_gb, jhi_D3_gb, 
     &                    D2,
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    ilo_fb, ihi_fb,
     &                    jlo_fb-1, jhi_fb+1,
     &                    order_3, y_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         { begin calculation of phi_y_plus
          phi_y_plus(i,j) = D1(i,j+1)

          if (abs(D2(i,j)).lt.abs(D2(i,j+1))) then
            phi_y_plus(i,j) = phi_y_plus(i,j) - half*D2(i,j) 
            if (abs(D3(i,j)).lt.abs(D3(i,j+1))) then
              phi_y_plus(i,j) = phi_y_plus(i,j) - sixth*D3(i,j)
            else
              phi_y_plus(i,j) = phi_y_plus(i,j) - sixth*D3(i,j+1)
            endif
          else
            phi_y_plus(i,j) = phi_y_plus(i,j) - half*D2(i,j+1) 
            if (abs(D3(i,j+1)).lt.abs(D3(i,j+2))) then
              phi_y_plus(i,j) = phi_y_plus(i,j) + third*D3(i,j+1)
            else
              phi_y_plus(i,j) = phi_y_plus(i,j) + third*D3(i,j+2)
            endif
          endif

c         divide phi_y_plus by dy
          phi_y_plus(i,j) = phi_y_plus(i,j)*inv_dy

c         } end calculation of phi_y_plus

c         { begin calculation of phi_y_minus
          phi_y_minus(i,j) = D1(i,j)

          if (abs(D2(i,j-1)).lt.abs(D2(i,j))) then
            phi_y_minus(i,j) = phi_y_minus(i,j) + half*D2(i,j-1) 
            if (abs(D3(i,j-1)).lt.abs(D3(i,j))) then
              phi_y_minus(i,j) = phi_y_minus(i,j) + third*D3(i,j-1)
            else
              phi_y_minus(i,j) = phi_y_minus(i,j) + third*D3(i,j)
            endif
          else
            phi_y_minus(i,j) = phi_y_minus(i,j) + half*D2(i,j) 
            if (abs(D3(i,j)).lt.abs(D3(i,j+1))) then
              phi_y_minus(i,j) = phi_y_minus(i,j) - sixth*D3(i,j)
            else
              phi_y_minus(i,j) = phi_y_minus(i,j) - sixth*D3(i,j+1)
            endif
          endif

c         divide phi_y_minus by dy
          phi_y_minus(i,j) = phi_y_minus(i,j)*inv_dy

c         } end calculation of phi_y_minus

        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dHJWENO5() computes the forward (plus) and backward (minus)
c  fifth-order Hamilton-Jacobi ENO approximations to the gradient of 
c  phi.
c
c  Arguments:
c    phi_*_plus (out):   components of grad(phi) in plus direction
c    phi_*_minus (out):  components of grad(phi) in minus direction
c    phi (in):           phi
c    D1 (in):            scratch space for holding undivided first-differences
c    dx, dy (in):        grid spacing
c    *_gb (in):          index range for ghostbox
c    *_fb (in):          index range for fillbox
c
c  NOTES:
c   - it is assumed that BOTH the plus AND minus derivatives have
c     the same fillbox
c
c***********************************************************************
      subroutine lsm2dHJWENO5(
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb, 
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_plus_gb refers to ghostbox for grad_phi plus data
c     _grad_phi_minus_gb refers to ghostbox for grad_phi minus data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb)
      real dx, dy
      real inv_dx, inv_dy

c     variables for WENO calculation 
      real v1,v2,v3,v4,v5
      real S1,S2,S3
      real a1,a2,a3, inv_sum_a
      real phi_x_1,phi_x_2,phi_x_3
      real phi_y_1,phi_y_2,phi_y_3
      real tiny_nonzero_number
      parameter (tiny_nonzero_number=1.d-99)
      real eps
      real one_third, seven_sixths, eleven_sixths
      real one_sixth, five_sixths
      real thirteen_twelfths, one_fourth
      parameter (one_third=1.d0/3.d0)
      parameter (seven_sixths=7.d0/6.d0)
      parameter (eleven_sixths=11.d0/6.d0) 
      parameter (one_sixth=1.d0/6.d0)
      parameter (five_sixths=5.d0/6.d0)
      parameter (thirteen_twelfths=13.d0/12.d0)
      parameter (one_fourth=0.25d0)

      integer i,j
      real zero
      parameter (zero=0.0d0)
      integer order_1
      parameter (order_1=1)
      integer x_dir, y_dir
      parameter (x_dir=1, y_dir=2)


c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDn(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb-2, ihi_fb+2,
     &                    jlo_fb, jhi_fb,
     &                    order_1, x_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         { begin calculation of phi_x_plus
  
c         extract v1,v2,v3,v4,v5 from D1
          v1 = D1(i+3,j)*inv_dx
          v2 = D1(i+2,j)*inv_dx
          v3 = D1(i+1,j)*inv_dx
          v4 = D1(i,j)*inv_dx
          v5 = D1(i-1,j)*inv_dx
  
c         WENO5 algorithm for current grid point using appropriate
c         upwind values for v1,...,v5
 
c         compute eps for current grid point
          eps = 1e-6*max(v1*v1,v2*v2,v3*v3,v4*v4,v5*v5)
     &        + tiny_nonzero_number

c         compute the phi_x_1, phi_x_2, phi_x_3
          phi_x_1 = one_third*v1 - seven_sixths*v2 + eleven_sixths*v3
          phi_x_2 = -one_sixth*v2 + five_sixths*v3 + one_third*v4
          phi_x_3 = one_third*v3 + five_sixths*v4 - one_sixth*v5
 
c         compute the smoothness measures
          S1 = thirteen_twelfths*(v1-2.d0*v2+v3)**2
     &       + one_fourth*(v1-4.d0*v2+3.d0*v3)**2
          S2 = thirteen_twelfths*(v2-2.d0*v3+v4)**2
     &       + one_fourth*(v2-v4)**2
          S3 = thirteen_twelfths*(v3-2.d0*v4+v5)**2
     &       + one_fourth*(3.d0*v3-4.d0*v4+v5)**2

c         compute normalized weights
          a1 = 0.1d0/(S1+eps)**2
          a2 = 0.6d0/(S2+eps)**2
          a3 = 0.3d0/(S3+eps)**2
          inv_sum_a = 1.0d0 / (a1 + a2 + a3)
          a1 = a1*inv_sum_a
          a2 = a2*inv_sum_a
          a3 = a3*inv_sum_a
  
c         compute phi_x_plus 
          phi_x_plus(i,j) = a1*phi_x_1 + a2*phi_x_2 + a3*phi_x_3
  
c         } end calculation of phi_x_plus

c         { begin calculation of phi_x_minus
  
c         extract v1,v2,v3,v4,v5 from D1
          v1 = D1(i-2,j)*inv_dx
          v2 = D1(i-1,j)*inv_dx
          v3 = D1(i,j)*inv_dx
          v4 = D1(i+1,j)*inv_dx
          v5 = D1(i+2,j)*inv_dx
 
c         WENO5 algorithm for current grid point using appropriate
c         upwind values for v1,...,v5
  
c         compute eps for current grid point
          eps = 1e-6*max(v1*v1,v2*v2,v3*v3,v4*v4,v5*v5)
     &        + tiny_nonzero_number

c         compute the phi_x_1, phi_x_2, phi_x_3
          phi_x_1 = one_third*v1 - seven_sixths*v2 + eleven_sixths*v3
          phi_x_2 = -one_sixth*v2 + five_sixths*v3 + one_third*v4
          phi_x_3 = one_third*v3 + five_sixths*v4 - one_sixth*v5
 
c         compute the smoothness measures
          S1 = thirteen_twelfths*(v1-2.d0*v2+v3)**2
     &       + one_fourth*(v1-4.d0*v2+3.d0*v3)**2
          S2 = thirteen_twelfths*(v2-2.d0*v3+v4)**2
     &       + one_fourth*(v2-v4)**2
          S3 = thirteen_twelfths*(v3-2.d0*v4+v5)**2
     &       + one_fourth*(3.d0*v3-4.d0*v4+v5)**2

c         compute normalized weights
          a1 = 0.1d0/(S1+eps)**2
          a2 = 0.6d0/(S2+eps)**2
          a3 = 0.3d0/(S3+eps)**2
          inv_sum_a = 1.0d0 / (a1 + a2 + a3)
          a1 = a1*inv_sum_a
          a2 = a2*inv_sum_a
          a3 = a3*inv_sum_a
 
c         compute phi_x_minus 
          phi_x_minus(i,j) = a1*phi_x_1 + a2*phi_x_2 + a3*phi_x_3
 
c         } end calculation of phi_x_minus

        enddo
      enddo
c     } end loop over grid 


c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDn(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb,  
     &                    jlo_fb-2, jhi_fb+2,
     &                    order_1, y_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         { begin calculation of phi_y_plus
  
c         extract v1,v2,v3,v4,v5 from D1
          v1 = D1(i,j+3)*inv_dy
          v2 = D1(i,j+2)*inv_dy
          v3 = D1(i,j+1)*inv_dy
          v4 = D1(i,j)*inv_dy
          v5 = D1(i,j-1)*inv_dy
 
c         WENO5 algorithm for current grid point using appropriate
c         upwind values for v1,...,v5
  
c         compute eps for current grid point
          eps = 1e-6*max(v1*v1,v2*v2,v3*v3,v4*v4,v5*v5)
     &        + tiny_nonzero_number

c         compute the phi_y_1, phi_y_2, phi_y_3
          phi_y_1 = one_third*v1 - seven_sixths*v2 + eleven_sixths*v3
          phi_y_2 = -one_sixth*v2 + five_sixths*v3 + one_third*v4
          phi_y_3 = one_third*v3 + five_sixths*v4 - one_sixth*v5
  
c         compute the smoothness measures
          S1 = thirteen_twelfths*(v1-2.d0*v2+v3)**2
     &       + one_fourth*(v1-4.d0*v2+3.d0*v3)**2
          S2 = thirteen_twelfths*(v2-2.d0*v3+v4)**2
     &       + one_fourth*(v2-v4)**2
          S3 = thirteen_twelfths*(v3-2.d0*v4+v5)**2
     &     + one_fourth*(3.d0*v3-4.d0*v4+v5)**2

c         compute normalized weights
          a1 = 0.1d0/(S1+eps)**2
          a2 = 0.6d0/(S2+eps)**2
          a3 = 0.3d0/(S3+eps)**2
          inv_sum_a = 1.0d0 / (a1 + a2 + a3)
          a1 = a1*inv_sum_a
          a2 = a2*inv_sum_a
          a3 = a3*inv_sum_a
 
c         compute phi_y_plus 
          phi_y_plus(i,j) = a1*phi_y_1 + a2*phi_y_2 + a3*phi_y_3
  
c         } end calculation of phi_y_plus

c         { begin calculation of phi_y_minus
  
c         extract v1,v2,v3,v4,v5 from D1
          v1 = D1(i,j-2)*inv_dy
          v2 = D1(i,j-1)*inv_dy
          v3 = D1(i,j)*inv_dy
          v4 = D1(i,j+1)*inv_dy
          v5 = D1(i,j+2)*inv_dy
 
c         WENO5 algorithm for current grid point using appropriate
c         upwind values for v1,...,v5
  
c         compute eps for current grid point
          eps = 1e-6*max(v1*v1,v2*v2,v3*v3,v4*v4,v5*v5)
     &        + tiny_nonzero_number

c         compute the phi_y_1, phi_y_2, phi_y_3
          phi_y_1 = one_third*v1 - seven_sixths*v2 + eleven_sixths*v3
          phi_y_2 = -one_sixth*v2 + five_sixths*v3 + one_third*v4
          phi_y_3 = one_third*v3 + five_sixths*v4 - one_sixth*v5
 
c         compute the smoothness measures
          S1 = thirteen_twelfths*(v1-2.d0*v2+v3)**2
     &       + one_fourth*(v1-4.d0*v2+3.d0*v3)**2
          S2 = thirteen_twelfths*(v2-2.d0*v3+v4)**2
     &       + one_fourth*(v2-v4)**2
          S3 = thirteen_twelfths*(v3-2.d0*v4+v5)**2
     &       + one_fourth*(3.d0*v3-4.d0*v4+v5)**2

c         compute normalized weights
          a1 = 0.1d0/(S1+eps)**2
          a2 = 0.6d0/(S2+eps)**2
          a3 = 0.3d0/(S3+eps)**2
          inv_sum_a = 1.0d0 / (a1 + a2 + a3)
          a1 = a1*inv_sum_a
          a2 = a2*inv_sum_a
          a3 = a3*inv_sum_a
  
c         compute phi_y_minus 
          phi_y_minus(i,j) = a1*phi_y_1 + a2*phi_y_2 + a3*phi_y_3
  
c         } end calculation of phi_y_minus

        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dUpwindHJENO1() computes the first-order Hamilton-Jacobi ENO 
c  upwind approximation to the gradient of phi.
c
c  Arguments:
c    phi_* (out):  components of grad(phi)
c    phi (in):     phi
c    vel_* (in):   components of the velocity
c    D1 (in):      scratch space for holding undivided first-differences
c    dx, dy (in):  grid spacing
c    *_gb (in):    index range for ghostbox
c    *_fb (in):    index range for fillbox
c
c***********************************************************************
      subroutine lsm2dUpwindHJENO1(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  vel_x, vel_y,
     &  ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _vel_gb refers to ghostbox for velocity data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb
      integer ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb)
      real dx, dy
      real inv_dx, inv_dy
      integer i,j
      real zero
      parameter (zero=0.0d0)
      real zero_tol
      parameter (zero_tol=1.d-11)
      integer order
      parameter (order=1)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)


c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute upwind phi_x
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDn(D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb,
     &                    jlo_fb, jhi_fb, 
     &                    order, x_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         phi_x
          if (abs(vel_x(i,j)) .lt. zero_tol) then
c           vel_x == 0
            phi_x(i,j) = zero
          elseif (vel_x(i,j) .gt. 0) then
c           vel_x > 0
            phi_x(i,j) = D1(i,j)*inv_dx
          else
c           vel_x < 0
            phi_x(i,j) = D1(i+1,j)*inv_dx
          endif

        enddo
      enddo
c     } end loop over grid 


c----------------------------------------------------
c    compute upwind phi_y
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDn(D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb, 
     &                    jlo_fb, jhi_fb, 
     &                    order, y_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         phi_y
          if (abs(vel_y(i,j)) .lt. zero_tol) then
c           vel_y == 0
            phi_y(i,j) = zero
          elseif (vel_y(i,j) .gt. 0) then
c           vel_y > 0
            phi_y(i,j) = D1(i,j)*inv_dy
          else
c           vel_y < 0
            phi_y(i,j) = D1(i,j+1)*inv_dy
          endif
   
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dUpwindHJENO2() computes the second-order Hamilton-Jacobi ENO 
c  upwind approximation to the gradient of phi.
c
c  Arguments:
c    phi_* (out):  components of grad(phi)
c    phi (in):     phi
c    vel_* (in):   components of the velocity
c    D1 (in):      scratch space for holding undivided first-differences
c    D2 (in):      scratch space for holding undivided second-differences
c    dx, dy (in):  grid spacing
c    *_gb (in):    index range for ghostbox
c    *_fb (in):    index range for fillbox
c
c***********************************************************************
      subroutine lsm2dUpwindHJENO2(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  vel_x, vel_y,
     &  ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb,
     &  D2,
     &  ilo_D2_gb, ihi_D2_gb, jlo_D2_gb, jhi_D2_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _vel_gb refers to ghostbox for velocity data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb
      integer ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
      integer ilo_D2_gb, ihi_D2_gb, jlo_D2_gb, jhi_D2_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb)
      real D2(ilo_D2_gb:ihi_D2_gb,
     &        jlo_D2_gb:jhi_D2_gb)
      real dx, dy
      real inv_dx, inv_dy
      integer i, j
      real zero, half
      parameter (zero=0.0d0, half=0.5d0)
      real zero_tol
      parameter (zero_tol=1.d-11)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)


c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute upwind phi_x
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDn(D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb-1, ihi_fb+1,
     &                    jlo_fb, jhi_fb, 
     &                    order_1, x_dir)

c     compute second undivided differences in x-direction
      call lsm2dComputeDn(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    ilo_fb-1, ihi_fb+1, 
     &                    jlo_fb, jhi_fb, 
     &                    order_2, x_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         phi_x
          if (abs(vel_x(i,j)) .lt. zero_tol) then

c           vel_x == 0
            phi_x(i,j) = zero

          elseif (vel_x(i,j) .gt. 0) then

c           vel_x > 0
            if (abs(D2(i-1,j)).lt.abs(D2(i,j))) then
              phi_x(i,j) = (D1(i,j) + half*D2(i-1,j))*inv_dx
            else
              phi_x(i,j) = (D1(i,j) + half*D2(i,j))*inv_dx
            endif

          else

c           vel_x < 0
            if (abs(D2(i,j)).lt.abs(D2(i+1,j))) then
              phi_x(i,j) = (D1(i+1,j) - half*D2(i,j))*inv_dx
            else
              phi_x(i,j) = (D1(i+1,j) - half*D2(i+1,j))*inv_dx
            endif

          endif

        enddo
      enddo
c     } end loop over grid 


c----------------------------------------------------
c    compute upwind phi_y
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDn(D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb, 
     &                    jlo_fb-1, jhi_fb+1,
     &                    order_1, y_dir)

c     compute second undivided differences in y-direction
      call lsm2dComputeDn(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    ilo_fb, ihi_fb, 
     &                    jlo_fb-1, jhi_fb+1, 
     &                    order_2, y_dir)


c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         phi_y
          if (abs(vel_y(i,j)) .lt. zero_tol) then

c           vel_y == 0
            phi_y(i,j) = zero

          elseif (vel_y(i,j) .gt. 0) then

c           vel_y > 0
            if (abs(D2(i,j-1)).lt.abs(D2(i,j))) then
              phi_y(i,j) = (D1(i,j) + half*D2(i,j-1))*inv_dy
            else
              phi_y(i,j) = (D1(i,j) + half*D2(i,j))*inv_dy
            endif

          else

c           vel_y < 0
            if (abs(D2(i,j)).lt.abs(D2(i,j+1))) then
              phi_y(i,j) = (D1(i,j+1) - half*D2(i,j))*inv_dy
            else
              phi_y(i,j) = (D1(i,j+1) - half*D2(i,j+1))*inv_dy
            endif

          endif

        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dUpwindHJENO3() computes the third-order Hamilton-Jacobi ENO
c  upwind approximation to the gradient of phi.
c
c  Arguments:
c    phi_* (out):  components of grad(phi)
c    phi (in):     phi
c    vel_* (in):   components of the velocity
c    D1 (in):      scratch space for holding undivided first-differences
c    D2 (in):      scratch space for holding undivided second-differences
c    D3 (in):      scratch space for holding undivided third-differences
c    dx, dy (in):  grid spacing
c    *_gb (in):    index range for ghostbox
c    *_fb (in):    index range for fillbox
c
c***********************************************************************
      subroutine lsm2dUpwindHJENO3(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  vel_x, vel_y,
     &  ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb,
     &  D2,
     &  ilo_D2_gb, ihi_D2_gb, jlo_D2_gb, jhi_D2_gb,
     &  D3,
     &  ilo_D3_gb, ihi_D3_gb, jlo_D3_gb, jhi_D3_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _vel_gb refers to ghostbox for velocity data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb
      integer ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
      integer ilo_D2_gb, ihi_D2_gb, jlo_D2_gb, jhi_D2_gb
      integer ilo_D3_gb, ihi_D3_gb, jlo_D3_gb, jhi_D3_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb)
      real D2(ilo_D2_gb:ihi_D2_gb,
     &        jlo_D2_gb:jhi_D2_gb)
      real D3(ilo_D3_gb:ihi_D3_gb,
     &        jlo_D3_gb:jhi_D3_gb)
      real dx, dy
      real inv_dx, inv_dy
      integer i,j
      real zero, half, third, sixth
      parameter (zero=0.0d0, half=0.5d0, third=1.d0/3.d0)
      parameter (sixth=1.d0/6.d0)
      real zero_tol
      parameter (zero_tol=1.d-11)
      integer order_1, order_2, order_3
      parameter (order_1=1,order_2=2,order_3=3)
      integer x_dir, y_dir
      parameter (x_dir=1, y_dir=2)


c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute upwind phi_x 
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDn(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb-2, ihi_fb+2,
     &                    jlo_fb, jhi_fb,
     &                    order_1, x_dir)

c     compute second undivided differences in x-direction
      call lsm2dComputeDn(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    ilo_fb-2, ihi_fb+2,
     &                    jlo_fb, jhi_fb,
     &                    order_2, x_dir)

c     compute third undivided differences in x-direction
      call lsm2dComputeDn(D3, 
     &                    ilo_D3_gb, ihi_D3_gb, 
     &                    jlo_D3_gb, jhi_D3_gb, 
     &                    D2,
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    ilo_fb-1, ihi_fb+1, 
     &                    jlo_fb, jhi_fb, 
     &                    order_3, x_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         { begin calculation of phi_x
          if (abs(vel_x(i,j)) .lt. zero_tol) then

c           vel_x == 0
            phi_x(i,j) = zero

          elseif (vel_x(i,j) .gt. 0) then

c           vel_x > 0
            phi_x(i,j) = D1(i,j)
            if (abs(D2(i-1,j)).lt.abs(D2(i,j))) then
              phi_x(i,j) = phi_x(i,j) + half*D2(i-1,j) 
              if (abs(D3(i-1,j)).lt.abs(D3(i,j))) then
                phi_x(i,j) = phi_x(i,j) + third*D3(i-1,j)
              else
                phi_x(i,j) = phi_x(i,j) + third*D3(i,j)
              endif
            else
              phi_x(i,j) = phi_x(i,j) + half*D2(i,j) 
              if (abs(D3(i,j)).lt.abs(D3(i+1,j))) then
                phi_x(i,j) = phi_x(i,j) - sixth*D3(i,j)
              else
                phi_x(i,j) = phi_x(i,j) - sixth*D3(i+1,j)
              endif
            endif

          else

c           vel_x < 0
            phi_x(i,j) = D1(i+1,j)
  
            if (abs(D2(i,j)).lt.abs(D2(i+1,j))) then
              phi_x(i,j) = phi_x(i,j) - half*D2(i,j) 
              if (abs(D3(i,j)).lt.abs(D3(i+1,j))) then
                phi_x(i,j) = phi_x(i,j) - sixth*D3(i,j)
              else
                phi_x(i,j) = phi_x(i,j) - sixth*D3(i+1,j)
              endif
            else
              phi_x(i,j) = phi_x(i,j) - half*D2(i+1,j) 
              if (abs(D3(i+1,j)).lt.abs(D3(i+2,j))) then
                phi_x(i,j) = phi_x(i,j) + third*D3(i+1,j)
              else
                phi_x(i,j) = phi_x(i,j) + third*D3(i+2,j)
              endif
            endif
  
          endif

c         divide phi_x by dx
          phi_x(i,j) = phi_x(i,j)*inv_dx

c         } end calculation of phi_x

        enddo
      enddo
c     } end loop over grid 


c----------------------------------------------------
c    compute upwind phi_y
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDn(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb,
     &                    jlo_fb-2, jhi_fb+2,
     &                    order_1, y_dir)

c     compute second undivided differences in y-direction
      call lsm2dComputeDn(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    ilo_fb, ihi_fb,
     &                    jlo_fb-2, jhi_fb+2,
     &                    order_2, y_dir)

c     compute third undivided differences in y-direction
      call lsm2dComputeDn(D3, 
     &                    ilo_D3_gb, ihi_D3_gb, 
     &                    jlo_D3_gb, jhi_D3_gb, 
     &                    D2,
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb, 
     &                    ilo_fb, ihi_fb, 
     &                    jlo_fb-1, jhi_fb+1,
     &                    order_3, y_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         { begin calculation of phi_y
          if (abs(vel_y(i,j)) .lt. zero_tol) then

c           vel_y == 0
            phi_y(i,j) = zero

          elseif (vel_y(i,j) .gt. 0) then

c           vel_y > 0
            phi_y(i,j) = D1(i,j)
            if (abs(D2(i,j-1)).lt.abs(D2(i,j))) then
              phi_y(i,j) = phi_y(i,j) + half*D2(i,j-1) 
              if (abs(D3(i,j-1)).lt.abs(D3(i,j))) then
                phi_y(i,j) = phi_y(i,j) + third*D3(i,j-1)
              else
                phi_y(i,j) = phi_y(i,j) + third*D3(i,j)
              endif
            else
              phi_y(i,j) = phi_y(i,j) + half*D2(i,j) 
              if (abs(D3(i,j)).lt.abs(D3(i,j+1))) then
                phi_y(i,j) = phi_y(i,j) - sixth*D3(i,j)
              else
                phi_y(i,j) = phi_y(i,j) - sixth*D3(i,j+1)
              endif
            endif

          else

c           vel_y < 0
            phi_y(i,j) = D1(i,j+1)
  
            if (abs(D2(i,j)).lt.abs(D2(i,j+1))) then
              phi_y(i,j) = phi_y(i,j) - half*D2(i,j) 
              if (abs(D3(i,j)).lt.abs(D3(i,j+1))) then
                phi_y(i,j) = phi_y(i,j) - sixth*D3(i,j)
              else
                phi_y(i,j) = phi_y(i,j) - sixth*D3(i,j+1)
              endif
            else
              phi_y(i,j) = phi_y(i,j) - half*D2(i,j+1) 
              if (abs(D3(i,j+1)).lt.abs(D3(i,j+2))) then
                phi_y(i,j) = phi_y(i,j) + third*D3(i,j+1)
              else
                phi_y(i,j) = phi_y(i,j) + third*D3(i,j+2)
              endif
            endif
  
          endif

c         divide phi_y by dy
          phi_y(i,j) = phi_y(i,j)*inv_dy

c         } end calculation of phi_y

        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dUpwindHJWENO5() computes the fifth-order Hamilton-Jacobi WENO
c  upwind approximation to the gradient of phi.  
c
c  Arguments:
c    phi_* (out):  components of grad(phi)
c    phi (in):     phi
c    vel_* (in):   components of the velocity
c    D1 (in):      scratch space for holding undivided first-differences
c    dx, dy (in):  grid spacing
c    *_gb (in):    index range for ghostbox
c    *_fb (in):    index range for fillbox
c
c***********************************************************************
      subroutine lsm2dUpwindHJWENO5(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb,
     &  vel_x, vel_y,
     &  ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _vel_gb refers to ghostbox for velocity data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb, jlo_phi_gb, jhi_phi_gb
      integer ilo_vel_gb, ihi_vel_gb, jlo_vel_gb, jhi_vel_gb
      integer ilo_D1_gb, ihi_D1_gb, jlo_D1_gb, jhi_D1_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &           jlo_vel_gb:jhi_vel_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb)
      real dx, dy
      real inv_dx, inv_dy

c     variables for WENO calculation 
      real v1,v2,v3,v4,v5
      real S1,S2,S3
      real a1,a2,a3, inv_sum_a
      real phi_x_1,phi_x_2,phi_x_3
      real phi_y_1,phi_y_2,phi_y_3
      real tiny_nonzero_number
      parameter (tiny_nonzero_number=1.d-99)
      real eps
      real one_third, seven_sixths, eleven_sixths
      real one_sixth, five_sixths
      real thirteen_twelfths, one_fourth
      parameter (one_third=1.d0/3.d0)
      parameter (seven_sixths=7.d0/6.d0)
      parameter (eleven_sixths=11.d0/6.d0) 
      parameter (one_sixth=1.d0/6.d0)
      parameter (five_sixths=5.d0/6.d0)
      parameter (thirteen_twelfths=13.d0/12.d0)
      parameter (one_fourth=0.25d0)

      integer i,j
      real zero
      parameter (zero=0.0d0)
      real zero_tol
      parameter (zero_tol=1.d-11)
      integer order_1
      parameter (order_1=1)
      integer x_dir, y_dir
      parameter (x_dir=1, y_dir=2)


c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute upwind phi_x 
c----------------------------------------------------

c     compute first undivided differences in x-direction
      call lsm2dComputeDn(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb-2, ihi_fb+2,
     &                    jlo_fb, jhi_fb,
     &                    order_1, x_dir)
c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         { begin upwind cases in x-direction
          if (abs(vel_x(i,j)) .lt. zero_tol) then
            phi_x(i,j) = zero
          else
            if (vel_x(i,j) .gt. 0) then
  
c             extract v1,v2,v3,v4,v5 from D1
              v1 = D1(i-2,j)*inv_dx
              v2 = D1(i-1,j)*inv_dx
              v3 = D1(i,j)*inv_dx
              v4 = D1(i+1,j)*inv_dx
              v5 = D1(i+2,j)*inv_dx
  
            else 
  
c             extract v1,v2,v3,v4,v5 from D1
              v1 = D1(i+3,j)*inv_dx
              v2 = D1(i+2,j)*inv_dx
              v3 = D1(i+1,j)*inv_dx
              v4 = D1(i,j)*inv_dx
              v5 = D1(i-1,j)*inv_dx
  
            endif
  
c           WENO5 algorithm for current grid point using appropriate
c           upwind values for v1,...,v5
  
c           compute eps for current grid point
            eps = 1e-6*max(v1*v1,v2*v2,v3*v3,v4*v4,v5*v5)
     &          + tiny_nonzero_number

c           compute the phi_x_1, phi_x_2, phi_x_3
            phi_x_1 = one_third*v1 - seven_sixths*v2 + eleven_sixths*v3
            phi_x_2 = -one_sixth*v2 + five_sixths*v3 + one_third*v4
            phi_x_3 = one_third*v3 + five_sixths*v4 - one_sixth*v5
  
c           compute the smoothness measures
            S1 = thirteen_twelfths*(v1-2.d0*v2+v3)**2
     &         + one_fourth*(v1-4.d0*v2+3.d0*v3)**2
            S2 = thirteen_twelfths*(v2-2.d0*v3+v4)**2
     &         + one_fourth*(v2-v4)**2
            S3 = thirteen_twelfths*(v3-2.d0*v4+v5)**2
     &         + one_fourth*(3.d0*v3-4.d0*v4+v5)**2

c           compute normalized weights
            a1 = 0.1d0/(S1+eps)**2
            a2 = 0.6d0/(S2+eps)**2
            a3 = 0.3d0/(S3+eps)**2
            inv_sum_a = 1.0d0 / (a1 + a2 + a3)
            a1 = a1*inv_sum_a
            a2 = a2*inv_sum_a
            a3 = a3*inv_sum_a
  
c           compute phi_x 
            phi_x(i,j) = a1*phi_x_1 + a2*phi_x_2 + a3*phi_x_3
  
          endif
c         } end upwind cases in x-direction

        enddo
      enddo
c     } end loop over grid 


c----------------------------------------------------
c    compute upwind phi_y
c----------------------------------------------------

c     compute first undivided differences in y-direction
      call lsm2dComputeDn(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    ilo_fb, ihi_fb, 
     &                    jlo_fb-2, jhi_fb+2, 
     &                    order_1, y_dir)

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         { begin upwind cases in y-direction
          if (abs(vel_y(i,j)) .lt. zero_tol) then
            phi_y(i,j) = zero
          else
            if (vel_y(i,j) .gt. 0) then
  
c             extract v1,v2,v3,v4,v5 from D1
              v1 = D1(i,j-2)*inv_dy
              v2 = D1(i,j-1)*inv_dy
              v3 = D1(i,j)*inv_dy
              v4 = D1(i,j+1)*inv_dy
              v5 = D1(i,j+2)*inv_dy
  
            else 
  
c             extract v1,v2,v3,v4,v5 from D1
              v1 = D1(i,j+3)*inv_dy
              v2 = D1(i,j+2)*inv_dy
              v3 = D1(i,j+1)*inv_dy
              v4 = D1(i,j)*inv_dy
              v5 = D1(i,j-1)*inv_dy
  
            endif
  
c           WENO5 algorithm for current grid point using appropriate
c           upwind values for v1,...,v5
  
c           compute eps for current grid point
            eps = 1e-6*max(v1*v1,v2*v2,v3*v3,v4*v4,v5*v5)
     &          + tiny_nonzero_number

c           compute the phi_y_1, phi_y_2, phi_y_3
            phi_y_1 = one_third*v1 - seven_sixths*v2 + eleven_sixths*v3
            phi_y_2 = -one_sixth*v2 + five_sixths*v3 + one_third*v4
            phi_y_3 = one_third*v3 + five_sixths*v4 - one_sixth*v5
  
c           compute the smoothness measures
            S1 = thirteen_twelfths*(v1-2.d0*v2+v3)**2
     &         + one_fourth*(v1-4.d0*v2+3.d0*v3)**2
            S2 = thirteen_twelfths*(v2-2.d0*v3+v4)**2
     &         + one_fourth*(v2-v4)**2
            S3 = thirteen_twelfths*(v3-2.d0*v4+v5)**2
     &         + one_fourth*(3.d0*v3-4.d0*v4+v5)**2

c           compute normalized weights
            a1 = 0.1d0/(S1+eps)**2
            a2 = 0.6d0/(S2+eps)**2
            a3 = 0.3d0/(S3+eps)**2
            inv_sum_a = 1.0d0 / (a1 + a2 + a3)
            a1 = a1*inv_sum_a
            a2 = a2*inv_sum_a
            a3 = a3*inv_sum_a
  
c           compute phi_y 
            phi_y(i,j) = a1*phi_y_1 + a2*phi_y_2 + a3*phi_y_3
  
          endif
c         } end upwind cases in y-direction

        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dCentralGradOrder2() computes the second-order, central, 
c  finite difference approximation to the gradient of phi.
c
c  Arguments:
c    phi_* (out):  components of grad(phi) 
c    phi (in):     phi
c    dx, dy (in):  grid spacing
c    *_gb (in):    index range for ghostbox
c    *_fb (in):    index range for fillbox
c
c***********************************************************************
      subroutine lsm2dCentralGradOrder2(
     &  phi_x, phi_y, 
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      real dx_factor, dy_factor
      integer i,j

c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          phi_x(i,j) = (phi(i+1,j) - phi(i-1,j))*dx_factor
          phi_y(i,j) = (phi(i,j+1) - phi(i,j-1))*dy_factor
 
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dCentralGradOrder4() computes the second-order, central,
c  finite difference approximation to the gradient of phi.
c
c  Arguments:
c    phi_* (out):  components of grad(phi) 
c    phi (in):     phi
c    dx, dy (in):  grid spacing
c    *_gb (in):    index range for ghostbox
c    *_fb (in):    index range for fillbox
c
c***********************************************************************
      subroutine lsm2dCentralGradOrder4(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      real dx_factor, dy_factor
      integer i,j
      real eight
      parameter (eight = 8.0d0)

c     compute denominator values
      dx_factor = 0.0833333333333333333333d0/dx
      dy_factor = 0.0833333333333333333333d0/dy

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          phi_x(i,j) = ( -phi(i+2,j) + eight*phi(i+1,j) 
     &                   +phi(i-2,j) - eight*phi(i-1,j) )
     &               * dx_factor
          phi_y(i,j) = ( -phi(i,j+2) + eight*phi(i,j+1) 
     &                   +phi(i,j-2) - eight*phi(i,j-1) )
     &               * dy_factor
   
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dLaplacianOrder2() computes the second-order, central,
c  finite difference approximation to the Laplacian of phi.
c
c  Arguments:
c    laplacian_phi (out):  Laplacian of phi
c    phi (in):             phi
c    dx (in):              grid spacing
c    *_gb (in):            index range for ghostbox
c    *_fb (in):            index range for fillbox
c
c***********************************************************************
      subroutine lsm2dLaplacianOrder2(
     &  laplacian_phi,
     &  ilo_laplacian_phi_gb, ihi_laplacian_phi_gb, 
     &  jlo_laplacian_phi_gb, jhi_laplacian_phi_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb, 
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _laplacian_phi_gb refers to ghostbox for laplacian_phi data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_laplacian_phi_gb, ihi_laplacian_phi_gb
      integer jlo_laplacian_phi_gb, jhi_laplacian_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real laplacian_phi(ilo_laplacian_phi_gb:ihi_laplacian_phi_gb,
     &                   jlo_laplacian_phi_gb:jhi_laplacian_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      integer i, j
      real inv_dx_sq
      real inv_dy_sq

c     compute denominator values
      inv_dx_sq = 1.0d0/dx/dx
      inv_dy_sq = 1.0d0/dy/dy

c     { begin loop over grid 
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          laplacian_phi(i,j) = 
     &      inv_dx_sq*( phi(i+1,j) - 2.0d0*phi(i,j) + phi(i-1,j) ) 
     &    + inv_dy_sq*( phi(i,j+1) - 2.0d0*phi(i,j) + phi(i,j-1) ) 
   
        enddo
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dPhiUpwindGradF() computes the "phi-upwind" gradient of a 
c  function, F, using the following "upwinding" scheme to compute 
c  the normal:
c
c    if phi > 0:  upwind direction is direction where phi is smaller
c
c    if phi < 0:  upwind direction is direction where phi is larger
c
c  Arguments:
c    F_* (out):       components of phi-upwinded grad(F)
c    F_*_plus (in):   components of grad(F) in plus direction
c    F_*_minus (in):  components of grad(F) in minus direction
c    phi (in):        level set function
c    *_gb (in):       index range for ghostbox
c    *_fb (in):       index range for fillbox
c
c  NOTES:
c   - phi is REQUIRED to have at least one ghost cell in each 
c     coordinate direction for upwinding
c
c***********************************************************************
      subroutine lsm2dPhiUpwindGradF(
     &  F_x, F_y,
     &  ilo_grad_F_gb, ihi_grad_F_gb,
     &  jlo_grad_F_gb, jhi_grad_F_gb,
     &  F_x_plus, F_y_plus, 
     &  ilo_grad_F_plus_gb, ihi_grad_F_plus_gb,
     &  jlo_grad_F_plus_gb, jhi_grad_F_plus_gb,
     &  F_x_minus, F_y_minus, 
     &  ilo_grad_F_minus_gb, ihi_grad_F_minus_gb,
     &  jlo_grad_F_minus_gb, jhi_grad_F_minus_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb,
     &  jlo_phi_gb, jhi_phi_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fill-box

      integer ilo_grad_F_gb, ihi_grad_F_gb
      integer jlo_grad_F_gb, jhi_grad_F_gb
      integer ilo_grad_F_plus_gb, ihi_grad_F_plus_gb
      integer jlo_grad_F_plus_gb, jhi_grad_F_plus_gb
      integer ilo_grad_F_minus_gb, ihi_grad_F_minus_gb
      integer jlo_grad_F_minus_gb, jhi_grad_F_minus_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real F_x(ilo_grad_F_gb:ihi_grad_F_gb,
     &         jlo_grad_F_gb:jhi_grad_F_gb)
      real F_y(ilo_grad_F_gb:ihi_grad_F_gb,
     &         jlo_grad_F_gb:jhi_grad_F_gb)
      real F_x_plus(ilo_grad_F_plus_gb:ihi_grad_F_plus_gb,
     &              jlo_grad_F_plus_gb:jhi_grad_F_plus_gb)
      real F_y_plus(ilo_grad_F_plus_gb:ihi_grad_F_plus_gb,
     &              jlo_grad_F_plus_gb:jhi_grad_F_plus_gb)
      real F_x_minus(ilo_grad_F_minus_gb:ihi_grad_F_minus_gb,
     &               jlo_grad_F_minus_gb:jhi_grad_F_minus_gb)
      real F_y_minus(ilo_grad_F_minus_gb:ihi_grad_F_minus_gb,
     &               jlo_grad_F_minus_gb:jhi_grad_F_minus_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real phi_cur
      real phi_neighbor_plus, phi_neighbor_minus
      integer i,j
      real zero
      parameter (zero=0.0d0)

c     compute "phi-upwind" derivatives
c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

c         cache current phi
          phi_cur = phi(i,j)

c         { begin computation of "upwind" derivative
          if (phi_cur .gt. 0) then

c           compute "upwind" derivative in x-direction

            phi_neighbor_minus = phi(i-1,j)
            phi_neighbor_plus = phi(i+1,j)
            if (phi_neighbor_minus .le. phi_cur) then
              if (phi_neighbor_plus .lt. phi_neighbor_minus) then
                F_x(i,j) = F_x_plus(i,j) 
              else
                F_x(i,j) = F_x_minus(i,j) 
              endif
            elseif (phi_neighbor_plus .le. phi_cur) then
              F_x(i,j) = F_x_plus(i,j) 
            else
              F_x(i,j) = zero
            endif

c           compute "upwind" derivative in y-direction

            phi_neighbor_minus = phi(i,j-1)
            phi_neighbor_plus = phi(i,j+1)
            if (phi_neighbor_minus .le. phi_cur) then
              if (phi_neighbor_plus .lt. phi_neighbor_minus) then
                F_y(i,j) = F_y_plus(i,j) 
              else
                F_y(i,j) = F_y_minus(i,j) 
              endif
            elseif (phi_neighbor_plus .le. phi_cur) then
              F_y(i,j) = F_y_plus(i,j) 
            else
              F_y(i,j) = zero
            endif

          elseif (phi_cur .lt. 0) then

c           compute "upwind" derivative in x-direction

            phi_neighbor_minus = phi(i-1,j)
            phi_neighbor_plus = phi(i+1,j)
            if (phi_neighbor_minus .ge. phi_cur) then
            if (phi_neighbor_plus .gt. phi_neighbor_minus) then
                F_x(i,j) = F_x_plus(i,j) 
              else
                F_x(i,j) = F_x_minus(i,j) 
              endif
            elseif (phi_neighbor_plus .ge. phi_cur) then
              F_x(i,j) = F_x_plus(i,j) 
            else
              F_x(i,j) = zero
            endif

c           compute "upwind" derivative in y-direction

            phi_neighbor_minus = phi(i,j-1)
            phi_neighbor_plus = phi(i,j+1)
            if (phi_neighbor_minus .ge. phi_cur) then
              if (phi_neighbor_plus .gt. phi_neighbor_minus) then
                F_y(i,j) = F_y_plus(i,j) 
              else
                F_y(i,j) = F_y_minus(i,j) 
              endif
            elseif (phi_neighbor_plus .ge. phi_cur) then
              F_y(i,j) = F_y_plus(i,j) 
            else
              F_y(i,j) = zero
            endif

c         } end computation of "upwind" derivative
          endif

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dAverageGradPhi() computes the average of the plus and minus 
c  derivatives:
c
c    phi_* = (phi_*_plus + phi_*_minus) / 2
c
c  Arguments:
c    phi_* (out):       components of average grad(phi)
c    phi_*_plus (in):   components of grad(phi) in plus direction
c    phi_*_minus (in):  components of grad(phi) in minus direction
c    *_gb (in):         index range for ghostbox
c    *_fb (in):         index range for fillbox
c
c***********************************************************************
      subroutine lsm2dAverageGradPhi(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb,
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb,
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb,
     &  ilo_fb, ihi_fb,
     &  jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fill-box

      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb
      integer jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer ilo_fb, ihi_fb
      integer jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_x_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                 jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      integer i,j
      real half
      parameter (half=0.5d0)

c     compute "phi-upwind" derivatives
c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          phi_x(i,j) = half * ( phi_x_plus(i,j) 
     &                        + phi_x_minus(i,j) )
          phi_y(i,j) = half * ( phi_y_plus(i,j) 
     &                        + phi_y_minus(i,j) )

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm2dGradientMagnitude() computes magnitude of the gradient of phi.
c
c  Arguments:
c    phi_* (in):          components of grad(phi)
c    grad_phi_mag (out):  gradient magnitude
c    *_gb (in):           index range for ghostbox
c    *_fb (in):           index range for fillbox
c
c***********************************************************************
      subroutine lsm2dGradientMagnitude(
     &  phi_x, phi_y, grad_phi_mag,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real grad_phi_mag(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                  jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_x_sq, phi_y_sq, zero_tol, zero, tmp
      parameter (zero = 0.d0)
      parameter (zero_tol=1.d-11)
      integer i,j

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb
          
          phi_x_sq = phi_x(i,j)*phi_x(i,j)
          phi_y_sq = phi_y(i,j)*phi_y(i,j)
          
          if( phi_x_sq .lt. zero_tol) then
              phi_x_sq = zero
          endif
          
          if( phi_y_sq .lt. zero_tol) then
               phi_y_sq = zero
          endif
          
          tmp = sqrt( phi_x_sq + phi_y_sq )
          if( tmp .lt. zero_tol ) then
            grad_phi_mag(i,j) = zero
          else
            grad_phi_mag(i,j) = tmp
          endif

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
      subroutine lsm2dDivergenceCentral(
     &  divF,
     &  ilo_divf_gb, ihi_divf_gb,
     &  jlo_divf_gb, jhi_divf_gb,
     &  FX, FY,
     &  ilo_gb, ihi_gb,jlo_gb, jhi_gb,
     &  ilo_fb, ihi_fb, jlo_fb, jhi_fb,
     &  dx, dy)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox for data
c     _fb refers to fill-box for data
      integer ilo_divf_gb, ihi_divf_gb
      integer jlo_divf_gb, jhi_divf_gb
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      integer ilo_fb, ihi_fb, jlo_fb, jhi_fb
      real divF(ilo_divf_gb:ihi_divf_gb,
     &          jlo_divf_gb:jhi_divf_gb)
      real FX(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real FY(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real dx, dy
      real dx_factor, dy_factor
      integer i,j

c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     { begin loop over grid
      do j=jlo_fb,jhi_fb
        do i=ilo_fb,ihi_fb

          divF(i,j) = (FX(i+1,j) - FX(i-1,j))*dx_factor +
     &                (FY(i,j+1) - FY(i,j-1))*dy_factor

        enddo
      enddo
c     } end loop over grid

      return
      end
c } end subroutine
c***********************************************************************
