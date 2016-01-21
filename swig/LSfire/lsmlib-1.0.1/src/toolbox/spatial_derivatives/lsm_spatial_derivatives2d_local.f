c***********************************************************************
c
c  File:        lsm_spatial_derivatives2d_local.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 routines for computing 2D ENO/WENO spatial 
c               derivatives on narrow-bands
c
c***********************************************************************

c***********************************************************************
c
c  lsm2dComputeDnLOCAL() computes the n-th undivided differences in the 
c  specified direction given the (n-1)-th undivided differences.  The 
c  subroutine assumes that valid data for the (n-1)-th undivided 
c  differences is only available n/2 or (n+1)/2 (depending on the parity 
c  of n) cells in from the boundary of the ghost-cell box.  The 
c  undivided differences in cells with insufficient data is set to a 
c  large number. 
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    Dn (out):           n-th undivided differences 
c    Dn_minus_one (in):  (n-1)-th undivided differences 
c    n (in):             order of undivided differences to compute
c    *_gb (in):          index range for ghostbox
c    index_[xy](in):     [xy] coordinates of local (narrow band) points
c    n*_index(in):       index range of points to loop over in index_*
c    narrow_band(in):    array that marks voxels outside desired fillbox
c    mark_fb(in):        upper limit narrow band value for voxels in 
c                        fillbox
c
c  NOTES:
c   - The index ranges for all ghostboxes and the fillbox should 
c     correspond to the index range for cell-centered data.
c   - The ghostbox for Dn_minus_one MUST be at least one ghostcell width
c     larger than the fillbox.
c
c***********************************************************************
      subroutine lsm2dComputeDnLOCAL(
     &  Dn,
     &  ilo_Dn_gb, ihi_Dn_gb, 
     &  jlo_Dn_gb, jhi_Dn_gb, 
     &  Dn_minus_one,
     &  ilo_Dn_minus_one_gb, ihi_Dn_minus_one_gb, 
     &  jlo_Dn_minus_one_gb, jhi_Dn_minus_one_gb, 
     &  n,
     &  dir,
     &  index_x, index_y,
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_Dn_gb, ihi_Dn_gb
      integer jlo_Dn_gb, jhi_Dn_gb
      integer ilo_Dn_minus_one_gb, ihi_Dn_minus_one_gb
      integer jlo_Dn_minus_one_gb, jhi_Dn_minus_one_gb
      real Dn(ilo_Dn_gb:ihi_Dn_gb,
     &        jlo_Dn_gb:jhi_Dn_gb)
      real Dn_minus_one(ilo_Dn_minus_one_gb:ihi_Dn_minus_one_gb,
     &                  jlo_Dn_minus_one_gb:jhi_Dn_minus_one_gb)
      integer n
      integer dir
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb     
            
      integer i,j,l
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

c     loop over indexed points only {
      do l= nlo_index, nhi_index      
        i = index_x(l) 
        j = index_y(l)
        if( narrow_band(i,j) .le. mark_fb ) then
              Dn(i,j) = sign_multiplier
     &          * ( Dn_minus_one(i,j)
     &            - Dn_minus_one(i-offset(1),j-offset(2)) )          
        else
          Dn(i,j) = big
        endif
      
      enddo
c     }  end loop over indexed points 

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dHJENO1LOCAL(
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, 
     &  jlo_D1_gb, jhi_D1_gb,
     &  dx, dy,
     &  index_x, index_y,
     &  nlo_index0, nhi_index0,
     &  nlo_index1, nhi_index1,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb,
     &  mark_D1)
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
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_D1_gb, ihi_D1_gb
      integer jlo_D1_gb, jhi_D1_gb
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
      integer nlo_index0, nhi_index0
      integer nlo_index1, nhi_index1
      integer index_x(nlo_index0:nhi_index1)
      integer index_y(nlo_index0:nhi_index1)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb, mark_D1     
      
      real inv_dx, inv_dy
      integer i,j,l
      integer order
      parameter (order=1)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)


c     compute inv_dx, inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    order, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1)

c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

            phi_x_plus(i,j) = D1(i+1,j)*inv_dx
            phi_x_minus(i,j) = D1(i,j)*inv_dx
   
        endif
      enddo
c     } end loop over narrow band points

c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    order, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1)

c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

            phi_y_minus(i,j) = D1(i,j)*inv_dy
            phi_y_plus(i,j)  = D1(i,j+1)*inv_dy
   
        endif
      enddo
c     } end loop over narrow band points 

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dHJENO2LOCAL(
     &  phi_x_plus, phi_y_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, 
     &  jlo_D1_gb, jhi_D1_gb,
     &  D2,
     &  ilo_D2_gb, ihi_D2_gb, 
     &  jlo_D2_gb, jhi_D2_gb,
     &  dx, dy,
     &  index_x,
     &  index_y, 
     &  nlo_index0, nhi_index0,
     &  nlo_index1, nhi_index1,
     &  nlo_index2, nhi_index2,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb,
     &  mark_D1,
     &  mark_D2)
     
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
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer ilo_D1_gb, ihi_D1_gb
      integer jlo_D1_gb, jhi_D1_gb
      integer ilo_D2_gb, ihi_D2_gb
      integer jlo_D2_gb, jhi_D2_gb
      integer nlo_index0, nhi_index0
      integer nlo_index1, nhi_index1
      integer nlo_index2, nhi_index2
      integer index_x(nlo_index0:nhi_index2)
      integer index_y(nlo_index0:nhi_index2)
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
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_D2
      integer*1 mark_D1
      integer*1 mark_fb
      
      real inv_dx, inv_dy
      integer i,j,l
      real half
      parameter (half=0.5d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)

c     compute inv_dx, inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    order_1, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1) 

c     compute second undivided differences x-direction
      call lsm2dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    order_2, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D2) 
c    loop over narrow band level 0 points {
      do l=nlo_index0, nhi_index0   
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

c             phi_x_plus
              if (abs(D2(i,j)).lt.abs(D2(i+1,j))) then
                phi_x_plus(i,j) = (D1(i+1,j) 
     &                          - half*D2(i,j))*inv_dx
              else
                phi_x_plus(i,j) = (D1(i+1,j) 
     &                          - half*D2(i+1,j))*inv_dx
              endif

c             phi_x_minus
              if (abs(D2(i-1,j)).lt.abs(D2(i,j))) then
                phi_x_minus(i,j) = (D1(i,j) 
     &                           + half*D2(i-1,j))*inv_dx
              else
                phi_x_minus(i,j) = (D1(i,j) 
     &                           + half*D2(i,j))*inv_dx
              endif
        endif      
      enddo
c     } end loop over indexed points


c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    order_1, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1) 
     
c     compute second undivided differences in y-direction
      call lsm2dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    order_2, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D2)
      
c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then
c             phi_y_plus
              if (abs(D2(i,j)).lt.abs(D2(i,j+1))) then
                phi_y_plus(i,j) = (D1(i,j+1) 
     &                          - half*D2(i,j))*inv_dy
              else
                phi_y_plus(i,j) = (D1(i,j+1) 
     &                          - half*D2(i,j+1))*inv_dy
              endif

c             phi_y_minus
              if (abs(D2(i,j-1)).lt.abs(D2(i,j))) then
                phi_y_minus(i,j) = (D1(i,j) 
     &                           + half*D2(i,j-1))*inv_dy
              else
                phi_y_minus(i,j) = (D1(i,j) 
     &                           + half*D2(i,j))*inv_dy
              endif
        endif      
      enddo
c     } end loop over narrow band points

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dHJENO3LOCAL(
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
     &  dx, dy, 
     &  nlo_index0, nhi_index0,
     &  nlo_index1, nhi_index1,
     &  nlo_index2, nhi_index2,
     &  nlo_index3, nhi_index3,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb,
     &  mark_D1,
     &  mark_D2,
     &  mark_D3)
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
      real phi_x_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &                    jlo_D1_gb:jhi_D1_gb)
      real D2(ilo_D2_gb:ihi_D2_gb,
     &                    jlo_D2_gb:jhi_D2_gb)
      real D3(ilo_D3_gb:ihi_D3_gb,
     &                    jlo_D3_gb:jhi_D3_gb)
      real dx, dy
      real inv_dx, inv_dy
      integer i,j,l
      real zero, half, third, sixth
      parameter (zero=0.0d0, half=0.5d0, third=1.d0/3.d0)
      parameter (sixth=1.d0/6.d0)
      real zero_tol
      parameter (zero_tol=1.d-11)
      integer order_1, order_2, order_3
      parameter (order_1=1,order_2=2,order_3=3)
      integer x_dir, y_dir
      parameter (x_dir=1, y_dir=2)

      integer nlo_index0, nhi_index0
      integer nlo_index1, nhi_index1
      integer nlo_index2, nhi_index2
      integer nlo_index3, nhi_index3
      integer index_x(nlo_index0:nhi_index3)
      integer index_y(nlo_index0:nhi_index3)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_D1, mark_D2, mark_D3
      integer*1 mark_fb

c     compute inv_dx and inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------

c     compute first undivided differences in x-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    order_1, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index3,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1)

c     compute second undivided differences x-direction
      call lsm2dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    order_2, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D2) 

c     compute third undivided differences in x-direction
      call lsm2dComputeDnLOCAL(D3, 
     &                    ilo_D3_gb, ihi_D3_gb, 
     &                    jlo_D3_gb, jhi_D3_gb,
     &                    D2,
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    order_3, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D3) 

c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

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

        endif
      enddo
c     } end loop over grid 


c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    order_1, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index3,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1) 

c     compute second undivided differences x-direction
      call lsm2dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    order_2, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D2) 

c     compute third undivided differences in x-direction
      call lsm2dComputeDnLOCAL(D3, 
     &                    ilo_D3_gb, ihi_D3_gb, 
     &                    jlo_D3_gb, jhi_D3_gb,
     &                    D2,
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    order_3, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D3) 

c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

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

        endif
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dHJWENO5LOCAL(
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
     &  dx, dy,
     &  index_x,
     &  index_y, 
     &  nlo_index0, nhi_index0,
     &  nlo_index1, nhi_index1,
     &  nlo_index2, nhi_index2,
     &  nlo_index3, nhi_index3,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb,
     &  mark_D1)
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
      real phi_x_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_y_plus(
     &                    ilo_grad_phi_plus_gb:ihi_grad_phi_plus_gb,
     &                    jlo_grad_phi_plus_gb:jhi_grad_phi_plus_gb)
      real phi_x_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi_y_minus(
     &                    ilo_grad_phi_minus_gb:ihi_grad_phi_minus_gb,
     &                    jlo_grad_phi_minus_gb:jhi_grad_phi_minus_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &                    jlo_D1_gb:jhi_D1_gb)
      real dx, dy
      real inv_dx, inv_dy

      integer nlo_index0, nhi_index0
      integer nlo_index1, nhi_index1
      integer nlo_index2, nhi_index2
      integer nlo_index3, nhi_index3
      integer index_x(nlo_index0:nhi_index3)
      integer index_y(nlo_index0:nhi_index3)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_D1
      integer*1 mark_fb
      
c     variables for WENO calculation 
      real v1,v2,v3,v4,v5
      real S1,S2,S3
      real a1,a2,a3, inv_sum_a
      real phi_x_1,phi_x_2,phi_x_3
      real phi_y_1,phi_y_2,phi_y_3
      real tiny_nonzero_number
      parameter (tiny_nonzero_number=1.d-35)
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

      integer i,j,l
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
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    order_1, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index3,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1) 

c    loop over narrow band level 0 points {
      do l=nlo_index0, nhi_index0   
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

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

        endif
      enddo
c     } end loop over narrow band level 0 points 


c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    order_1, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index3,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1) 

c    loop over narrow band level 0 points {
      do l=nlo_index0, nhi_index0   
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then
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

        endif
      enddo
c     } end loop over narrow band points 

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dUpwindHJENO2LOCAL(
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
     &  dx, dy,
     &  index_x,
     &  index_y, 
     &  nlo_index0, nhi_index0,
     &  nlo_index1, nhi_index1,
     &  nlo_index2, nhi_index2,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb,
     &  mark_D1,
     &  mark_D2)
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
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                       jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                       jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &                     jlo_phi_gb:jhi_phi_gb)
      real vel_x(ilo_vel_gb:ihi_vel_gb,
     &                       jlo_vel_gb:jhi_vel_gb)
      real vel_y(ilo_vel_gb:ihi_vel_gb,
     &                       jlo_vel_gb:jhi_vel_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &                    jlo_D1_gb:jhi_D1_gb)
      real D2(ilo_D2_gb:ihi_D2_gb,
     &                    jlo_D2_gb:jhi_D2_gb)
      real dx, dy
      
      integer nlo_index0, nhi_index0
      integer nlo_index1, nhi_index1
      integer nlo_index2, nhi_index2
      integer index_x(nlo_index0:nhi_index2)
      integer index_y(nlo_index0:nhi_index2)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_D2
      integer*1 mark_D1
      integer*1 mark_fb
      
      real inv_dx, inv_dy
      integer i,j,l
      real zero_tol, zero
      parameter (zero_tol=1.d-11, zero = 0.d0)
      real half
      parameter (half=0.5d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir
      parameter (x_dir=1,y_dir=2)

c     compute inv_dx, inv_dy
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    order_1, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1) 

c     compute second undivided differences x-direction
      call lsm2dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    order_2, x_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D2) 
c    loop over narrow band level 0 points {
      do l=nlo_index0, nhi_index0   
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

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

        endif
      enddo
c     } end loop over indexed points


c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm2dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    order_1, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D1) 
     
c     compute second undivided differences in y-direction
      call lsm2dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    order_2, y_dir,
     &                    index_x, index_y,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    mark_D2)
      
c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

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

        endif
      enddo
c     } end loop over narrow band points

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dCentralGradOrder2LOCAL(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  dx, dy,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
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
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb

c     local variables      
      integer i,j,l
      real dx_factor, dy_factor

c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then
         
          phi_x(i,j) = (phi(i+1,j) - phi(i-1,j))*dx_factor
          phi_y(i,j) = (phi(i,j+1) - phi(i,j-1))*dy_factor

        endif  
      enddo
c     } end loop over indexed points
      
      
      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dCentralGradOrder4LOCAL(
     &  phi_x, phi_y,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  dx, dy,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
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
      real phi_x(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &           jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb
      
      integer i,j,l
      real dx_factor, dy_factor
      real eight
      parameter (eight = 8.0d0)

c     compute denominator values
      dx_factor = 0.0833333333333333333332d0/dx
      dy_factor = 0.0833333333333333333332d0/dy

c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

            phi_x(i,j) = ( -phi(i+2,j) + eight*phi(i+1,j) 
     &                     +phi(i-2,j) - eight*phi(i-1,j) )
     &                   * dx_factor
            phi_y(i,j) = ( -phi(i,j+2) + eight*phi(i,j+1) 
     &                     +phi(i,j-2) - eight*phi(i,j-1) )
     &                   * dy_factor
          
        endif
      enddo
c     } end loop over indexed points 

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
      subroutine lsm2dLaplacianOrder2LOCAL(
     &  laplacian_phi,
     &  ilo_laplacian_phi_gb, ihi_laplacian_phi_gb, 
     &  jlo_laplacian_phi_gb, jhi_laplacian_phi_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb, 
     &  dx, dy,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
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
      real laplacian_phi(ilo_laplacian_phi_gb:ihi_laplacian_phi_gb,
     &                   jlo_laplacian_phi_gb:jhi_laplacian_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb
            
      integer i, j, l
      real inv_dx_sq
      real inv_dy_sq

c     compute denominator values
      inv_dx_sq = 1.0d0/dx/dx
      inv_dy_sq = 1.0d0/dy/dy

c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

            laplacian_phi(i,j) = 
     &        inv_dx_sq *
     &        ( phi(i+1,j) - 2.0d0*phi(i,j) + phi(i-1,j) ) 
     &      + inv_dy_sq *
     &        ( phi(i,j+1) - 2.0d0*phi(i,j) + phi(i,j-1) ) 
           
        endif
      enddo
c     } end loop over indexed points 

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dComputeAveGradPhiLOCAL(
     &  grad_phi_ave, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  dx, dy,
     &  index_x,
     &  index_y, 
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      real grad_phi_ave
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb)
      real dx, dy
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb

c     local variables      
      integer i,j,l,count
      real dx_factor, dy_factor
      real phi_x, phi_y
      

c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

      grad_phi_ave = 0.d0
      count = 0
c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then
      
          phi_x = (phi(i+1,j) - phi(i-1,j))*dx_factor
          phi_y = (phi(i,j+1) - phi(i,j-1))*dy_factor

          grad_phi_ave = grad_phi_ave + sqrt(phi_x*phi_x + phi_y*phi_y)
          count = count + 1
        endif
      enddo
c     } end loop over indexed points
 
      if ( count .gt. 0 ) then
        grad_phi_ave = grad_phi_ave / (count)
      endif

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dGradientMagnitudeLocal(
     &  phi_x, phi_y, grad_phi_mag,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  index_x,
     &  index_y,
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      real phi_x(
     &                    ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                    jlo_grad_phi_gb:jhi_grad_phi_gb)
      real phi_y(
     &                    ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                    jlo_grad_phi_gb:jhi_grad_phi_gb)
      real grad_phi_mag(
     &                    ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                    jlo_grad_phi_gb:jhi_grad_phi_gb)
      
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb

c     local variables      
      real phi_x_sq, phi_y_sq, zero_tol, zero, tmp
      parameter (zero_tol=1.d-11, zero = 0.d0)
      integer i,j,l

c     { begin loop over indexed points
       do l=nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

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

        endif
      enddo
c     } end loop over indexed points 

      return
      end
c } end subroutine
c***********************************************************************



c***********************************************************************
      subroutine lsm2dDivergenceCentralLocal(
     &  divF, 
     &  ilo_divf_gb, ihi_divf_gb, 
     &  jlo_divf_gb, jhi_divf_gb,
     &  FX, FY,
     &  ilo_gb, ihi_gb,jlo_gb, jhi_gb,
     &  dx, dy,
     &  index_x,
     &  index_y,
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox for data
c     _fb refers to fill-box for data
      integer ilo_divf_gb, ihi_divf_gb
      integer jlo_divf_gb, jhi_divf_gb
      integer ilo_gb, ihi_gb
      integer jlo_gb, jhi_gb
      real divF(ilo_divf_gb:ihi_divf_gb,
     &          jlo_divf_gb:jhi_divf_gb)
      real FX(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real FY(ilo_gb:ihi_gb,jlo_gb:jhi_gb)
      real dx, dy
      
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb)
      integer*1 mark_fb

c     local variables      
      real dx_factor, dy_factor
      integer i,j,l

c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy

c     { begin loop over indexed points
       do l=nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j) .le. mark_fb ) then

          divF(i,j) = (FX(i+1,j) - FX(i-1,j))*dx_factor +
     &                (FY(i,j+1) - FY(i,j-1))*dy_factor
 
        endif
      enddo
c     } end loop over over indexed points 

      return
      end
c } end subroutine
c***********************************************************************