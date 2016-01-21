c***********************************************************************
c
c  File:        lsm_spatial_derivatives3d_local.f
c  Copyrights:  (c) 2005 The Trustees of Princeton University and Board of
c                   Regents of the University of Texas.  All rights reserved.
c               (c) 2009 Kevin T. Chu.  All rights reserved.
c  Revision:    $Revision$
c  Modified:    $Date$
c  Description: F77 routines for computing 3D ENO/WENO spatial 
c               derivatives on narrow-bands
c
c***********************************************************************

c***********************************************************************
c
c  lsm3dComputeDnLOCAL() computes the n-th undivided differences in the 
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
c    index_[xyz](in):    [xyz] coordinates of local (narrow band) points
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
      subroutine lsm3dComputeDnLOCAL(
     &  Dn,
     &  ilo_Dn_gb, ihi_Dn_gb, 
     &  jlo_Dn_gb, jhi_Dn_gb, 
     &  klo_Dn_gb, khi_Dn_gb,
     &  Dn_minus_one,
     &  ilo_Dn_minus_one_gb, ihi_Dn_minus_one_gb, 
     &  jlo_Dn_minus_one_gb, jhi_Dn_minus_one_gb, 
     &  klo_Dn_minus_one_gb, khi_Dn_minus_one_gb,
     &  n,
     &  dir,
     &  index_x, index_y, index_z, 
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _gb refers to ghostbox 
c     _fb refers to fillbox 
      integer ilo_Dn_gb, ihi_Dn_gb
      integer jlo_Dn_gb, jhi_Dn_gb
      integer klo_Dn_gb, khi_Dn_gb
      integer ilo_Dn_minus_one_gb, ihi_Dn_minus_one_gb
      integer jlo_Dn_minus_one_gb, jhi_Dn_minus_one_gb
      integer klo_Dn_minus_one_gb, khi_Dn_minus_one_gb
      real Dn(ilo_Dn_gb:ihi_Dn_gb,
     &        jlo_Dn_gb:jhi_Dn_gb,
     &        klo_Dn_gb:khi_Dn_gb)
      real Dn_minus_one(ilo_Dn_minus_one_gb:ihi_Dn_minus_one_gb,
     &                  jlo_Dn_minus_one_gb:jhi_Dn_minus_one_gb,
     &                  klo_Dn_minus_one_gb:khi_Dn_minus_one_gb)
      integer n
      integer dir
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb     
            
      integer i,j,k,l
      integer offset(1:3)
      integer fillbox_shift(1:3)
      real sign_multiplier
      real big
      parameter (big=1.d10)
      

c     calculate offsets, fillbox shifts, and sign_multiplier used 
c     when computing undivided differences.
c     NOTE:  even and odd undivided differences are taken in
c            opposite order because of the discrepancy between
c            face- and cell-centered data.  the sign discrepancy 
c            is taken into account by sign_multiplier
      do i=1,3
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
        k = index_z(l)
        if( narrow_band(i,j,k) .le. mark_fb ) then
          Dn(i,j,k) = sign_multiplier
     &      * ( Dn_minus_one(i,j,k)
     &      - Dn_minus_one(i-offset(1),j-offset(2),k-offset(3)))          
        else
          Dn(i,j,k) = big
        endif
      
      enddo
c     }  end loop over indexed points 

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm3dHJENO1LOCAL() computes the forward (plus) and backward (minus)
c  first-order Hamilton-Jacobi ENO approximations to the gradient of 
c  phi.
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    phi_*_plus (out):   components of grad(phi) in plus direction 
c    phi_*_minus (out):  components of grad(phi) in minus direction
c    phi (in):           phi
c    D1 (in):            scratch space for holding undivided first-differences
c    dx, dy, dz (in):    grid spacing
c    *_gb (in):          index range for ghostbox
c    index_*(in):        coordinates of local (narrow band) points
c    n*_index[01](in):   index range of points in index_* that are in
c                        level [01] of the narrow band
c    narrow_band(in):    array that marks voxels outside desired fillbox
c    mark_*(in):         upper limit narrow band value for voxels in 
c                        the appropriate fillbox
c
c  NOTES:
c   - it is assumed that BOTH the plus AND minus derivatives have
c     the same fillbox
c   - index_[xyz] arrays range at minimum from nlo_index0 to nhi_index1
c
c***********************************************************************
      subroutine lsm3dHJENO1LOCAL(
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, 
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, 
     &  jlo_D1_gb, jhi_D1_gb,
     &  klo_D1_gb, khi_D1_gb,
     &  dx, dy, dz,
     &  index_x, index_y, index_z, 
     &  nlo_index0, nhi_index0,
     &  nlo_index1, nhi_index1,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
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
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_D1_gb, ihi_D1_gb
      integer jlo_D1_gb, jhi_D1_gb
      integer klo_D1_gb, khi_D1_gb
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
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb,
     &        klo_D1_gb:khi_D1_gb)
      real dx, dy, dz
      integer nlo_index0, nhi_index0
      integer nlo_index1, nhi_index1
      integer index_x(nlo_index0:nhi_index1)
      integer index_y(nlo_index0:nhi_index1)
      integer index_z(nlo_index0:nhi_index1)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb, mark_D1     
      
      real inv_dx, inv_dy, inv_dz
      integer i,j,k,l
      integer order
      parameter (order=1)
      integer x_dir, y_dir, z_dir
      parameter (x_dir=1,y_dir=2,z_dir=3)


c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      inv_dz = 1.0d0/dz

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------
c     compute first undivided differences in x-direction
      call lsm3dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    klo_D1_gb, khi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    klo_phi_gb, khi_phi_gb, 
     &                    order, x_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D1)

c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)
        k = index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then

            phi_x_plus(i,j,k) = D1(i+1,j,k)*inv_dx
            phi_x_minus(i,j,k) = D1(i,j,k)*inv_dx
   
         
        endif
      enddo
c     } end loop over narrow band points

c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm3dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    klo_D1_gb, khi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    klo_phi_gb, khi_phi_gb, 
     &                    order, y_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D1)

c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)
        k = index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then

            phi_y_minus(i,j,k) = D1(i,j,k)*inv_dy
            phi_y_plus(i,j,k) = D1(i,j+1,k)*inv_dy
   
        endif
      enddo
c     } end loop over narrow band points 

c----------------------------------------------------
c    compute phi_z_plus and phi_z_minus
c----------------------------------------------------
c     compute first undivided differences in z-direction
      call lsm3dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb, 
     &                    klo_D1_gb, khi_D1_gb, 
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb, 
     &                    klo_phi_gb, khi_phi_gb,
     &                    order, z_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D1)

c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)
        k = index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then

            phi_z_plus(i,j,k) = D1(i,j,k+1)*inv_dz
            phi_z_minus(i,j,k) = D1(i,j,k)*inv_dz
   
        endif
      enddo
c     } end loop over narrow band points

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm3dHJENO2LOCAL() computes the forward (plus) and backward (minus)
c  second-order Hamilton-Jacobi ENO approximations to the gradient of 
c  phi.
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    phi_*_plus (out):   components of grad(phi) in plus direction 
c    phi_*_minus (out):  components of grad(phi) in minus direction
c    phi (in):           phi
c    D1 (in):            scratch space for holding undivided first-differences
c    D2 (in):            scratch space for holding undivided second-differences
c    dx, dy, dz (in):    grid spacing
c    *_gb (in):          index range for ghostbox
c    index_*(in):        coordinates of local (narrow band) points
c    n*_index[012](in):  index range of points in index_* that are in
c                        level [012] of the narrow band
c    narrow_band(in):    array that marks voxels outside desired fillbox
c    mark_*(in):         upper limit narrow band value for voxels in 
c                        the appropriate fillbox
c
c  NOTES:
c   - it is assumed that BOTH the plus and minus derivatives have
c     the same fillbox
c   - index_* arrays range at minimum from nlo_index0 to nhi_index2
c
c***********************************************************************
      subroutine lsm3dHJENO2LOCAL(
     &  phi_x_plus, phi_y_plus, phi_z_plus,
     &  ilo_grad_phi_plus_gb, ihi_grad_phi_plus_gb, 
     &  jlo_grad_phi_plus_gb, jhi_grad_phi_plus_gb,
     &  klo_grad_phi_plus_gb, khi_grad_phi_plus_gb,
     &  phi_x_minus, phi_y_minus, phi_z_minus,
     &  ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb, 
     &  jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb, 
     &  klo_grad_phi_minus_gb, khi_grad_phi_minus_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  D1,
     &  ilo_D1_gb, ihi_D1_gb, 
     &  jlo_D1_gb, jhi_D1_gb,
     &  klo_D1_gb, khi_D1_gb,
     &  D2,
     &  ilo_D2_gb, ihi_D2_gb, 
     &  jlo_D2_gb, jhi_D2_gb,
     &  klo_D2_gb, khi_D2_gb,
     &  dx, dy, dz,
     &  index_x,
     &  index_y, 
     &  index_z, 
     &  nlo_index0, nhi_index0,
     &  nlo_index1, nhi_index1,
     &  nlo_index2, nhi_index2,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
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
      integer klo_grad_phi_plus_gb, khi_grad_phi_plus_gb
      integer ilo_grad_phi_minus_gb, ihi_grad_phi_minus_gb
      integer jlo_grad_phi_minus_gb, jhi_grad_phi_minus_gb
      integer klo_grad_phi_minus_gb, khi_grad_phi_minus_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      integer ilo_D1_gb, ihi_D1_gb
      integer jlo_D1_gb, jhi_D1_gb
      integer klo_D1_gb, khi_D1_gb
      integer ilo_D2_gb, ihi_D2_gb
      integer jlo_D2_gb, jhi_D2_gb
      integer klo_D2_gb, khi_D2_gb
      integer nlo_index0, nhi_index0
      integer nlo_index1, nhi_index1
      integer nlo_index2, nhi_index2
      integer index_x(nlo_index0:nhi_index2)
      integer index_y(nlo_index0:nhi_index2)
      integer index_z(nlo_index0:nhi_index2)
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
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real D1(ilo_D1_gb:ihi_D1_gb,
     &        jlo_D1_gb:jhi_D1_gb,
     &        klo_D1_gb:khi_D1_gb)
      real D2(ilo_D2_gb:ihi_D2_gb,
     &        jlo_D2_gb:jhi_D2_gb,
     &        klo_D2_gb:khi_D2_gb)      
      real dx, dy, dz
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_D2
      integer*1 mark_D1
      integer*1 mark_fb
      
      real inv_dx, inv_dy, inv_dz
      integer i,j,k,l
      real half
      parameter (half=0.5d0)
      integer order_1, order_2
      parameter (order_1=1,order_2=2)
      integer x_dir, y_dir, z_dir
      parameter (x_dir=1,y_dir=2,z_dir=3)

c     compute inv_dx, inv_dy, and inv_dz
      inv_dx = 1.0d0/dx
      inv_dy = 1.0d0/dy
      inv_dz = 1.0d0/dz

c----------------------------------------------------
c    compute phi_x_plus and phi_x_minus
c----------------------------------------------------
c     compute first undivided differences in x-direction
c     for now, these are computed everywhere (and not only in narrow band)
      call lsm3dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    klo_D1_gb, khi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    klo_phi_gb, khi_phi_gb,
     &                    order_1, x_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D1) 

c     compute second undivided differences x-direction
      call lsm3dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    klo_D2_gb, khi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    klo_D1_gb, khi_D1_gb,
     &                    order_2, x_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D2) 
c    loop over narrow band level 0 points {
      do l=nlo_index0, nhi_index0   
        i=index_x(l)
        j=index_y(l)
        k=index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then

c             phi_x_plus
              if (abs(D2(i,j,k)).lt.abs(D2(i+1,j,k))) then
                phi_x_plus(i,j,k) = (D1(i+1,j,k) 
     &                            - half*D2(i,j,k))*inv_dx
              else
                phi_x_plus(i,j,k) = (D1(i+1,j,k) 
     &                            - half*D2(i+1,j,k))*inv_dx
              endif

c             phi_x_minus
              if (abs(D2(i-1,j,k)).lt.abs(D2(i,j,k))) then
                phi_x_minus(i,j,k) = (D1(i,j,k) 
     &                             + half*D2(i-1,j,k))*inv_dx
              else
                phi_x_minus(i,j,k) = (D1(i,j,k) 
     &                             + half*D2(i,j,k))*inv_dx
              endif
        endif      
      enddo
c     } end loop over indexed points


c----------------------------------------------------
c    compute phi_y_plus and phi_y_minus
c----------------------------------------------------
c     compute first undivided differences in y-direction
      call lsm3dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    klo_D1_gb, khi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    klo_phi_gb, khi_phi_gb,
     &                    order_1, y_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D1) 
     
c     compute second undivided differences in y-direction
      call lsm3dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    klo_D2_gb, khi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    klo_D1_gb, khi_D1_gb,
     &                    order_2, y_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D2)
      
c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0       
        i = index_x(l)
        j = index_y(l)
        k = index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then
c             phi_y_plus
              if (abs(D2(i,j,k)).lt.abs(D2(i,j+1,k))) then
                phi_y_plus(i,j,k) = (D1(i,j+1,k) 
     &                            - half*D2(i,j,k))*inv_dy
              else
                phi_y_plus(i,j,k) = (D1(i,j+1,k) 
     &                            - half*D2(i,j+1,k))*inv_dy
              endif

c             phi_y_minus
              if (abs(D2(i,j-1,k)).lt.abs(D2(i,j,k))) then
                phi_y_minus(i,j,k) = (D1(i,j,k) 
     &                             + half*D2(i,j-1,k))*inv_dy
              else
                phi_y_minus(i,j,k) = (D1(i,j,k) 
     &                             + half*D2(i,j,k))*inv_dy
              endif
        endif      
      enddo
c     } end loop over narrow band points


c----------------------------------------------------
c    compute phi_z_plus and phi_z_minus
c----------------------------------------------------
c     compute first undivided differences in z-direction
      call lsm3dComputeDnLOCAL(D1, 
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    klo_D1_gb, khi_D1_gb,
     &                    phi,
     &                    ilo_phi_gb, ihi_phi_gb, 
     &                    jlo_phi_gb, jhi_phi_gb,
     &                    klo_phi_gb, khi_phi_gb,
     &                    order_1, z_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index2,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D1)
     
c     compute second undivided differences in z-direction
      call lsm3dComputeDnLOCAL(D2, 
     &                    ilo_D2_gb, ihi_D2_gb, 
     &                    jlo_D2_gb, jhi_D2_gb,
     &                    klo_D2_gb, khi_D2_gb,
     &                    D1,
     &                    ilo_D1_gb, ihi_D1_gb, 
     &                    jlo_D1_gb, jhi_D1_gb,
     &                    klo_D1_gb, khi_D1_gb,
     &                    order_2, z_dir,
     &                    index_x, index_y, index_z,
     &                    nlo_index0, nhi_index1,
     &                    narrow_band,     
     &                    ilo_nb_gb, ihi_nb_gb, 
     &                    jlo_nb_gb, jhi_nb_gb,
     &                    klo_nb_gb, khi_nb_gb,
     &                    mark_D2) 
     
c    loop over  narrow band level 0 points only {
      do l = nlo_index0, nhi_index0     
        i = index_x(l)
        j = index_y(l)
        k = index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then

c             phi_z_plus
              if (abs(D2(i,j,k)).lt.abs(D2(i,j,k+1))) then
                phi_z_plus(i,j,k) = (D1(i,j,k+1) 
     &                            - half*D2(i,j,k))*inv_dz
              else
                phi_z_plus(i,j,k) = (D1(i,j,k+1) 
     &                            - half*D2(i,j,k+1))*inv_dz
              endif

c             phi_z_minus
              if (abs(D2(i,j,k-1)).lt.abs(D2(i,j,k))) then
                phi_z_minus(i,j,k) = (D1(i,j,k) 
     &                             + half*D2(i,j,k-1))*inv_dz
              else
                phi_z_minus(i,j,k) = (D1(i,j,k) 
     &                             + half*D2(i,j,k))*inv_dz
              endif
        endif      
      enddo
c     } end loop over grid 

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm3dCentralGradOrder2LOCAL() computes the second-order central 
c  approximation to the gradient of phi.
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    phi_* (out):      components of grad(phi) 
c    phi (in):         phi
c    dx, dy, dz (in):  grid spacing
c    *_gb (in):        index range for ghostbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):     index range of points to loop over in index_*
c    narrow_band(in):  array that marks voxels outside desired fillbox
c    mark_fb(in):      upper limit narrow band value for voxels in 
c                      fillbox
c
c***********************************************************************
      subroutine lsm3dCentralGradOrder2LOCAL(
     &  phi_x, phi_y, phi_z,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  klo_grad_phi_gb, khi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  dx, dy, dz,
     &  index_x,
     &  index_y, 
     &  index_z, 
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
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
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb

c     local variables      
      integer i,j,k,l
      real dx_factor, dy_factor, dz_factor

c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy
      dz_factor = 0.5d0/dz

c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)
        k=index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then
          phi_x(i,j,k) = (phi(i+1,j,k) - phi(i-1,j,k))*dx_factor
          phi_y(i,j,k) = (phi(i,j+1,k) - phi(i,j-1,k))*dy_factor
          phi_z(i,j,k) = (phi(i,j,k+1) - phi(i,j,k-1))*dz_factor
        endif  
      enddo
c     } end loop over indexed points
      
      
      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dCentralGradOrder4LOCAL() computes the second-order, central, 
c  finite difference approximation to the gradient of phi.
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    phi_* (out):      components of grad(phi) 
c    phi (in):         phi
c    dx, dy, dz (in):  grid spacing
c    *_gb (in):        index range for ghostbox
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):     index range of points to loop over in index_*
c    narrow_band(in):  array that marks voxels outside desired fillbox
c    mark_fb(in):      upper limit narrow band value for voxels in 
c                      fillbox
c
c***********************************************************************
      subroutine lsm3dCentralGradOrder4LOCAL(
     &  phi_x, phi_y, phi_z,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb, 
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  klo_grad_phi_gb, khi_grad_phi_gb,
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  dx, dy, dz,
     &  index_x,
     &  index_y, 
     &  index_z, 
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
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
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb
      
      integer i,j,k,l
      real dx_factor, dy_factor, dz_factor
      real eight
      parameter (eight = 8.0d0)

c     compute denominator values
      dx_factor = 0.0833333333333333333333d0/dx
      dy_factor = 0.0833333333333333333333d0/dy
      dz_factor = 0.0833333333333333333333d0/dz

c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)
        k=index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then

            phi_x(i,j,k) = ( -phi(i+2,j,k) + eight*phi(i+1,j,k) 
     &                       +phi(i-2,j,k) - eight*phi(i-1,j,k) )
     &                   * dx_factor
            phi_y(i,j,k) = ( -phi(i,j+2,k) + eight*phi(i,j+1,k) 
     &                       +phi(i,j-2,k) - eight*phi(i,j-1,k) )
     &                   * dy_factor
            phi_z(i,j,k) = ( -phi(i,j,k+2) + eight*phi(i,j,k+1) 
     &                       +phi(i,j,k-2) - eight*phi(i,j,k-1) )
     &                   * dz_factor
   
          
        endif
      enddo
c     } end loop over indexed points 

      return
      end
c } end subroutine
c***********************************************************************

c***********************************************************************
c
c  lsm3dLaplacianOrder2LOCAL() computes the second-order, central, 
c  finite difference approximation to the Laplacian of phi.
c  The routine loops only over local (narrow band) points.
c
c  Arguments:
c    laplacian_phi (out):  Laplacian of phi
c    phi (in):             phi
c    dx (in):              grid spacing
c    *_gb (in):            index range for ghostbox
c    index_[xyz](in):     [xyz] coordinates of local (narrow band) points
c    n*_index(in):        index range of points to loop over in index_*
c    narrow_band(in):     array that marks voxels outside desired fillbox
c    mark_fb(in):         upper limit narrow band value for voxels in 
c                         fillbox
c
c***********************************************************************
      subroutine lsm3dLaplacianOrder2LOCAL(
     &  laplacian_phi,
     &  ilo_laplacian_phi_gb, ihi_laplacian_phi_gb, 
     &  jlo_laplacian_phi_gb, jhi_laplacian_phi_gb, 
     &  klo_laplacian_phi_gb, khi_laplacian_phi_gb, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb, 
     &  klo_phi_gb, khi_phi_gb,
     &  dx, dy, dz,
     &  index_x,
     &  index_y, 
     &  index_z, 
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _laplacian_phi_gb refers to ghostbox for laplacian_phi data
c     _phi_gb refers to ghostbox for phi data
c     _fb refers to fill-box for grad_phi data
      integer ilo_laplacian_phi_gb, ihi_laplacian_phi_gb
      integer jlo_laplacian_phi_gb, jhi_laplacian_phi_gb
      integer klo_laplacian_phi_gb, khi_laplacian_phi_gb
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      real laplacian_phi(ilo_laplacian_phi_gb:ihi_laplacian_phi_gb,
     &                   jlo_laplacian_phi_gb:jhi_laplacian_phi_gb,
     &                   klo_laplacian_phi_gb:khi_laplacian_phi_gb)
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real dx, dy, dz
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb
            
      integer i, j, k, l
      real inv_dx_sq
      real inv_dy_sq
      real inv_dz_sq

c     compute denominator values
      inv_dx_sq = 1.0d0/dx/dx
      inv_dy_sq = 1.0d0/dy/dy
      inv_dz_sq = 1.0d0/dz/dz

c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)
        k=index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then

            laplacian_phi(i,j,k) = 
     &        inv_dx_sq *
     &        ( phi(i+1,j,k) - 2.0d0*phi(i,j,k) + phi(i-1,j,k) ) 
     &      + inv_dy_sq *
     &        ( phi(i,j+1,k) - 2.0d0*phi(i,j,k) + phi(i,j-1,k) ) 
     &      + inv_dz_sq *
     &        ( phi(i,j,k+1) - 2.0d0*phi(i,j,k) + phi(i,j,k-1) ) 
           
        endif
      enddo
c     } end loop over indexed points 

      return
      end
c } end subroutine
c***********************************************************************


c***********************************************************************
c
c  lsm3dComputeAveGradPhiLocal() computes the average of the 
c  second-order, central, finite difference approximation to the 
c  gradient of phi within the narrow band.
c  The routine loops only over local (narrow band) points. 
c
c  Arguments:
c    phi (in):          phi
c    grad_phi_ave(out): average of the gradient 
c    dx, dy, dz (in):   grid spacing
c    index_[xyz](in):  [xyz] coordinates of local (narrow band) points
c    n*_index(in):     index range of points to loop over in index_*
c    narrow_band(in):  array that marks voxels outside desired fillbox
c    mark_fb(in):      upper limit narrow band value for voxels in 
c                      fillbox
c    *_gb (in):        index range for ghostbox
c
c***********************************************************************
      subroutine lsm3dComputeAveGradPhiLocal(
     &  grad_phi_ave, 
     &  phi,
     &  ilo_phi_gb, ihi_phi_gb, 
     &  jlo_phi_gb, jhi_phi_gb,
     &  klo_phi_gb, khi_phi_gb,
     &  dx, dy, dz,
     &  index_x,
     &  index_y, 
     &  index_z, 
     &  nlo_index, nhi_index,
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb, 
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

      real grad_phi_ave
      integer ilo_phi_gb, ihi_phi_gb
      integer jlo_phi_gb, jhi_phi_gb
      integer klo_phi_gb, khi_phi_gb
      real phi(ilo_phi_gb:ihi_phi_gb,
     &         jlo_phi_gb:jhi_phi_gb,
     &         klo_phi_gb:khi_phi_gb)
      real dx, dy, dz
      integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb

c     local variables      
      integer i,j,k,l,count
      real dx_factor, dy_factor, dz_factor
      real phi_x, phi_y, phi_z
      

c     compute denominator values
      dx_factor = 0.5d0/dx
      dy_factor = 0.5d0/dy
      dz_factor = 0.5d0/dz

      grad_phi_ave = 0.d0
      count = 0
c     { begin loop over indexed points
      do l= nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)
        k=index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then
      
          phi_x = (phi(i+1,j,k) - phi(i-1,j,k))*dx_factor
          phi_y = (phi(i,j+1,k) - phi(i,j-1,k))*dy_factor
          phi_z = (phi(i,j,k+1) - phi(i,j,k-1))*dz_factor

          grad_phi_ave = grad_phi_ave + sqrt(phi_x*phi_x + phi_y*phi_y +
     &                                       phi_z*phi_z)
          count = count+1
     
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
c
c  lsm3dGradientMagnitudeLocal() computes magnitude of the gradient of phi.
c
c  Arguments:
c    phi_* (in):         components of grad(phi) 
c    grad_phi_mag (out): gradient magnitude
c    *_gb (in):          index range for ghostbox
c    index_*(in):        coordinates of local (narrow band) points
c    n*_index(in):       index range of points to loop over in index_*
c    narrow_band(in):    array that marks voxels outside desired fillbox
c    mark_fb(in):        upper limit narrow band value for voxels in 
c                        fillbox
c
c***********************************************************************
      subroutine lsm3dGradientMagnitudeLocal(
     &  phi_x, phi_y, phi_z, grad_phi_mag,
     &  ilo_grad_phi_gb, ihi_grad_phi_gb,
     &  jlo_grad_phi_gb, jhi_grad_phi_gb,
     &  klo_grad_phi_gb, khi_grad_phi_gb,
     &  index_x,
     &  index_y,
     &  index_z,
     &  nlo_index, nhi_index, 
     &  narrow_band,
     &  ilo_nb_gb, ihi_nb_gb, 
     &  jlo_nb_gb, jhi_nb_gb,     
     &  klo_nb_gb, khi_nb_gb,
     &  mark_fb)
c***********************************************************************
c { begin subroutine
      implicit none

c     _grad_phi_gb refers to ghostbox for grad_phi data
      integer ilo_grad_phi_gb, ihi_grad_phi_gb
      integer jlo_grad_phi_gb, jhi_grad_phi_gb
      integer klo_grad_phi_gb, khi_grad_phi_gb
       real phi_x(
     &                    ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                    jlo_grad_phi_gb:jhi_grad_phi_gb,
     &                    klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_y(
     &                    ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                    jlo_grad_phi_gb:jhi_grad_phi_gb,
     &                    klo_grad_phi_gb:khi_grad_phi_gb)
      real phi_z(
     &                    ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                    jlo_grad_phi_gb:jhi_grad_phi_gb,
     &                    klo_grad_phi_gb:khi_grad_phi_gb)
      real grad_phi_mag(ilo_grad_phi_gb:ihi_grad_phi_gb,
     &                              jlo_grad_phi_gb:jhi_grad_phi_gb,
     &                              klo_grad_phi_gb:khi_grad_phi_gb)
      
       integer nlo_index, nhi_index
      integer index_x(nlo_index:nhi_index)
      integer index_y(nlo_index:nhi_index)
      integer index_z(nlo_index:nhi_index)
      integer ilo_nb_gb, ihi_nb_gb
      integer jlo_nb_gb, jhi_nb_gb
      integer klo_nb_gb, khi_nb_gb
      integer*1 narrow_band(ilo_nb_gb:ihi_nb_gb,
     &                      jlo_nb_gb:jhi_nb_gb,
     &                      klo_nb_gb:khi_nb_gb)
      integer*1 mark_fb

c     local variables      
      real phi_x_sq, phi_y_sq, phi_z_sq
      real zero_tol, zero, tmp
      parameter (zero_tol=1.d-11, zero = 0.d0)
      integer i,j,k,l

c     { begin loop over indexed points
       do l=nlo_index, nhi_index      
        i=index_x(l)
        j=index_y(l)
        k=index_z(l)

c       include only fill box points (marked appropriately)
        if( narrow_band(i,j,k) .le. mark_fb ) then

          phi_x_sq = phi_x(i,j,k)*phi_x(i,j,k)
          phi_y_sq = phi_y(i,j,k)*phi_y(i,j,k)
          phi_z_sq = phi_z(i,j,k)*phi_z(i,j,k)

          if( phi_x_sq .lt. zero_tol) then
            phi_x_sq = zero
          endif    

          if( phi_y_sq .lt. zero_tol) then
            phi_y_sq = zero
          endif 

          if( phi_z_sq .lt. zero_tol) then
            phi_z_sq = zero
          endif   
 
          tmp = sqrt( phi_x_sq + phi_y_sq + phi_z_sq)
          if( tmp .lt. zero_tol ) then
            grad_phi_mag(i,j,k) = zero
          else
            grad_phi_mag(i,j,k) = tmp
          endif

        endif
      enddo
c     } end loop over indexed points 

      return
      end
c } end subroutine
c***********************************************************************
