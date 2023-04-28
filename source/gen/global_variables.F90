      module global_variables 
      use kind_parameters      
      implicit none

      integer(ikind), parameter :: npar=9999999  !! Only used in source/gen/datclass.F90 (up to npar in a slice)
      integer(ikind) :: np, npfb, nb,nbio

      real(rkind) :: dx,dx0,dxb,dxio,dx_in,dx_out,dx_wall
      real(rkind) :: grx, gry
      real(rkind) :: xb_min, xb_max, yb_min, yb_max, xl
      real(rkind), dimension(:), allocatable :: xp, yp,thta,xnorm,ynorm,dxp
      integer(ikind),dimension(:),allocatable :: node_type

      !! jack's boundary condition framework
      real(rkind),dimension(:,:),allocatable, target :: b_node,b_edge
      integer(ikind),dimension(:),allocatable,target :: b_type,b_periodic_parent
      integer(ikind) :: nb_patches
      integer(ikind) :: nb_circles,nb_blobs
      real(rkind),dimension(:,:),allocatable :: c_centre
      real(rkind),dimension(:),allocatable :: c_radius,b_theta,c_theta
      real(rkind),dimension(:,:),allocatable :: blob_centre,blob_coeffs
      real(rkind),dimension(:),allocatable :: blob_rotation
      integer(ikind),dimension(:),allocatable :: blob_ellipse
      integer(ikind) :: n_blob_coefs
      
      !! Blob-perimeter
      real(rkind),dimension(:),allocatable :: Sblob


      end module global_variables 
