module common_2d

  use iso_c_binding
  use common_parameter
  use kind_parameters      

  implicit none

  !! Evolved fluid properties
  real(rkind), dimension(:), allocatable, target :: u,v,w,lnro,roE,Y0  
  
  !! Secondary fluid properties
  real(rkind), dimension(:), allocatable, target :: p,visc,T
  real(rkind), dimension(:), allocatable :: divvel  
  
  !! Discretisation properties
  real(rkind), dimension(:,:), allocatable, target :: rp,rnorm
  real(rkind), dimension(:), allocatable, target   :: h,filter_coeff,s
  integer(ikind),dimension(:),allocatable :: node_type !! Identify whether node is boundary, fluid etc...
  integer(ikind),dimension(:),allocatable :: zlayer_index
  integer(ikind),dimension(:),allocatable :: boundary_list,internal_list
  real(rkind) :: dz   !! FD spacing in third dimension
  integer(ikind) :: nz
  
  !! Numbers of nodes and neighbour lists
  integer(ikind) :: np,npfb,npfb_esti,nb,nplink  !! THESE ARE ALL LOCAL!!
  integer(ikind) :: np_global,npfb_global,nb_global !! THESE ARE GLOBAL
  integer(ikind) :: np_layer,npfb_layer,nb_layer  !! THESE ARE ALL LOCAL!!
  integer(ikind) :: np_layer_global,npfb_layer_global,nb_layer_global !! THESE ARE GLOBAL


  !! Variables related to stencil sizes 
  real(rkind) :: h0,sup_size,h3,hovs,ss,h2,hovs_bound

  !! Parameters related to time and some forces etc
  real(rkind) :: time,dt,time_end,dt_out,dt_mout
  real(rkind) :: umax,smax                    !! maximum velocity and node spacing   
  integer(ikind) :: itime
  real(rkind) :: emax_nm1,emax_n,emax_np1  !! errors for PID controller
  real(rkind) :: eflow_nm1,eflow_n,sum_eflow !! errors for PID to control <u> (constant-ish flow rate)
  real(rkind), dimension(dims) :: driving_force
  real(rkind) :: mean_int_energy0

  !! Neighbour numbers and lists
  integer(ikind),dimension(:),allocatable :: ij_count
  integer(ikind),dimension(:,:),allocatable :: ij_link
  integer(ikind),dimension(:,:),allocatable :: ij_link_fd

  !! LABFM weightings for derivative operators
  real(rkind),dimension(:,:,:),allocatable :: ij_w_grad,ij_w_grad2
  real(rkind),dimension(:,:),allocatable :: ij_w_hyp,ij_w_lap
  real(rkind),dimension(:,:),allocatable :: ij_w_grad_sum,ij_w_grad2_sum
  real(rkind),dimension(:),allocatable :: ij_w_hyp_sum,ij_w_lap_sum
  
  !! Finite Difference weightings 
  integer(ikind) :: ij_count_fd ! Size of FD stencil  
  real(rkind),dimension(:),allocatable :: ij_fd_grad,ij_fd_grad2,ij_fd_hyp
  
  !! Parents and boundaries... 
  integer(ikind),dimension(:),allocatable :: irelation,vrelation  ! used for periodic and symmetric boundaries
  real(rkind) :: xmin,xmax,ymin,ymax  !! Global domain size (required for NRBCs)
  integer(ikind) :: xbcond,ybcond !! BC flags for "simple" geometries without boundary nodes...
  real(rkind),dimension(:,:),allocatable :: L  !! The "L" in NSCBC formulation
  integer(ikind),dimension(:),allocatable :: btype !! What type of BC is node i?
  integer(ikind),dimension(:),allocatable :: fd_parent !! pointer to the boundary node which is parent 
  real(rkind),dimension(:),allocatable :: T_bound
  
  !! Profiling and openMP parallelisation
  real(rkind) ts_start,ts_end,t_run,t_per_dt,t_last_X
  integer(ikind) :: n_threads  
  real(rkind) :: segment_tstart,segment_tend
  real(rkind),dimension(10) :: segment_time_local,segment_time_global
  real(rkind) :: cputimecheck
  
  !! MPI decomposition related variables
  integer(ikind) :: nprocs,iproc,ierror  !! processes, this process id, error int
  integer(ikind) :: nprocsX,nprocsY,nprocsZ,iprocX,iprocY,iprocZ     !! decomposition grid sizes, and indices
  integer(ikind) :: np_nohalo !! nodes with no halos  
  real(rkind) :: XL_thisproc,XR_thisproc,YU_thisproc,YD_thisproc
  integer(ikind),dimension(:),allocatable :: iproc_S_LR,iproc_R_LR,iproc_S_UD,iproc_R_UD !! Neighbouring processors
  integer(ikind),dimension(:,:),allocatable :: halo_lists_LR,halo_lists_UD  !! Lists of halo nodes 
  integer(ikind),dimension(:),allocatable :: nhalo_LR,nhalo_UD,inhalo_LR,inhalo_UD  !! Halo sizes, outgoing, incoming
  integer(ikind),dimension(:),allocatable :: nrecstart_UDLR  !! Indexing for halos
  
          
end module common_2d
