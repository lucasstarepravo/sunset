module common_vars
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains variables common to all modules within the sunset code

  use iso_c_binding
  use common_parameter
  use kind_parameters      

  implicit none

  !! Evolved fluid quantities
  real(rkind), dimension(:), allocatable, target :: u,v,w,lnro,roE
  real(rkind), dimension(:,:), allocatable :: Yspec  
  real(rkind), dimension(:),allocatable :: alpha_out
  
  !! Number of species
  integer(ikind) :: nspec
     
  !! Secondary fluid quantities
  real(rkind), dimension(:), allocatable, target :: p,T
  real(rkind), dimension(:), allocatable :: divvel  
  
  !! Transport properties
  real(rkind), dimension(:), allocatable :: Rgas_mix,cp,visc,lambda_th
  real(rkind), dimension(:,:), allocatable :: Mdiff
  
  !! Transport properties - arrays covering the species
  real(rkind), dimension(:), allocatable :: molar_mass,one_over_Lewis_number
  real(rkind), dimension(:,:),allocatable :: coef_cp,coef_h !! indexing: ispec,j-exponent
  integer(ikind) :: polyorder_cp,ncoefs_cp  !! polynomial order,number of coefs
  real(rkind) :: T_low,T_high
  
  !! Chemical kinetics control data
  integer(ikind) :: nsteps
  real(rkind),dimension(:,:), allocatable :: Arrhenius_coefs
  
  !! Right-hand-sides
  real(rkind),dimension(:),allocatable :: rhs_lnro,rhs_u,rhs_v,rhs_w,rhs_roE
  real(rkind),dimension(:,:),allocatable :: rhs_Yspec
    
  !! Discretisation properties
  real(rkind), dimension(:,:), allocatable, target :: rp,rnorm
  real(rkind), dimension(:), allocatable, target   :: h,filter_coeff,s
  integer(ikind),dimension(:),allocatable :: node_type !! Identify whether node is boundary, fluid etc...
  integer(ikind),dimension(:),allocatable :: zlayer_index_global,ilayer_index !! Identify where in the z-stack the node is
  integer(ikind),dimension(:),allocatable :: boundary_list,internal_list !! Lists for quick looping
  real(rkind) :: dz   !! FD spacing in third dimension
  integer(ikind) :: nz,nz_global
  
  !! Numbers of nodes and neighbour lists
  integer(ikind) :: np,npfb,nb,nplink  !! THESE ARE ALL LOCAL
  integer(ikind) :: np_global,npfb_global,nb_global !! THESE ARE GLOBAL
  integer(ikind) :: npfb_layer  !! THESE ARE ALL LOCAL
  integer(ikind) :: npfb_layer_global !! THESE ARE GLOBAL


  !! Variables related to stencil sizes 
  real(rkind) :: h0,sup_size,h3,hovs,ss,h2,hovs_bound

  !! Parameters related to time and some forces etc
  real(rkind) :: time,dt,dt_previous,time_end,dt_out
  real(rkind) :: time_star !! Dimensionless time (for outputs...)
  real(rkind) :: umax,smax,cmax,smin                  !! maximum velocity,node-spacing,sound speed
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
  real(rkind),dimension(:,:,:),allocatable :: ij_w_grad,ij_wb_grad2
  real(rkind),dimension(:,:),allocatable :: ij_w_hyp,ij_w_lap
  real(rkind),dimension(:,:),allocatable :: ij_w_grad_sum,ij_wb_grad2_sum
  real(rkind),dimension(:),allocatable :: ij_w_hyp_sum,ij_w_lap_sum
  
  !! Finite Difference weightings 
  real(rkind),dimension(:),allocatable :: ij_fd_grad,ij_fd_grad2,ij_fd_hyp         
#define FDORDER 8               
  !! Size of Stencil
#if FDORDER==4
  integer(ikind),parameter :: ij_count_fd = 5
#elif FDORDER==6
  integer(ikind),parameter :: ij_count_fd = 7
#elif FDORDER==8
  integer(ikind),parameter :: ij_count_fd = 9
#elif FDORDER==10
  integer(ikind),parameter :: ij_count_fd = 11
#elif FDORDER==12
  integer(ikind),parameter :: ij_count_fd = 13
#endif    
  
  !! Parents and boundaries... 
  integer(ikind),dimension(:),allocatable :: irelation,vrelation  ! used for periodic and symmetric boundaries
  real(rkind) :: xmin,xmax,ymin,ymax  !! Global domain size (required for NRBCs)
  real(rkind) :: L_domain_x,L_domain_y
  integer(ikind) :: xbcond,ybcond !! BC flags for "simple" geometries without boundary nodes...
  integer(ikind),dimension(:),allocatable :: btype !! What type of BC is node i?
  integer(ikind),dimension(:),allocatable :: fd_parent !! pointer to the boundary node which is parent 
  real(rkind),dimension(:),allocatable :: T_bound
  real(rkind) :: p_outflow   !! Desired pressure on outflow boundary
  
  !! Profiling and openMP parallelisation
  real(rkind) ts_start,ts_end,t_run,t_per_dt,t_last_X
  integer(ikind) :: n_threads  
  real(rkind) :: segment_tstart,segment_tend
  real(rkind),dimension(10) :: segment_time_local,segment_time_global
  real(rkind) :: cputimecheck
  
  !! MPI decomposition related variables
  integer(ikind) :: nprocs,iproc,ierror,iproc_in_sheet  !! processes, this process id, error int,process id in sheet
  integer(ikind) :: nprocsX,nprocsY,nprocsZ,iprocX,iprocY,iprocZ     !! decomposition grid sizes, and indices
  integer(ikind) :: np_nohalo !! nodes with no halos  
  real(rkind) :: XL_thisproc,XR_thisproc,YU_thisproc,YD_thisproc,ZF_thisproc,ZB_thisproc
  real(rkind),dimension(:),allocatable :: XL,XR,YU,YD,ZF,ZB
  integer(ikind),dimension(:),allocatable :: iproc_S_LR,iproc_R_LR,iproc_S_UD,iproc_R_UD !! Neighbouring processors
  integer(ikind),dimension(:),allocatable :: iproc_S_FB,iproc_R_FB 
  integer(ikind),dimension(:,:),allocatable :: halo_lists_LR,halo_lists_UD,halo_lists_FB  !! Lists of halo nodes 
  integer(ikind),dimension(:),allocatable :: nhalo_LR,nhalo_UD,inhalo_LR,inhalo_UD  !! Halo sizes, outgoing, incoming
  integer(ikind),dimension(:),allocatable :: nhalo_FB,inhalo_FB
  integer(ikind),dimension(:),allocatable :: nrecstart  !! Indexing for halos
  
  
          
end module common_vars
