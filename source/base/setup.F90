module setup
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to read in node distributions, calls to MPI domain decomposion
  !! routines, read in chemistry, and specify initial conditions.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  use neighbours
  use output
#ifdef mp
  use mpi
  use mpi_transfers
#endif    
  implicit none
  
contains
!! ------------------------------------------------------------------------------------------------
  subroutine initial_setup       
     !! Initialises key simulation parameters, and loads control data
     integer(ikind) :: i
    
     !! Set up integers for processor and number of processors
#ifdef mp
     call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
     call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierror)
#endif

     !! Load data from the control file
     call load_control_data
     

     !! Time begins at zero
     time = zero;itime=0
     dt_out = 0.001d0*Time_char         !! Frequency to output fields
     time_end = 1.0d0*Time_char
  
     !! Particles per smoothing length and supportsize/h
     hovs = 2.7   !! 2.4 for 6th order, 2.7 for 8th?
     hovs_bound = 2.4  !! Boundaries get unstable if h is too big
     ss = 2.0
     nplink = 2*4*ceiling(ss*hovs)**2  !! # square stencil w/ side length 2*dr*hovs*ss

     !! Initial data for profiling
#ifdef mp
     !$omp parallel  
     n_threads=omp_get_num_threads()
     !$omp end parallel
     write(6,*) "nprocs,iproc,n_threads:",nprocs,iproc,n_threads  
#else  
     !$omp parallel
     n_threads = omp_get_num_threads()
     !$omp end parallel
#endif  
     t_run = zero;t_last_X=zero

     !! Open a files for outputting
     open(unit=21,file='./data_out/time.out')
     open(unit=191,file='data_out/statistics/cputime.out')
     open(unit=192,file='data_out/statistics/dt.out')
     open(unit=193,file='data_out/statistics/masscheck.out')
     open(unit=194,file='data_out/statistics/liftdrag.out')
     open(unit=195,file='data_out/statistics/velcheck.out')  
     open(unit=196,file='data_out/statistics/l2.out')
     open(unit=197,file='data_out/statistics/energy_mean.out')  
     open(unit=199,file='data_out/statistics/species.out')
  
     !! Profiling:
     segment_time_local = zero
     cputimecheck = zero
  
  end subroutine initial_setup
!! ------------------------------------------------------------------------------------------------
  subroutine setup_domain
     !! Reads in boundary patches, builds domain, calls decomposition and 
     !! boundary setup routines
     use mirror_boundaries
     integer(ikind) i,j,ii,jj,npfb_tmp,k,dummy_int,dummy_int2,nm
     real(rkind) :: ns,dummy,rad,x,y,xn,yn
     real(rkind),dimension(dims) :: rij
     real(rkind),dimension(:,:),allocatable :: tmp_vec
     integer(ikind) :: nl_ini,nl_end,nl_iniC,nl_endC,nl_ini2,nl_end2

     !! STEP 1: Load IPART (some params, plus list of nodes + boundary normals)
     !! =======================================================================
     open(13,file='IPART')
     read(13,*) nb,npfb,dummy      !! dummy is largest s(i) in domain...
     read(13,*) xmin,xmax,ymin,ymax
     !! Set the domain lengths
     L_domain_x = (xmax - xmin)*L_char
     L_domain_y = (ymax - ymin)*L_char
     
     read(13,*) xbcond,ybcond
     !! Calculate some useful constants
     h0 = hovs*dummy;sup_size = ss*h0;h2=h0*h0;h3=h2*h0
         
#ifdef mp

     read(13,*) nprocsX,nprocsY,nprocsZ
#ifndef dim3
     if(nprocsZ.ne.1) then
        write(6,*) "ERROR: 2D simulation, but 3D domain decomposition. Stopping."
        call MPI_Abort(MPI_COMM_WORLD, ii, ierror)              
     end if
#endif

     !! Domain decomposition: How many processors in X and Y, and build lists of neighbour processors 
     call processor_mapping

     !! Total numbers of particles...
     nb_global = nb
     npfb_global = npfb !! without 4*nb in FD stencils
     npfb_global = npfb_global + 4*nb_global         

     !! Read the start and end indices of the decomposition schedule from ipart  
     read(13,*) dummy_int
     if(dummy_int.ne.nprocs) then
        write(6,*) "ERROR: nprocs doesn't match ishift decomposition schedule. Stopping."
        stop
     else
        do i=1,nprocs
           if(iproc.eq.i-1) then
              read(13,*) nl_iniC,nl_endC
           else
              read(13,*) dummy_int,dummy_int
           end if                        
        end do
        npfb = nl_endC - nl_iniC + 1
     end if

     nl_ini=nl_iniC;nl_end=nl_endC

     
!     write(6,*) "process",iproc,"start,end",nl_ini,nl_end             
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)
#else
     !! No MPI, so ignore the decomposition schedule in IPART
     read(13,*) dummy_int,dummy_int
     read(13,*) dummy_int
     do i=1,dummy_int
        read(13,*) dummy_int2,dummy_int2
     end do
     nl_ini = 1;nl_end = npfb
#endif    

     !! Allocate local node arrays
     nm = 10
     allocate(rp(nm*npfb,dims),rnorm(nm*npfb,dims),h(nm*npfb),s(nm*npfb));rp=zero;rnorm=zero
     allocate(node_type(nm*npfb));node_type=0
     allocate(fd_parent(nm*npfb));fd_parent=0
          
     !! Load all nodes. Build FD stencils near boundaries on the fly.
     npfb_tmp = npfb
     nb = 0;ii = 0
#ifdef mp
     !! Skip particles left or below this block
     if(nl_ini.ne.1) then 
        do i=1,nl_ini-1
           read(13,*) dummy,dummy,jj,dummy,dummy,dummy
        end do
     end if
#endif     
     do i=nl_ini,nl_end
        ii = ii + 1
        read(13,*) rp(ii,1),rp(ii,2),jj,rnorm(ii,1),rnorm(ii,2),dummy
        h(ii) = dummy*hovs
        s(ii) = dummy
        node_type(ii) = jj
        if(jj.ge.0.and.jj.le.2) then !! If it is a boundary node
           h(ii) = s(ii)*hovs_bound        
           k = ii !! k is the index of the parent node
           nb = nb + 1
           do j=1,4  !! Make 4 additional nodes
              ii = ii + 1
              rp(ii,:) = rp(k,:) + rnorm(k,:)*dble(j)*s(k)   !! Moving along an FD stencil
              rnorm(ii,:)=rnorm(k,:)          !! Copy normals
              h(ii)=h(k);s(ii)=s(k)          !! length-scales
              node_type(ii) = -j           !! and node type
              fd_parent(ii) = k            !! and lineage
              npfb = npfb + 1           
           end do
        end if
     end do


#ifdef mp   
     !! Set the spatial limits of this processor block
     XL_thisproc = minval(rp(1:npfb,1))
     XR_thisproc = maxval(rp(1:npfb,1))    
     YU_thisproc = maxval(rp(floor(0.75*npfb):npfb,2))  !! Adjustments for cyclical columns. 
     YD_thisproc = minval(rp(1:floor(0.25*npfb),2))       !! Should be oK given node ordering
     
!write(6,*) "iproc",iproc,"YU,YD",YU_thisproc,YD_thisproc
     
    
     write(6,*) "process",iproc,"npfb,nb",npfb,nb
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)     
     write(6,*) iproc,nb,npfb     
#else
     write(6,*) nb,npfb     
#endif     
     close(13)     

     !! Construct the 3rd dimension if required (multiple slices of the 2d domain)
#ifdef dim3
     call build_3rd_dimension
#else
     !! Layer sizes (local and global)
     npfb_layer = npfb;npfb_layer_global = npfb_global  
     dz = h0/hovs  
     nz=1    
     nz_global = 1
#endif     

     !! STEP 2: build mirrors, halos and neighbours
     !! =======================================================================
     call create_mirror_particles

#ifdef mp
     ZF_thisproc = maxval(rp(1:npfb,3))
     ZB_thisproc = minval(rp(1:npfb,3))
  
     !! Initial halo build - much too big (based on kernel size x 1.2?)
     call build_halos

#ifdef dim3
     !! Transfer information about z-layer
     call halo_exchange_int(zlayer_index_global)  
     call halo_exchange_int(ilayer_index)   
#endif     
           
     write(6,*) "Proc",iproc,"with",nb,npfb,np_nohalo,np
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)    

     !! Refine halos: compile lists of halo nodes which are used, and send them back to originator
     !! processors, then re-build halos
     call find_neighbours
     call refine_halos

     !! Shrink arrays to fit number of nodes
     call reduce_arrays

     deallocate(ij_link,ij_count)
     
     !! Transfer discretisation properties   
     call halo_exchange(h)
     call halo_exchange(s)
     call halo_exchange_int(node_type)    
     call halo_exchange(rnorm(:,1))
     call halo_exchange(rnorm(:,2)) 
#ifdef dim3
     !! Transfer information about z-layer
     call halo_exchange_int(zlayer_index_global)     
     call halo_exchange_int(ilayer_index)
     ilayer_index(npfb+1:nrecstart(2+2*nprocsY+1)-1)=0 !! Zero out local mirrors and UDLR halo nodes
#endif
     
#else
     np_nohalo = np  
#endif     
     
     !! Build link lists for boundary and internal nodes
     if(nb.ne.0) then
        allocate(boundary_list(nb));boundary_list=0
     end if
     allocate(internal_list(npfb-nb));internal_list=0    
     ii=0;jj=0
     do i=1,npfb
        if(node_type(i).lt.0.or.node_type(i).eq.999) then
           ii=ii+1
           internal_list(ii) = i
        else
           jj=jj+1
           boundary_list(jj) = i
        end if
     end do    
             

     !! Set the global number of boundary nodes
#ifdef mp
     call MPI_ALLREDUCE(nb,nb_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(np,np_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)                    
     call MPI_ALLREDUCE(npfb,npfb_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)                         
#else
     nb_global = nb
#endif
     
     !! We may wish to do a test output here whilst debugging MPI stuff... 
!!     call find_neighbours         
!     call output_layer(1)
!     call output_layer(2)
!     call output_layer(3)     
!     call MPI_BARRIER( MPI_COMM_WORLD, ierror)     
!     call MPI_Abort(MPI_COMM_WORLD, ii, ierror)
           
     !! Initialise PID controller variables
     emax_np1=5.0d-4;emax_n=5.0d-4;emax_nm1=5.0d-4

     !! Initialise PID controller variables for <|u|>
     eflow_nm1 = one
     sum_eflow = zero   
       
           
write(6,*) "sizes",iproc,npfb,np_nohalo,np   
         
                 
     return
  end subroutine setup_domain
!! ------------------------------------------------------------------------------------------------
  subroutine initial_solution
     use mirror_boundaries
     use derivatives
     use thermodynamics
     !! Temporary subroutine whilst developing. Initialises all fields
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro
     
     
     !! Allocate arrays for properties - primary
     allocate(u(np),v(np),w(np),lnro(np),roE(np),divvel(np))
     allocate(Yspec(np,nspec))
     u=zero;v=zero;w=zero;lnro=zero;roE=one;Yspec=one;divvel=zero

     !! Secondary properties
     allocate(alpha_out(np));alpha_out = zero
     allocate(T(np));T=T_ref
     allocate(p(np));p=zero
     
     !! Transport properties
     allocate(visc(npfb));visc = visc_ref
#ifndef isoT
     allocate(lambda_th(npfb))
#endif
#ifdef ms     
     allocate(Mdiff(npfb,nspec))
#endif     
     allocate(cp(np),Rgas_mix(np))
     
     !! Allocate the boundary temperatures
     if(nb.ne.0) then
        allocate(T_bound(nb));T_bound = T_ref                 
     end if
     
     !! =======================================================================
     !! Choose initial conditions
#ifndef restart     
#ifdef react
!     call make_1d_1step_flame
     call load_flame_file    
#else
     call hardcode_initial_conditions     
#endif
#else    
     !! RESTART OPTION. Ask for input number (just hard-coded for now...)
     call load_restart_file(2)
#endif
     !! =======================================================================
            
     !! Set energy from lnro,u,Y,T
     call initialise_energy         
   
     !! Mirrors and halos                        
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchanges_all
#endif          
     
     !! Set the initial velocity divergence
     call calc_divergence(u,v,w,divvel(1:npfb))
     
     !! Mirrors and halos for divvel                   
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchange_divvel
#endif         
     
     !! Set the initial forcing to zero. It will only be changed if PID controller in velcheck is used.
     driving_force(:) = zero

     !! Profiling - re-zero time accumulators
     segment_time_local = zero
     cputimecheck = zero         
     
     !! Initialise the time-stepping (necessary for PID controlled stepping)
     dt = 1.0d-8
     
         
         
     return
  end subroutine initial_solution   
!! ------------------------------------------------------------------------------------------------
  subroutine build_3rd_dimension
     integer(ikind) :: nm,i,k,iz,ilayerm1
     real(rkind) :: dz_local,Lz_dimensionless
     real(rkind),dimension(:,:),allocatable :: rptmp,rnormtmp
     real(rkind),dimension(:),allocatable :: htmp,stmp
     integer(ikind),dimension(:),allocatable :: node_typetmp,fd_parenttmp
     
     !! z spacing
     dz_local = half*(minval(s(1:npfb))+maxval(s(1:npfb)))  !! Match mean resolution
     dz_local = dsqrt(half*(minval(s(1:npfb))**2.0 + maxval(s(1:npfb))**2.0)) !! Match mean resolution (area weighted)
#ifdef mp
     call MPI_ALLREDUCE(dz_local,dz,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)    
     dz = dz/dble(nprocs)                
#else
     dz = dz_local
#endif     

     !! Set the dimensionless extent of the z-domain
     Lz_dimensionless = Lz/L_char

     !! Extent of domain in Z dimension = Lz
     nz = ceiling(Lz_dimensionless/dz/dble(nprocsZ))
     dz = Lz_dimensionless/dble(nz)/dble(nprocsZ)
     nz_global = nz*nprocsZ
                        
     !! Minimum number in z is ij_count_fd/2 + 2
     if(nz.lt.ij_count_fd/2 + 2) then
        nz = ij_count_fd/2 + 2
        dz = Lz_dimensionless/dble(ij_count_fd/2 + 2)
     end if
     
     write(6,*) "iproc",iproc,"Z-domain number and spacing",nz_global,dz
     
     !! Temporary arrays
     allocate(rptmp(npfb,dims),rnormtmp(npfb,dims),htmp(npfb),stmp(npfb))
     allocate(node_typetmp(npfb),fd_parenttmp(npfb))
     !$omp parallel do
     do i=1,npfb
        rptmp(i,:) = rp(i,:)
        rnormtmp(i,:) = rnorm(i,:)
        htmp(i) = h(i)
        stmp(i) = s(i)
        node_typetmp(i) = node_type(i)
        fd_parenttmp(i) = fd_parent(i)
     end do
     !$omp end parallel do

     !! Layer sizes (local and global)
     npfb_layer = npfb;npfb_layer_global = npfb_global

     !! New sizes (local and global)
     npfb = npfb*nz;nb = nb*nz;npfb_global = npfb_global*nz_global
     
     !! Deallocate and reallocate arrays
     nm = 10
     deallocate(rp,rnorm,h,s,node_type,fd_parent)
     allocate(rp(nm*npfb,dims),rnorm(nm*npfb,dims),h(nm*npfb),s(nm*npfb));rp=zero
     allocate(node_type(nm*npfb));node_type=0
     allocate(fd_parent(nm*npfb));fd_parent=0
     allocate(zlayer_index_global(nm*npfb))
     allocate(ilayer_index(nm*npfb));ilayer_index=0
     
     !! Build layers
     k=0
     do iz=1,nz
        ilayerm1 = k
        do i=1,npfb_layer
           k=k+1
           rp(k,1:2) = rptmp(i,1:2)
           rnorm(k,:) = rnormtmp(i,:)
           h(k) = htmp(i)
           s(k) = stmp(i)
           node_type(k) = node_typetmp(i)
           fd_parent(k) = ilayerm1 + fd_parenttmp(i)
           zlayer_index_global(k) = iz + iprocZ*nz  !! z-layer within global stack
           ilayer_index(k) = i                      !! index within layer
           rp(k,3) = dble(iz + iprocZ*nz - 1)*dz
        end do
     end do
        
        
     !! Deallocate temporary arrays
     deallocate(rptmp,rnormtmp,htmp,stmp,fd_parenttmp,node_typetmp)    
     
  
     return
  end subroutine build_3rd_dimension  
!! ------------------------------------------------------------------------------------------------
  subroutine reduce_arrays
     !! Reduces the sizes of arrays for position,s,h,node_type,fd_parent,zlayer_index_global,ilayer_index
     integer(ikind) :: newsize
     real(rkind),dimension(:),allocatable :: tmp_array_real
     real(rkind),dimension(:,:),allocatable :: tmp_array2_real
     integer(ikind),dimension(:),allocatable :: tmp_array_int
     
     !! Set the new-size and allocate temporary arrays
     newsize = np     
     allocate(tmp_array_real(newsize),tmp_array_int(newsize))
     allocate(tmp_array2_real(newsize,dims))
     
     !! Copy position
     tmp_array2_real(1:newsize,1:dims)=rp(1:newsize,1:dims)
     deallocate(rp);allocate(rp(newsize,dims))
     rp = tmp_array2_real

     !! Copy normals
     tmp_array2_real(1:newsize,1:dims)=rnorm(1:newsize,1:dims)
     deallocate(rnorm);allocate(rnorm(newsize,dims))
     rnorm = tmp_array2_real

     !! Copy s
     tmp_array_real(1:newsize)=s(1:newsize)
     deallocate(s);allocate(s(newsize))
     s = tmp_array_real

     !! Copy h
     tmp_array_real(1:newsize)=h(1:newsize)
     deallocate(h);allocate(h(newsize))
     h = tmp_array_real
     
     !! Copy fd_parent
     tmp_array_int(1:newsize)=fd_parent(1:newsize)
     deallocate(fd_parent);allocate(fd_parent(newsize))
     fd_parent = tmp_array_int
     
#ifdef dim3     
     !! Copy zlayer_index_global
     tmp_array_int(1:newsize)=zlayer_index_global(1:newsize)
     deallocate(zlayer_index_global);allocate(zlayer_index_global(newsize))
     zlayer_index_global = tmp_array_int     

     !! Copy ilayer_index
     tmp_array_int(1:newsize)=ilayer_index(1:newsize)
     deallocate(ilayer_index);allocate(ilayer_index(newsize))
     ilayer_index = tmp_array_int
#endif

     deallocate(tmp_array_real,tmp_array_int,tmp_array2_real)
     
     return
  end subroutine reduce_arrays
!! ------------------------------------------------------------------------------------------------
  subroutine load_control_data
     integer(ikind) :: dummy_int
     real(rkind) :: dummy_real
     

     !! Load data from the control file
     open(unit=12,file='control.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)

     !! Length-scale
     read(12,*)
     read(12,*) dummy_real
     read(12,*)
     
     
     
     
     close(12)
     
     write(6,*) iproc,dummy_real  
     
     return
  end subroutine load_control_data
!! ------------------------------------------------------------------------------------------------  
  subroutine load_chemistry_data
     integer(ikind) :: ispec,iorder
     
     !! Load data from the thermochemistry control file
     open(unit=12,file='thermochem.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)

     !! Number of chemical species
     read(12,*)
     read(12,*) nspec
     read(12,*)
#ifndef ms
     nspec = 1  !! If not multispecies, nspec must be 1
#endif          

     !! Number of coefs and order of polynomial for cp(T)
     read(12,*)
     read(12,*) ncoefs_cp
     read(12,*)
     !! ncoefs = (polyorder + 1) + 1 + 1: terms in polynomial,coef for h, coef for s
     polyorder_cp = ncoefs_cp - 3
     
     !! Temperature limits
     read(12,*)
     read(12,*) T_low,T_high
     read(12,*)

     !! Allocate space for molar mass, Lewis number, and polynomial fitting for cp(T)   
     allocate(molar_mass(nspec),one_over_Lewis_number(nspec))
     allocate(coef_cp(nspec,ncoefs_cp),coef_h(nspec,ncoefs_cp),coef_dcpdT(nspec,ncoefs_cp))
          
     !! Load molar mass, Lewis, and polynomial fits   
     read(12,*) !! Read comment line  
     do ispec = 1,nspec
        read(12,*) !! Read species identifier comment line
        read(12,*) molar_mass(ispec)
        read(12,*) one_over_Lewis_number(ispec)
        one_over_Lewis_number(ispec) = one/one_over_Lewis_number(ispec)

        read(12,*)  !! Comment line
        do iorder=1,ncoefs_cp
           read(12,*) coef_cp(ispec,iorder)
        end do
        read(12,*) !! Blank line                            
                           
        !! Convert cp coefficients from molar to mass based
        coef_cp(ispec,:) = coef_cp(ispec,:)*Rgas_universal/molar_mass(ispec)

        !! Pre-divide coefs by iorder for h.
        do iorder = 1,polyorder_cp + 1
           coef_h(ispec,iorder) = coef_cp(ispec,iorder)/dble(iorder)
        end do
        coef_h(ispec,polyorder_cp+2) = coef_cp(ispec,polyorder_cp+2)
        
        !! Pre-multiply coefs by iorder-1 for dcp/dT
        do iorder = 1,polyorder_cp + 1
           coef_dcpdT(ispec,iorder) = coef_cp(ispec,iorder)*dble(iorder-1)
        end do
                
     end do     

     read(12,*) !! Read comment line    
     !! Number of steps
     read(12,*)
     read(12,*) nsteps
     read(12,*)
     
     !! Space for rate constants
     allocate(Arrhenius_coefs(nsteps,3))
     
     read(12,*) !! Blank lines TBC
     read(12,*)
     read(12,*)     
     read(12,*)
     read(12,*) Arrhenius_coefs(1,1:3)

     !! Take logarithm of pre-exponential factor
     Arrhenius_coefs(1,1) = log(Arrhenius_coefs(1,1))
     Arrhenius_coefs(1,3) = Arrhenius_coefs(1,3)/(Rgas_universal)
     
     close(12)               

     
     return
  end subroutine load_chemistry_data    
!! ------------------------------------------------------------------------------------------------
  subroutine make_1d_1step_flame
     integer(ikind) :: i,ispec,j
     real(rkind) :: flame_location,flame_thickness
     real(rkind) :: T_reactants,T_products
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z

     !! Position and scale     
     flame_location = zero
     flame_thickness = 5.0d-4/L_char !! Scale thickness because position vectors are scaled...

     !! Temperatures
     T_reactants = 3.0d2
     T_products = 2.3d3
     
     !! Pressure through flame
     P_flame = rho_char*Rgas_universal*T_reactants/molar_mass(1) 
     
     !! Inflow speed
     u_reactants = u_inflow
     
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Error function based progress variable
        c = half*(one + erf((x-flame_location)/flame_thickness))
        
        !! Temperature profile
        T(i) = T_reactants + (T_products - T_reactants)*c
        
        !! Composition
        Yspec(i,1) = one - c
        Yspec(i,2) = c
        
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)/molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal
        
        !! Density
        lnro(i) = log(P_flame) - log(Rmix_local) - log(T(i))
        
        !! Velocity
        u(i) = u_reactants*rho_char/exp(lnro(i))
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              T_bound(j) = T(i)
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
              T_bound(j) = T(i) !! Inflow temperature is T_cold
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              T_bound(j) = T(i)  !! Outflow temperature is T_hot
           end if
        end do
     end if
 
  
  
     return
  end subroutine make_1d_1step_flame
!! ------------------------------------------------------------------------------------------------
  subroutine load_flame_file
     use thermodynamics
     integer(ikind) :: i,ispec,j,nflamein
     real(rkind) :: flame_location,flame_thickness
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z
     real(rkind),dimension(:),allocatable :: flamein_ro,flamein_u,flamein_T,flamein_Y,flamein_roE
     real(rkind),dimension(:),allocatable :: flamein_p
     real(rkind) :: dx_flamein

     !! Open a file containing flame profile
     open(unit=19,file='init_flame.in')
     
     !! Specify the flame pressure
     P_flame = 1.0d5

     !! Read file, then close it.     
     !! Format of the file is ro,u,T,Y_reactants over 1001 evenly spaced steps covering 
     !! a 0.01m long domain.
     nflamein = 1001
     dx_flamein = 0.01/(nflamein-1)
     allocate(flamein_ro(nflamein),flamein_u(nflamein),flamein_T(nflamein),flamein_Y(nflamein))
     allocate(flamein_roE(nflamein),flamein_p(nflamein))
     do i=1,nflamein
        read(19,*) flamein_ro(i),flamein_u(i),flamein_T(i),flamein_Y(i), &
                   flamein_roE(i),flamein_p(i)
     end do
     close(19)
     
       
     !! Loop through all particles. Find the "cell" the particle resides in. Copy data.     
     !$omp parallel do private(j,x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = (rp(i,1)+half)*L_char  !! Scale x - this requires L_char to match the length
        !! also needs shifting depending on set up.
        
        !! Nearest index in flame-in data
        j = floor(x/dx_flamein) + 1
        
        !! Copy data
        lnro(i) = log(flamein_ro(j))
        u(i) = flamein_u(j)
        v(i) = zero
        w(i) = zero
        T(i) = flamein_T(j)
        Yspec(i,1) = flamein_Y(j)
        Yspec(i,2) = one - Yspec(i,1)    
        roE(i) = flamein_roE(j)   
        p(i) = flamein_p(j)   
        
     end do
     !$omp end parallel do
     
     !! Free up space
     deallocate(flamein_ro,flamein_u,flamein_T,flamein_Y,flamein_roE,flamein_p)
     
     !! Temporarily copy some energy data to halos and mirrors (it will be later overwritten)
     !$omp parallel do
     do i=npfb+1,np
        roE(i) = roE(1)
     end do
     !$omp end parallel do
     
     !! Re-evaluate temperature from energy.  
     call evaluate_temperature_and_pressure

     
     !! Re-evaluate density from T and read-from-file-Pressure
!     !$omp parallel do 
!     do i=1,npfb
!        lnro(i) = log(p(i)/(Rgas_mix(i)*T(i)))
!     end do
!     !$omp end parallel do
                
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              T_bound(j) = T(i)
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
              T_bound(j) = T(i) !! Inflow temperature is T_cold
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              T_bound(j) = T(i)  !! Outflow temperature is T_hot
           end if
        end do
     end if
 
  
  
     return
  end subroutine load_flame_file
!! ------------------------------------------------------------------------------------------------
  subroutine hardcode_initial_conditions
     !! Temporary routine to generate initial conditions from some hard-coded functions.
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro
     
     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(x,y,z,tmp,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        u(i) = u_char!-cos(two*pi*x)*sin(two*pi*y)!*cos(two*pi*z/Lz)!*oosqrt2
        v(i) = zero!sin(two*pi*x)*cos(two*pi*y)!*cos(two*pi*z/Lz)    !!c c
        w(i) = zero!u(i);u(i)=zero
!        tmp = -half*half*(cos(4.0d0*pi*x) + cos(4.0d0*pi*y))/csq  !! Modify for not(isoT)
        lnro(i) = log(rho_char)!log(rho_char + tmp)
!if(x.le.zero) lnro(i) = log(1.2*rho_char)             

        tmp = half*(one+erf(5.0d0*x))
        T(i) = T_ref    
        tmp = rho_char*T_ref/T(i)
        lnro(i) = log(tmp)                                            
              
#ifdef ms    
!        do ispec=1,nspec      
           tmp = one - half*(one + erf(5.0d0*x))
           Yspec(i,1) = tmp
           Yspec(i,2) = one - tmp

!        end do
#endif         
                    
     end do
     !$OMP END PARALLEL DO
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero

              tmp = T_ref*(one + 0.01*sin(two*pi*rp(i,3)/Lz))
              T_bound(j) = tmp !! Might want changing                
              T(i) = T_bound(j)
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              u(i)=u_char
           end if
        end do
     end if   
  
     return
  end subroutine hardcode_initial_conditions  
!! ------------------------------------------------------------------------------------------------
  subroutine load_restart_file(n_restart)
     !! Load initial conditions from a dump file
     integer(ikind),intent(in) :: n_restart
     integer(ikind) :: k,i,j
     real(rkind) :: tmp,tmpro
     character(70) :: fname  

#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     
     if( n_restart .lt. 10 ) then 
        write(fname,'(A17,I5,A1,I1)') './data_out/layer_',k,'_',n_restart
     else if( n_restart .lt. 100 ) then 
        write(fname,'(A17,I5,A1,I2)') './data_out/layer_',k,'_',n_restart        
     else if( n_restart .lt. 1000 ) then
        write(fname,'(A17,I5,A1,I3)') './data_out/layer_',k,'_',n_restart        
     else
        write(fname,'(A17,I5,A1,I4)') './data_out/layer_',k,'_',n_restart        
     end if 
     !! Open the file
     open(14,file=fname)
     read(14,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart"
     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(14,*) tmp,tmp,tmp,tmp,h(i),k,tmpro,u(i),v(i),w(i),tmp,T(i),Yspec(i,1:nspec)
#else
        read(14,*) tmp,tmp,tmp,h(i),k,tmpro,u(i),v(i),tmp,T(i),Yspec(i,1:nspec)
#endif        
        if(k.ne.node_type(i)) then
           write(6,*) "ERROR: Problem in restart file. STOPPING."
#ifdef mp
           call MPI_Abort(MPI_COMM_WORLD, k, ierror)
#else
           stop
#endif
        end if
        lnro(i) = log(tmpro)
     end do
     
     !! Re-specify the boundary temperatures
     if(nb.ne.0) then
        do j=1,nb
           i=internal_list(j)
           T_bound(j) = T(i)
        end do
     end if
      
     close(14)
  
     return
  end subroutine load_restart_file
!! ------------------------------------------------------------------------------------------------
end module setup
