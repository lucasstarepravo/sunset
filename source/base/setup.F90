module setup
  !! This module contains routines to read in node distributions, modify them with shifting 
  !! pre-process to obtain boundary normal vectors, and write field data out to file.
  use kind_parameters
  use common_parameter
  use common_2d
  use omp_lib
  use neighbours
  use output
#ifdef mp
  use mpi
  use mpi_transfers
#endif    
  implicit none
  
  real(rkind),dimension(:,:),allocatable :: rndshift_x,rndshift_y 
  real(ikind) :: n_wvnmbrs
  real(rkind) :: wn_max,wn_min  

contains
!! ------------------------------------------------------------------------------------------------
  subroutine initial_setup  
     !! Initialises some key simulation parameters

#ifdef mp
     call MPI_COMM_SIZE(MPI_COMM_WORLD, nprocs, ierror)
     call MPI_COMM_RANK(MPI_COMM_WORLD, iproc, ierror)
#endif

     !! Time begins at zero
     time = zero;itime=0
     dt_out = 0.1d0         !! Frequency to output fields
     dt_mout = 0.01d0       !! Frequency to output mean stats
     time_end = 1.0d2
  
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
     open(unit=198,file='data_out/statistics/intenergy_mean.out')
  
     !! Profiling:
     segment_time_local = zero
     cputimecheck = zero
  
  end subroutine initial_setup
!! ------------------------------------------------------------------------------------------------
  subroutine setup_domain
     !! Reads in boundary patches, builds domain, calls decomposition and 
     !! boundary setup routines
     use boundaries
     integer(ikind) i,j,ii,jj,npfb_tmp,k,dummy_int,dummy_int2,nblock,nblock0,nm
     real(rkind) :: ns,dummy,prox,rad,radmin,xmin_local,xmax_local,x,y,xn,yn
     real(rkind),dimension(dims) :: rij
     real(rkind),dimension(:,:),allocatable :: tmp_vec
     integer(ikind) :: nl_ini,nl_end,nl_iniC,nl_endC,nl_ini2,nl_end2

     !! STEP 1: Load IPART (some params, plus list of nodes + boundary normals)
     !! =======================================================================
     open(13,file='IPART')
     read(13,*) nb,npfb,dummy      !! dummy is largest s(i) in domain...
     read(13,*) xmin,xmax,ymin,ymax
     read(13,*) xbcond,ybcond
     !! Calculate some useful constants
     h0 = hovs*dummy;sup_size = ss*h0;h2=h0*h0;h3=h2*h0
         
#ifdef mp

     read(13,*) nprocsX,nprocsY

     !! Domain decomposition: How many processors in X and Y 
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
     YU_thisproc = maxval(rp(10+floor(0.5*npfb):npfb,2))  !! Adjustments for cyclical columns. 
     YD_thisproc = minval(rp(1:floor(0.5*npfb)-10,2))       !! Should be oK given node ordering
     
    
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
     npfb_layer = npfb;nb_layer=nb;nb_layer_global = nb_global;npfb_layer_global = npfb_global  
     dz = h0/hovs  
     nz=1     
#endif     

     !! STEP 2: build mirrors, halos and neighbours
     !! =======================================================================
     call create_mirror_particles

#ifdef mp
     !! Initial halo build - much too big (based on kernel size x 1.2?)
     call build_halos

#ifdef dim3
     !! Transfer information about z-layer
     call halo_exchange_int(zlayer_index)
#endif     
           
     write(6,*) "Proc",iproc,"with",nb,npfb,np_nohalo,np
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)    

     !! Refine halos: compile lists of halo nodes which are used, and send them back to originator
     !! processors, then re-build halos
     call find_neighbours

     call refine_halos

     deallocate(ij_link,ij_count)
     
     !! Transfer discretisation properties   
     call halo_exchange(h)
     call halo_exchange(s)
     call halo_exchange_int(node_type)     
#ifdef dim3
     !! Transfer information about z-layer
     call halo_exchange_int(zlayer_index)
#endif
     
#else
     np_nohalo = np  
#endif     

     !! Allocate arrays for properties
     allocate(u(np),v(np),w(np),lnro(np),roE(np),Y0(np))
     u=zero;v=zero;w=zero;lnro=zero;roE=one;Y0=zero

#ifdef mp
     call halo_exchanges_all
     call halo_exchange_int(node_type)     
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
     use boundaries
     !! Temporary subroutine whilst building mfcomp. Initialises all fields
     integer(ikind) :: i,j,k,n_restart
     real(rkind) :: x,y,z,tmp,tmpro
     character(70) :: fname
     
     call initialise_grf
 

#ifndef restart     
     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(x,y,z,tmp)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        u(i) = u_char!-cos(two*pi*x)*sin(two*pi*y)!*oosqrt2
        v(i) = zero!sin(two*pi*x)*cos(two*pi*y)!*cos(two*pi*z/dble(nz)/dz)    !!c c
        w(i) = zero!u(i);u(i)=zero
        tmp = -half*half*(cos(4.0d0*pi*x) + cos(4.0d0*pi*y))/csq
        lnro(i) = log(rho_char)!log(rho_char + tmp)!log(rho_char)       
        roE(i) = exp(lnro(i))*(T0*Rs0/gammagasm1 + half*u(i)*u(i) + half*v(i)*v(i) + half*w(i)*w(i))
        
!        call evaluate_grf(x,y,z,tmp)
!        if(tmp.le.zero) then 
!           lnro(i) = log(rho_char)
!        else
!           lnro(i) = log(onethird*rho_char)
!        end if

        
        !! Hydrostatic energy gradient 
!        tmp = T0*(one + 0.01*sin(two*pi*z/dble(nz)/dz))
!        tmp = (rho_char/exp(lnro(i)))*tmp*Rs0/gammagasm1 + dot_product(grav,rp(i,:))/gammagasm1
!        roE(i) = (tmp + half*u(i)*u(i) + half*v(i)*v(i) + half*w(i)*w(i))*exp(lnro(i))
              
        !! Lamb-Oseen vortex
!        tmp = half ! vortex strength, vortex size is 0.1 (so 2R^2 = 0.02)
!        u(i) = u_char + tmp*two*y*exp(-(x*x+y*y)/0.02d0)/0.02
!        v(i) = -tmp*two*x*exp(-(x*x+y*y)/0.02d0)/0.02
!        tmp = -tmp*tmp*exp(-(x*x+y*y)/0.01d0)/0.01/csq
!        lnro(i) = log(one+tmp)         
!         Y0(i) = exp(-100.0d0*((x+0.3d0)**two + y*y))!verysmall
!         Y0(i) = half*erf(-200.0d0*y) + half
     end do
     !$OMP END PARALLEL DO
     
     !! Values on boundaries
     if(nb_global.ne.0)then
        if(nb.ne.0)then
           allocate(T_bound(nb));T_bound = T0
           do j=1,nb
              i=boundary_list(j)
              if(node_type(i).eq.0) then !! wall initial conditions
                 u(i)=zero;v(i)=zero;w(i)=zero

                 tmp = T0*(one + 0.01*sin(two*pi*rp(i,3)/dble(nz)/dz))
                 T_bound(j) = tmp*rho_char/exp(lnro(i)) + dot_product(grav,rp(i,:))/Rs0                
!                 T_bound(j) = (rho_char/exp(lnro(i)))*T0
                 tmp = T_bound(j)*Rs0/gammagasm1
                 roE(i) = tmp*exp(lnro(i)) 
              end if                 
              if(node_type(i).eq.1) then !! inflow initial conditions
                 u(i)=u_char
              end if
              if(node_type(i).eq.2) then !! outflow initial conditions
                 u(i)=u_char
              end if
           end do
        end if
     end if
#else    
     !! RESTART OPTION. Ask for input number (just hard-coded for now...)
     n_restart = 200
#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     
     if( n_restart .lt. 10 ) then 
        write(fname,'(A14,I5,A1,I1)') './data_out/uv_',k,'_',n_restart
     else if( n_restart .lt. 100 ) then 
        write(fname,'(A14,I5,A1,I2)') './data_out/uv_',k,'_',n_restart        
     else if( n_restart .lt. 1000 ) then
        write(fname,'(A14,I5,A1,I3)') './data_out/uv_',k,'_',n_restart        
     else
        write(fname,'(A14,I5,A1,I4)') './data_out/uv_',k,'_',n_restart        
     end if 
     !! Open the file
     open(14,file=fname)
     read(14,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart"
     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(14,*) tmp,tmp,tmp,h(i),k,tmpro,u(i),v(i),w(i),tmp,roE(i),tmp
#else
        read(14,*) tmp,tmp,h(i),k,tmpro,u(i),v(i),tmp,roE(i),tmp
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
        roE(i) = roE(i)*tmpro
     end do
     close(14)
#endif
   
     !! Mirrors and halos                        
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchanges_all
#endif          
     
     
     !! Set the initial forcing to zero. It will only be changed if PID controller in velcheck is used.
     driving_force(:) = zero
     
     return
  end subroutine initial_solution   
!! ------------------------------------------------------------------------------------------------
  subroutine build_3rd_dimension
     integer(ikind) :: nm,i,k,iz,ilayerm1
     real(rkind) :: dz_local
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

     !! Z dimension = L_char
!     nz = ceiling(L_char/dz)
!     dz = L_char/dble(nz)
     
     nz = 10
     
     !! Minimum number in z is 6
     if(nz.lt.6) then
        nz = 6
        dz = L_char/6.0d0
     end if
     
     write(6,*) dz,nz         
     
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
     npfb_layer = npfb;nb_layer=nb;nb_layer_global = nb_global;npfb_layer_global = npfb_global

     !! New sizes (local and global)
     npfb = npfb*nz;nb = nb*nz;nb_global = nb_global*nz;npfb_global = npfb_global*nz
     
     !! Deallocate and reallocate arrays
     nm = 10
     deallocate(rp,rnorm,h,s,node_type,fd_parent)
     allocate(rp(nm*npfb,dims),rnorm(nm*npfb,dims),h(nm*npfb),s(nm*npfb));rp=zero
     allocate(node_type(nm*npfb));node_type=0
     allocate(fd_parent(nm*npfb));fd_parent=0
     allocate(zlayer_index(nm*npfb))
     
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
           zlayer_index(k) = iz
           rp(k,3) = dble(iz-1)*dz
        end do
     end do
        
        
     !! Deallocate temporary arrays
     deallocate(rptmp,rnormtmp,htmp,stmp,fd_parenttmp,node_typetmp)    
     
  
     return
  end subroutine build_3rd_dimension  
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_grf
     integer(ikind) :: i,j,k
     
     n_wvnmbrs = 2
     allocate(rndshift_x(n_wvnmbrs,n_wvnmbrs))
     allocate(rndshift_y(n_wvnmbrs,n_wvnmbrs))     
     
     do i=1,n_wvnmbrs
        do j=1,n_wvnmbrs
           rndshift_x(i,j) = rand()*two*pi
           rndshift_y(i,j) = rand()*two*pi           
        end do
     end do
     
     wn_max = 5.0d0*two*pi/(xmax-xmin)
     wn_min = 1.0d0*two*pi/(xmax-xmin)
          
     return
  end subroutine initialise_grf
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_grf(x,y,z,grf)
     real(rkind),intent(in) :: x,y,z
     real(rkind),intent(out) :: grf
     integer(ikind) :: i,j,k
     real(rkind) :: wn_i,wn_j
         
     grf = zero
     do i=1,n_wvnmbrs
        wn_i = dble(i)*two*pi/(xmax-xmin)!wn_min + dble(i)*(wn_max-wn_min)/dble(n_wvnmbrs)
        do j=1,n_wvnmbrs
           wn_j = dble(j)*two*pi/(xmax-xmin)!wn_min + dble(j)*(wn_max-wn_min)/dble(n_wvnmbrs)
           grf = grf + cos(wn_i*x + rndshift_x(i,j))* &
                       cos(wn_j*y + rndshift_y(i,j))
        end do
     end do
  
  
     return
  end subroutine evaluate_grf  
!! ------------------------------------------------------------------------------------------------
end module setup
