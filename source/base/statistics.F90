module statistics
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2022 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to perform reduction calculations to obtain and output global 
  !! statistics, and write data to standard-out.
  use kind_parameters
  use common_parameter
  use common_vars
  use omp_lib
  use neighbours
#ifdef mp
  use mpi
  use mpi_transfers
#endif    
  implicit none

contains
!! ------------------------------------------------------------------------------------------------
  subroutine open_stats_files
     !! Opens units for statistic outputting files

     !! Main time out-file, contains time, dt, etc.
     open(unit=21,file='./data_out/time.out')
     
     !! Total cpu time per step
     open(unit=191,file='data_out/statistics/cputime.out')
     
     !! Time, value of time-step, value of time-step divided by CFL-based timestep (if using PID),
     !! and maximum sound speed in domain
     open(unit=192,file='data_out/statistics/dt.out')
     
     !! Time, total mass within the domain
     open(unit=193,file='data_out/statistics/masscheck.out')
     
     !! Time, net forces on wall boundaries. Needs modifying on case-by-case basis
     open(unit=194,file='data_out/statistics/liftdrag.out')
     
     !! Time, mean (or L2) velocity, and components thereof
     open(unit=195,file='data_out/statistics/velcheck.out')
     
     !! Time, L2 error of velocity field relative to analytic solution for Taylor-Green flow  
     open(unit=196,file='data_out/statistics/taylor_green_l2.out')
     
     !! Time, total energy within domain.
     open(unit=197,file='data_out/statistics/energy_sum.out')  
     
     !! Time, enstrophy
     open(unit=198,file='data_out/statistics/enstrophy.out')

     !! Time, L2 error in mass fraction sum, total mass of each species in domain
     open(unit=199,file='data_out/statistics/species.out')
     
     !! Time, heat release rate
     open(unit=200,file='data_out/statistics/heat_release_rate.out')

     return     
  end subroutine open_stats_files
!! ------------------------------------------------------------------------------------------------
  subroutine statistics_control
     !! This routine controls the statistics calculation and output routines.
     integer(ikind) :: i
     real(rkind) :: tot_mass,tot_vol,tmpro,dVi,tot_mass_tmp,tot_vol_tmp,tot_roE,tot_roE_tmp
     integer(ikind),parameter :: istats_freq = 10
     
     
     !! Evaluate the mean velocity and adjust pressure gradient if required
     !! This should be done every step if using PID for pressure gradient
#ifdef pgrad     
     call velocity_control
#else
     if(itime.eq.0.or.mod(itime,istats_freq).eq.0) call velocity_control     
#endif  
     
     !! Other 
     if(itime.eq.0.or.mod(itime,istats_freq).eq.0) then    
        !! Check conservation of mass and energy
        call mass_and_energy_check
     
        !! Error evaluation for Taylor Green vortices?
!        call error_TG

        !! Calculate the lift and drag on all solid obstacles
!        call liftdrag

        !! Calculate how well balanced the MPI decomposition is
!        call check_load_balance

        !! Check the conservation of the species equations
        call species_check
        
        !! Check heat release
        call heat_release_rate
     endif

     return
  end subroutine statistics_control  
!! ------------------------------------------------------------------------------------------------   
  subroutine check_load_balance  
     real(rkind),dimension(:),allocatable :: prof_tmp,prof_tmp_local,load_n_n
     integer(ikind),dimension(:),allocatable :: sum_n_n,sum_n_n_local
     integer(ikind) :: sum_sum_n_n
#ifdef mp     
     !! Calculate relative amounts of work being done by each processor
     allocate(prof_tmp_local(nprocs),prof_tmp(nprocs));prof_tmp_local=zero
     prof_tmp_local(iproc+1) = sum(segment_time_local(1:8))
     call MPI_ALLREDUCE(prof_tmp_local,prof_tmp,nprocs,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     
     !! Sum of the total neighbours of all nodes on each processor
     allocate(sum_n_n(nprocs),sum_n_n_local(nprocs));sum_n_n_local=0
     sum_n_n_local(iproc+1) = sum(ij_count(1:npfb))
     call MPI_ALLREDUCE(sum_n_n_local,sum_n_n,nprocs,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)
     sum_sum_n_n = sum(sum_n_n(:))
     allocate(load_n_n(nprocs))
     load_n_n = sum_n_n/dble(npfb)!dble(nprocs)*dble(sum_n_n(:))/dble(sum_sum_n_n)
     
     
     !! Output the amounts of work being done
     if(iproc.eq.0) write(6,*) dble(nprocs)*prof_tmp(:)/sum(prof_tmp(:))

     !! Output the expected load (nodes*neighbours)
     if(iproc.eq.0) write(6,*) load_n_n       
     
     deallocate(prof_tmp_local,prof_tmp)
     deallocate(sum_n_n,sum_n_n_local,load_n_n)
#endif
  
     return
  end subroutine check_load_balance   
!! ------------------------------------------------------------------------------------------------
  subroutine liftdrag
     !! TO DO: Update for 3 dimensional simulations  
     !! TO DO: Update for non-isothermal simulations
     use derivatives
#ifdef isoT     
     integer(ikind) :: i,j
     real(rkind),dimension(dims) :: gradu0,gradv0,Fn,force,force_tmp
     real(rkind),dimension(dims,dims) :: Jinv,sigma
     real(rkind) :: xn,yn
   
  
     !! Calculate the velocity gradient 
     allocate(gradu(npfb,dims),gradv(npfb,dims),gradw(npfb,dims))
     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)
     call calc_gradient(w,gradw)     

     force = zero
     !$omp parallel do private(i,Jinv,gradu0,gradv0,xn,yn,sigma,Fn) reduction(+:force)
     do j=1,nb
        i=boundary_list(j)
        if(node_type(i).eq.0)then
           xn = rnorm(i,1);yn=rnorm(i,2)
           Jinv(1,1)=xn;Jinv(1,2)=-yn;Jinv(2,1)=yn;Jinv(2,2)=xn   !! Jacobian for normal-tangent to x-y  
           Jinv(3,:)=zero;Jinv(:,3)=zero;Jinv(3,3)=one         
           gradu0(:) = matmul(Jinv,gradu(i,:))  !! Velocity gradients in x-y FoR
           gradv0(:) = matmul(Jinv,gradv(i,:))           
           
           !! Total stress on surface
           sigma(1,1) = visc(i)*(fourthirds*gradu0(1)-twothirds*gradv0(2)) - csq*(ro(i)-one)
           sigma(1,2) = visc(i)*(gradu0(2)+gradv0(1))
           sigma(2,1) = visc(i)*(gradu0(2)+gradv0(1))           
           sigma(2,2) = visc(i)*(fourthirds*gradv0(2)-twothirds*gradu0(1)) - csq*(ro(i)-one)
          
           Fn(:) = matmul(sigma,rnorm(i,:))   !! Force on surface (sigma.n)
                     
           force(:) = force(:) + Fn(:)*s(i)*L_char  !! Integrate over surface... s(i) is dimensionless node spacing
        end if
     end do
     !$omp end parallel do
     deallocate(gradu,gradv,gradw)

#ifdef mp
     force_tmp = force
     call MPI_ALLREDUCE(force_tmp,force,dims,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)     
     if(iproc.eq.0)then
        write(194,*) time/Time_char,force
        flush(194)
     end if        
#else
     write(194,*) time/Time_char,force
     flush(194)
#endif     
     
#endif     
     return
  end subroutine liftdrag
!! ------------------------------------------------------------------------------------------------  
  subroutine velocity_control
     !! Output the L2 of velocity over the domain
     integer(ikind) :: i
     real(rkind) :: tot_vel,tot_vol,tmpro,dVi,tmpvel
     real(rkind) :: tot_vel_tmp,tot_vol_tmp
     real(rkind),dimension(dims) :: tot_u,tot_u_tmp
     real(rkind) :: facA,facB,facC,facT,deflowdt
       
     tot_vel = zero
     tot_vol = zero
     tot_u = zero
     !$omp parallel do private(tmpro,tmpvel,dVi) reduction(+:tot_vel,tot_vol,tot_u)
     do i=1,npfb
        tmpro = ro(i)
        dVi = s(i)*s(i) !! assume square nodes for now...
#ifdef dim3
        dVi = dVi*dz
#endif               
        tmpvel = (u(i)*u(i)+v(i)*v(i)+w(i)*w(i))
        tot_vel = tot_vel + tmpvel*dVi
        tot_vol = tot_vol + dVi
        tot_u = tot_u + (/u(i),v(i),w(i)/)*dVi
     end do
     !$omp end parallel do

#ifdef mp
     tot_vel_tmp = tot_vel;tot_vol_tmp = tot_vol;tot_u_tmp = tot_u
     call MPI_ALLREDUCE(tot_vel_tmp,tot_vel,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_vol_tmp,tot_vol,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_u_tmp,tot_u,3,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
#endif         
    
     !! Normalise over volume
     tot_vel = sqrt(tot_vel/tot_vol)
     tot_u(:) = tot_u(:)/tot_vol

     !! If we want to P.I.D. control over the mean velocity
#ifdef pgrad     
     !! New error     
     eflow_n = u_char - tot_u(1)!tot_vel
          
     !! Integral term
     sum_eflow = sum_eflow + eflow_n*dt
     
     !! Derivative term
     deflowdt = (eflow_n-eflow_nm1)/dt
    
     !! P, I and D factors..  
     facA = two*one
     facB = facA/1.0d-1
     facC = facA*0.02d0
         
     driving_force(1) = facA*eflow_n + facB*sum_eflow + facC*deflowdt
     !! Impose some upper and lower limits
     driving_force(1) = min(5.0d0,driving_force(1))
     driving_force(1) = max(-5.0d0,driving_force(1))
                       
     !! Pass new eflow to old eflow
     eflow_nm1=eflow_n
#else
     driving_force = zero
#endif
      
#ifdef mp
     if(iproc.eq.0)then
        write(195,*) time/Time_char,tot_vel,tot_u,driving_force(1)
        flush(195)
     end if
#else
     write(195,*) time/Time_char,tot_vel,tot_u,driving_force(1)    
     flush(195)
#endif       
     
     
     return
  end subroutine velocity_control   
!! ------------------------------------------------------------------------------------------------  
  subroutine heat_release_rate
     !! This subroutine calculates the total heat release rate integrated over the domain
     integer(ikind) :: i
     real(rkind) :: tot_hrr,dVi,tot_hrr_tmp,tot_vol_tmp
#ifdef react    
   
     tot_hrr = zero
     !$omp parallel do private(dVi) reduction(+:tot_hrr)
     do i=1,npfb
        dVi = s(i)*s(i) !! assume square nodes for now...
#ifdef dim3
        dVi = dVi*dz
#endif        
        tot_hrr = tot_hrr + hrr(i)*dVi
     end do
     !$omp end parallel do
     
#ifdef mp
     tot_hrr_tmp = tot_hrr
     call MPI_ALLREDUCE(tot_hrr_tmp,tot_hrr,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
#endif     

     !! Scale by length-scale squared or cubed
#ifdef dim3
     tot_hrr = tot_hrr*L_char**3.0d0
#else
     tot_hrr = tot_hrr*L_char**two
#endif     
        
#ifdef mp
     if(iproc.eq.0)then
        write(200,*) time/Time_char,tot_hrr
        flush(200)
     end if
#else
     write(200,*) time/Time_char,tot_hrr
     flush(200)
#endif

#endif

     return
  end subroutine heat_release_rate 
!! ------------------------------------------------------------------------------------------------
  subroutine mass_and_energy_check
     !! This subroutine calculates the total mass and total energy in the domain.
     integer(ikind) :: i
     real(rkind) :: tot_mass,tot_vol,tmpro,dVi,tot_mass_tmp,tot_vol_tmp,tot_roE,tot_roE_tmp
    
   
     tot_mass = zero;tot_vol = zero;tot_roE = zero
     !$omp parallel do private(tmpro,dVi) reduction(+:tot_mass,tot_vol,tot_roE)
     do i=1,npfb
        tmpro = ro(i)
        dVi = s(i)*s(i) !! assume square nodes for now...
#ifdef dim3
        dVi = dVi*dz
#endif        
        tot_mass = tot_mass + tmpro*dVi
        tot_vol = tot_vol + dVi
#ifndef isoT
        tot_roE = tot_roE + dVi*roE(i)
#endif        
     end do
     !$omp end parallel do
     
     !! Scale to make dimensional
#ifdef dim3
     tot_mass = tot_mass*L_char**3.0d0
     tot_vol = tot_vol*L_char**3.0d0
     tot_roE = tot_roE*L_char**3.0d0          
#else
     tot_mass = tot_mass*L_char**two
     tot_vol = tot_vol*L_char**two
     tot_roE = tot_roE*L_char**two          
#endif     
     
#ifdef mp
     tot_mass_tmp = tot_mass;tot_vol_tmp = tot_vol;tot_roE_tmp = tot_roE
     call MPI_ALLREDUCE(tot_mass_tmp,tot_mass,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_vol_tmp,tot_vol,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_roE_tmp,tot_roE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
#endif     
              
#ifdef mp
     if(iproc.eq.0)then
        write(193,*) time/Time_char,tot_mass,tot_vol
        flush(193)
        write(197,*) time/Time_char,tot_roE
        flush(197)
     end if
#else
     write(193,*) time/Time_char,tot_mass,tot_vol
     flush(193)
     write(197,*) time/Time_char,tot_roE
     flush(197)
#endif

     return
  end subroutine mass_and_energy_check 
!! ------------------------------------------------------------------------------------------------   
  subroutine species_check
     !! This subroutine calculates total quantity of each species in the domain
     integer(ikind) :: i,ispec
     real(rkind) :: dVi,tmpY,sumY,tot_error,tot_error_tmp,tmpro,tot_vol,tot_vol_tmp
     real(rkind),dimension(:),allocatable :: tot_Yspec,tot_Yspec_tmp
#ifdef ms     
     
     !! Allocate and zero accumulators
     allocate(tot_Yspec(nspec),tot_Yspec_tmp(nspec))   
     tot_Yspec = zero     
     tot_error = zero
     tot_vol=zero
     !$omp parallel do private(tmpro,dVi,ispec,tmpY,sumY) reduction(+:tot_Yspec,tot_error,tot_vol)
     do i=1,npfb
        !! Extract one/ro
        tmpro = one/ro(i)
     
        !! size of volume element
        dVi = s(i)*s(i) !! assume square nodes for now...
#ifdef dim3
        dVi = dVi*dz
#endif        

        !! Sum over all species
        sumY = -one
        do ispec=1,nspec
           
           tmpY = Yspec(i,ispec)
           
           !! Total mass of species in domain - integral of roY dV
           tot_Yspec(ispec) = tot_Yspec(ispec) + dVi*tmpY
           
           !! 1- sum_over_species(Y)
           sumY = sumY + tmpY*tmpro
        end do
        
        !! Volume integral
        tot_vol = tot_vol + dVi
        
        !! L2 error in sumY
        tot_error = tot_error + sumY*sumY*dVi
     end do
     !$omp end parallel do
     
     !! Rescale integral quantities to be dimensional
#ifdef dim3
     tot_error = tot_error*L_char**3.0d0
     tot_Yspec = tot_Yspec*L_char**3.0d0
#else
     tot_error = tot_error*L_char**two
     tot_Yspec = tot_Yspec*L_char**two
#endif     
     
#ifdef mp
     tot_Yspec_tmp = tot_Yspec
     tot_error_tmp = tot_error
     tot_vol_tmp = tot_vol
     call MPI_ALLREDUCE(tot_Yspec_tmp,tot_Yspec,nspec,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_error_tmp,tot_error,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_vol_tmp,tot_vol,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)     
     
     !! L2 norm for error           
     tot_error = sqrt(tot_error/tot_vol)
          
     if(iproc.eq.0)then
        write(199,*) time/Time_char,tot_error,tot_Yspec(:)
        flush(199)
     end if
#else
     !! Normalise mass of species by volume
     write(199,*) time/Time_char,sqrt(tot_error/tot_vol),tot_Yspec(:)
     flush(199)
#endif

#endif

     return
  end subroutine species_check    
!! ------------------------------------------------------------------------------------------------    
  subroutine poiseuille_l2norm
     !! Output the L2norm of the velocity compared with Poiseuille flow
     !! Umax=1, unit width domain...
     !! N.B. This routine needs updating for multi-processor simulations
     integer(ikind) :: i,j,jj
     real(rkind) :: y,uexact,local_error,sum_e,sum_exact,dVi
     real(rkind) :: N,X1,X2,y1
       
     sum_e = zero
     sum_exact = zero
     !$omp parallel do private(y,uexact,local_error,dVi,j,jj,N,X1,X2,y1) &
     !$omp reduction(+:sum_e,sum_exact)
     do i=1,npfb
        y = rp(i,2)
        uexact = (half-y)*(half+y)*4.0d0
        y1 = y+half


        jj = 30!floor(10.0d0/sqrt(time+0.1d0)); % vary the degree of expansion to avoid NaNs...
        do j=1,jj
           N=(2*j-1)*pi
           X1 = sin(N*y1)/N**3.0d0
           X2 = exp(-N*N*visc(i)*time/Time_char)
           uexact = uexact - 32.0d0*X1*X2
        end do       
        
        dVi = h(i)*h(i) !! assume square nodes for now...
        local_error = u(i)-uexact
        sum_e = sum_e + local_error**two
        sum_exact = sum_exact + uexact**two
     end do
     !$omp end parallel do
     sum_e = dsqrt(sum_e/sum_exact)

     write(196,*) time/Time_char,sum_e
     flush(196)  
     return
  end subroutine poiseuille_l2norm  
!! ------------------------------------------------------------------------------------------------
  subroutine error_TG
    !! Error relative to analytic solution for 2D Taylor-Green vortex decay
    implicit none
    integer(ikind) :: i
    real(rkind) :: u_exact,v_exact,x,y,expo
    real(rkind) :: U_ex,U_num,error_sum,L2error,Ex_sum,error_sum_local,ex_sum_local
  
    error_sum = zero;Ex_sum =zero
    !$omp parallel do private(x,y,u_exact,v_exact,U_ex,U_num) reduction(+:Ex_sum,error_sum)
    do i=1,npfb
       x=rp(i,1)
       y=rp(i,2)
       expo = exp(-8.0d0*pi*pi*time/time_char/200.0d0)   !! Hard-coded for Re=200. Need to modify.
       u_exact = -expo*cos(2.0*pi*x)*sin(2.0*pi*y)
       v_exact = expo*sin(2.0*pi*x)*cos(2.0*pi*y)
       U_ex = dsqrt(u_exact**2. + v_exact**2.)
       U_num = dsqrt(u(i)**2. + v(i)**2. + w(i)**2.)
     
       error_sum = error_sum + (U_num-U_ex)**2.
       Ex_sum = Ex_sum + (U_ex)**2.
    end do
    !$omp end parallel do
#ifdef mp
    error_sum_local = error_sum;ex_sum_local = ex_sum
    call MPI_ALLREDUCE(error_sum_local,error_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
    call MPI_ALLREDUCE(ex_sum_local,ex_sum,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
    L2error = dsqrt(error_sum/Ex_sum)
    if(iproc.eq.0)then
!       write(6,*) time,L2error,expo,maxval(u(1:npfb))
       write(196,*) time,L2error
       flush(196) 
    end if
#else   
    L2error = dsqrt(error_sum/Ex_sum)
!    write(6,*) time,L2error,expo,maxval(u(1:npfb))
    write(196,*) time/Time_char,L2error
    flush(196) 
#endif    
  
  
    return
  end subroutine error_TG   
!! ------------------------------------------------------------------------------------------------    
  subroutine check_enstrophy
     !! Evaluate volume averaged enstrophy and also the volume averaged kinetic energy
     !! This routine is designed specifically for 3D Taylor-Green vortex tests.
     integer(ikind) :: i
     real(rkind) :: srtnorm,sum_enstrophy,vol_i,sum_vol,sum_ke
     real(rkind) :: sum_enstrophy_local,sum_vol_local,sum_ke_local
     
     sum_enstrophy = zero
     sum_vol = zero
     sum_ke = zero
     !$omp parallel do private(srtnorm,vol_i) reduction(+:sum_enstrophy,sum_vol,sum_ke)
     do i=1,npfb
     
        vol_i = s(i)*s(i)
#ifdef dim3
        vol_i = vol_i*dz
#endif        
        srtnorm = dot_product(gradu(i,:),gradu(i,:)) + &
                  dot_product(gradv(i,:),gradv(i,:)) + &
                  dot_product(gradw(i,:),gradw(i,:))
        sum_enstrophy = sum_enstrophy + visc(i)*srtnorm*vol_i
        sum_vol = sum_vol + vol_i
        
        !!
        sum_ke = sum_ke + ro(i)*(u(i)**two + v(i)**two + w(i)**two)*vol_i
     end do
     !$omp end parallel do
     
#ifdef mp
     sum_enstrophy_local = sum_enstrophy
     sum_vol_local = sum_vol
     sum_ke_local = sum_ke
     call MPI_ALLREDUCE(sum_enstrophy_local,sum_enstrophy,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(sum_vol_local,sum_vol,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(sum_ke_local,sum_ke,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)  
     if(iproc.eq.0)then
        write(198,*) time/Time_char,sum_enstrophy/sum_vol,sum_ke/sum_vol
        flush(198)
     end if        
#else
     write(198,*) time/Time_char,sum_enstrophy/sum_vol,sum_ke/sum_vol
     flush(198)
#endif     
           
     return
  end subroutine check_enstrophy  
!! ------------------------------------------------------------------------------------------------ 
end module statistics
