module step
  !! This module contains time-stepping routines, and a routine to set the time-step.
  !! Time-stepping routines start with "step_" and perform one time step. They are
  !! called from the main loop, and they themselves call routines in the rhs module.

  use kind_parameters
  use common_parameter
  use common_2d
  use boundaries
  use rhs
  use mpi_transfers
#ifdef mp  
  use mpi
#endif  
  implicit none

  private
  public step_rk4,set_tstep,step_rk3_2N,step_rk3_4S_2R, &
         set_tstep_PID
  
  !! Error norms for RK3(2)4S[2R+]C scheme
  real(rkind) :: enrm_ro,enrm_u,enrm_v,enrm_E,enrm_Y0,enrm_w

contains
!! ------------------------------------------------------------------------------------------------
  subroutine step_rk4
     !! Classical 4th order Runge Kutta (i.e. the one on Wikipedia)
     integer(ikind) :: i,k
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: u0,v0,lnro0,roE0,Y00,w0
     real(rkind),dimension(:,:),allocatable :: rhsu,rhsv,rhslnro,rhsroE,rhsY0,rhsw

     allocate(u0(npfb),v0(npfb),w0(npfb))
     allocate(lnro0(npfb),roE0(npfb),Y00(npfb))
     allocate(rhsu(npfb,4),rhsv(npfb,4),rhsw(npfb,4))
     allocate(rhslnro(npfb,4))
     allocate(rhsroE(npfb,4),rhsY0(npfb,4))
        
     !! Temporary storage of u and time
     time0=time
     !$OMP PARALLEL DO
     do i=1,npfb
        lnro0(i) = lnro(i)
        u0(i) = u(i)
        v0(i) = v(i)
        w0(i) = w(i)
        roE0(i) = roE(i)
     end do
     !$OMP END PARALLEL DO
       
     !! RK step 1
     call calc_rhs_lnro(rhslnro(:,1))
     call calc_rhs_vel(rhsu(:,1),rhsv(:,1),rhsw(:,1))
     call calc_rhs_roE(rhsroE(:,1))
     call calc_rhs_Y0(rhsY0(:,1))
     if(nb.ne.0) call calc_rhs_nscbc(rhslnro(1:nb,1),rhsu(1:nb,1),rhsv(1:nb,1),rhsw(1:nb,1),rhsroE(1:nb,1),rhsY0(1:nb,1))

     !! RK steps 2,3,4
     do k=2,4
        !! Set the time for this substep
        time = time0 + rk4_c(k)*dt
           
        !! Find the temporary u for next substep
        !$OMP PARALLEL DO
        do i=1,npfb
           lnro(i) = lnro0(i) + rk4_a(k)*rhslnro(i,k-1)*dt
           u(i) = u0(i) + rk4_a(k)*rhsu(i,k-1)*dt
           v(i) = v0(i) + rk4_a(k)*rhsv(i,k-1)*dt
           w(i) = w0(i) + rk4_a(k)*rhsw(i,k-1)*dt           
           roE(i) = roE0(i) + rk4_a(k)*rhsroE(i,k-1)*dt
           Y0(i) = Y00(i) + rk4_a(k)*rhsY0(i,k-1)*dt
        end do
        !$OMP END PARALLEL DO
        
        !! Apply BCs and update halos
        call reapply_mirror_bcs
        call halo_exchanges_all
   
        !! calc the RHS and store
        call calc_rhs_lnro(rhslnro(:,k))
        call calc_rhs_vel(rhsu(:,k),rhsv(:,k),rhsw(:,k))
        call calc_rhs_roE(rhsroE(:,k))
        call calc_rhs_Y0(rhsY0(:,k))
        if(nb.ne.0) call calc_rhs_nscbc(rhslnro(:,k),rhsu(:,k),rhsv(:,k),rhsw(:,k),rhsroE(:,k),rhsY0(:,k))        
     end do

     !! Final weighted sum of substeps
     time = time0 + dt
     !$OMP PARALLEL DO
     do i=1,npfb
        lnro(i) = lnro0(i) + dt*dot_product(rk4_b(:),rhslnro(i,:))       
        u(i) = u0(i) + dt*dot_product(rk4_b(:),rhsu(i,:))
        v(i) = v0(i) + dt*dot_product(rk4_b(:),rhsv(i,:))
        w(i) = w0(i) + dt*dot_product(rk4_b(:),rhsw(i,:))        
        roE(i) = roE0(i) + dt*dot_product(rk4_b(:),rhsroE(i,:))
        Y0(i) = Y00(i) + dt*dot_product(rk4_b(:),rhsY0(i,:))
     end do
     !$OMP END PARALLEL DO
     deallocate(rhslnro,rhsu,rhsv,rhsroE,u0,v0,lnro0,roE0,rhsY0,Y00,rhsw,w0)

     !! Apply BCs and update halos
     call reapply_mirror_bcs  
     call halo_exchanges_all
     
     !! Filter the solution
     call filter_variables

     !! Apply BCs and update halos     
     call reapply_mirror_bcs       
     call halo_exchanges_all        


     return
  end subroutine step_rk4
!! ------------------------------------------------------------------------------------------------  
  subroutine step_rk3_2N
     !! Third order low storage Runge Kutta
     integer(ikind) :: i,k
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: u_W,v_W,lnro_W,tmpu,tmpv,tmplnro,roE_W,tmproE   
     real(rkind),dimension(:),allocatable :: Y0_W,tmpY0,w_W,tmpw

     allocate(u_W(npfb),v_W(npfb),lnro_W(npfb),roE_W(npfb),w_W(npfb));u_W=zero;v_W=zero;lnro_W=zero;roE_W=zero;w_W=zero
     allocate(tmpu(npfb),tmpv(npfb),tmplnro(npfb),tmproE(npfb),tmpw(npfb))
     allocate(Y0_W(npfb),tmpY0(npfb));Y0_W=zero
        
     !! Temporary storage of time
     time0=time

     do k=1,3

        !! Calculate the RHS
        call calc_rhs_lnro(tmplnro)
        call calc_rhs_vel(tmpu,tmpv,tmpw)
        call calc_rhs_roE(tmproE)
        call calc_rhs_Y0(tmpY0)

        if(nb.ne.0) call calc_rhs_nscbc(tmplnro,tmpu,tmpv,tmpw,tmproE,tmpY0)

        !! Set the intermediate time
        time = time0 + rk3_2N_c(k)*dt        

        !! Set w_i and new u,v
        !$omp parallel do
        do i=1,npfb
           lnro_W(i) = rk3_2n_a(k)*lnro_W(i) + dt*tmplnro(i)
           u_W(i) = rk3_2n_a(k)*u_W(i) + dt*tmpu(i)
           v_W(i) = rk3_2n_a(k)*v_W(i) + dt*tmpv(i)
           w_W(i) = rk3_2n_a(k)*w_W(i) + dt*tmpw(i)           
#ifndef isoT           
           roE_W(i) = rk3_2n_a(k)*roE_W(i) + dt*tmproE(i)
#endif           
#ifdef ms
           Y0_W(i) = rk3_2n_a(k)*Y0_W(i) + dt*tmpY0(i)
#endif
           
           lnro(i) = lnro(i) + rk3_2n_b(k)*lnro_W(i)
           u(i) = u(i) + rk3_2n_b(k)*u_W(i)
           v(i) = v(i) + rk3_2n_b(k)*v_W(i)
           w(i) = w(i) + rk3_2n_b(k)*w_W(i)           
#ifndef isoT
           roE(i) = roE(i) + rk3_2n_b(k)*roE_W(i)
#endif           
#ifdef ms
           Y0(i) = Y0(i) + rk3_2n_b(k)*Y0_W(i)
#endif
        end do
        !$omp end parallel do
       
        !! Apply BCs and update halos
        call reapply_mirror_bcs        
        call halo_exchanges_all
       
     end do
     
     !! Deallocation
     deallocate(u_W,v_W,lnro_W,roE_W,tmpu,tmpv,tmplnro,tmproE,Y0_W,tmpY0,w_W,tmpw) 

     !! Filter the solution
     call filter_variables
     
     !! Apply BCs and update halos     
     call reapply_mirror_bcs
     call halo_exchanges_all
         
     !! Set the new time   
     time = time0 + dt
 
     return
  end subroutine step_rk3_2N  
!! ------------------------------------------------------------------------------------------------  
  subroutine step_rk3_4S_2R
     !! 3rd order 4step 2 register Runge Kutta
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is lnro_W,u_W,v_W. Register 2 is lnro,u,v.
     !! Register 3 is tmplnro,tmpu,tmpv (only used for RHS)
     !! Register 4 is e_acc_lnro,e_acc_u,e_acc_v - error accumulator
     integer(ikind) :: i,k
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: u_W,v_W,w_W,lnro_W,roE_W,Y0_W,tmpu,tmpv,tmplnro,tmproE,tmpY0,tmpw
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb
         
     !! Set RKa,RKb with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)

     allocate(u_W(npfb),v_W(npfb),lnro_W(npfb),roE_W(npfb),Y0_W(npfb),w_W(npfb))
     allocate(tmpu(npfb),tmpv(npfb),tmplnro(npfb),tmproE(npfb),tmpY0(npfb),tmpw(npfb))

     !! Store prim vars in register 1 (w-register)
     !$omp parallel do
     do i=1,npfb
        lnro_W(i)=lnro(i);u_W(i)=u(i);v_W(i)=v(i);w_W(i)=w(i);roE_W(i)=roE(i);Y0_W(i)=Y0(i)
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time

     do k=1,3
        !! Calculate the RHS
        call calc_rhs_lnro(tmplnro)
        call calc_rhs_vel(tmpu,tmpv,tmpw)
        call calc_rhs_roE(tmproE)
        call calc_rhs_Y0(tmpY0)
        if(nb.ne.0) call calc_rhs_nscbc(tmplnro,tmpu,tmpv,tmpw,tmproE,tmpY0)        
     
        !! Set w_i and new u,v
        !$omp parallel do
        do i=1,npfb
           !! Store next U in register 2
           lnro(i) = lnro_W(i) + RKa(k)*tmplnro(i)
           u(i) = u_W(i) + RKa(k)*tmpu(i)
           v(i) = v_W(i) + RKa(k)*tmpv(i)
           w(i) = w_W(i) + RKa(k)*tmpw(i)           
#ifndef isoT           
           roE(i) = roE_W(i) + RKa(k)*tmproE(i)
#endif
#ifdef ms           
           Y0(i) = Y0_W(i) + RKa(k)*tmpY0(i)
#endif

           !! Store next S in register 1
           lnro_W(i) = lnro_W(i) + RKb(k)*tmplnro(i)
           u_W(i) = u_W(i) + RKb(k)*tmpu(i)
           v_W(i) = v_W(i) + RKb(k)*tmpv(i) 
           w_W(i) = w_W(i) + RKb(k)*tmpw(i)            
#ifndef isoT
           roE_W(i) = roE_W(i) + RKb(k)*tmproE(i)
#endif
#ifdef ms           
           Y0_W(i) = Y0_W(i) + RKb(k)*tmpY0(i)
#endif
        end do
        !$omp end parallel do
       
        !! Apply BCs and update halos
        call reapply_mirror_bcs
        call halo_exchanges_all
        
     end do
     
     !! Final substep: returns solution straight to lnro,u,v,E (register 2)
     !! and doesn't update S
     call calc_rhs_lnro(tmplnro)
     call calc_rhs_vel(tmpu,tmpv,tmpw)
     call calc_rhs_roE(tmproE)
     call calc_rhs_Y0(tmpY0)
     if(nb.ne.0) call calc_rhs_nscbc(tmplnro,tmpu,tmpv,tmpw,tmproE,tmpY0)     
     
     !$omp parallel do
     do i=1,npfb
        !! Final values of prim vars
        lnro(i) = lnro_W(i) + RKb(4)*tmplnro(i)
        u(i) = u_W(i) + RKb(4)*tmpu(i)
        v(i) = v_W(i) + RKb(4)*tmpv(i)
        w(i) = w_W(i) + RKb(4)*tmpw(i)        
#ifndef isoT
        roE(i) = roE_W(i) + RKb(4)*tmproE(i)
#endif
#ifdef ms 
        Y0(i) = Y0_W(i) + RKb(4)*tmpY0(i)
#endif        
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(u_W,v_W,w_W,lnro_W,roE_W,tmpu,tmpv,tmpw,tmplnro,tmproE,Y0_W,tmpY0)     

     !! Set the new time   
     time = time0 + dt
     
     !! Apply BCs and update halos
     call reapply_mirror_bcs
     call halo_exchanges_all
          
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
     call reapply_mirror_bcs
     call halo_exchanges_all

     return
  end subroutine step_rk3_4S_2R
!! ------------------------------------------------------------------------------------------------
  subroutine step_rk3_4S_2R_EE
     !! 3rd order 4step 2 register Runge Kutta, with embedded 2nd order
     !! scheme for error estimation.     
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is lnro_W,u_W,v_W. Register 2 is lnro,u,v.
     !! Register 3 is tmplnro,tmpu,tmpv (only used for RHS)
     !! Register 4 is e_acc_lnro,e_acc_u,e_acc_v - error accumulator
     integer(ikind) :: i,k
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: u_W,v_W,lnro_W,roE_W,Y0_W,tmpu,tmpv,tmplnro,tmproE,tmpY0,w_W,tmpw
     real(rkind),dimension(:),allocatable :: e_acc_lnro,e_acc_u,e_acc_v,e_acc_E,e_acc_Y0,e_acc_w
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb,RKbmbh
     
     !! Push the max error storage back one
     emax_nm1 = emax_n;emax_n=emax_np1
     
     !! Set RKa,RKb,RKbmbh with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)
     RKbmbh(:) = dt*rk3_4s_2r_bmbh(:)

     allocate(u_W(npfb),v_W(npfb),lnro_W(npfb),roE_W(npfb),Y0_W(npfb),w_W(npfb))
     allocate(tmpu(npfb),tmpv(npfb),tmplnro(npfb),tmproE(npfb),tmpY0(npfb),tmpw(npfb))
     allocate(e_acc_lnro(npfb),e_acc_u(npfb),e_acc_v(npfb),e_acc_E(npfb),e_acc_Y0(npfb),e_acc_w(npfb))
     e_acc_lnro=zero;e_acc_u=zero;e_acc_v=zero;e_acc_E=zero;e_acc_Y0=zero;e_acc_w=zero

     !! Store prim vars in register 1 (w-register)
     !$omp parallel do
     do i=1,npfb
        lnro_W(i)=lnro(i);u_W(i)=u(i);v_W(i)=v(i);w_W(i)=w(i);roE_W(i)=roE(i);Y0_W(i)=Y0(i)
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time

     do k=1,3
        !! Calculate the RHS
        call calc_rhs_lnro(tmplnro)
        call calc_rhs_vel(tmpu,tmpv,tmpw)
        call calc_rhs_roE(tmproE)
        call calc_rhs_Y0(tmpY0)
        if(nb.ne.0) call calc_rhs_nscbc(tmplnro,tmpu,tmpv,tmpw,tmproE,tmpY0)        

        !! Set the intermediate time
!        time = time0 + RK_c(k)*dt        

        !! Set w_i and new u,v
        !$omp parallel do
        do i=1,npfb
           !! Store next U in register 2
           lnro(i) = lnro_W(i) + RKa(k)*tmplnro(i)
           u(i) = u_W(i) + RKa(k)*tmpu(i)
           v(i) = v_W(i) + RKa(k)*tmpv(i)
           w(i) = w_W(i) + RKa(k)*tmpw(i)           
#ifndef isoT           
           roE(i) = roE_W(i) + RKa(k)*tmproE(i)
#endif
#ifdef ms           
           Y0(i) = Y0_W(i) + RKa(k)*tmpY0(i)
#endif

           !! Store next S in register 1
           lnro_W(i) = lnro_W(i) + RKb(k)*tmplnro(i)
           u_W(i) = u_W(i) + RKb(k)*tmpu(i)
           v_W(i) = v_W(i) + RKb(k)*tmpv(i) 
           w_W(i) = w_W(i) + RKb(k)*tmpw(i)            
#ifndef isoT
           roE_W(i) = roE_W(i) + RKb(k)*tmproE(i)
#endif
#ifdef ms           
           Y0_W(i) = Y0_W(i) + RKb(k)*tmpY0(i)
#endif
           
           !! Error accumulation
           e_acc_lnro(i) = e_acc_lnro(i) + RKbmbh(k)*tmplnro(i)       
           e_acc_u(i) = e_acc_u(i) + RKbmbh(k)*tmpu(i)
           e_acc_v(i) = e_acc_v(i) + RKbmbh(k)*tmpv(i)  
           e_acc_w(i) = e_acc_w(i) + RKbmbh(k)*tmpw(i)             
#ifndef isoT
           e_acc_E(i) = e_acc_E(i) + RKbmbh(k)*tmproE(i)                    
#endif
#ifdef ms           
           e_acc_Y0(i) = e_acc_Y0(i) + RKbmbh(k)*tmpY0(i)
#endif
        end do
        !$omp end parallel do
       
        !! Apply BCs and update halos
        call reapply_mirror_bcs
        call halo_exchanges_all
        
     end do
     
     !! Final substep: returns solution straight to lnro,u,v,E (register 2)
     !! and doesn't update S
     call calc_rhs_lnro(tmplnro)
     call calc_rhs_vel(tmpu,tmpv,tmpw)
     call calc_rhs_roE(tmproE)
     call calc_rhs_Y0(tmpY0)
     if(nb.ne.0) call calc_rhs_nscbc(tmplnro,tmpu,tmpv,tmpw,tmproE,tmpY0)     
     
     enrm_ro=zero;enrm_u=zero;enrm_v=zero;enrm_E=zero;enrm_Y0=zero;enrm_w=zero
     !$omp parallel do reduction(max:enrm_ro,enrm_u,enrm_v,enrm_E,enrm_Y0,enrm_w)
     do i=1,npfb
        !! Final values of prim vars
        lnro(i) = lnro_W(i) + RKb(4)*tmplnro(i)
        u(i) = u_W(i) + RKb(4)*tmpu(i)
        v(i) = v_W(i) + RKb(4)*tmpv(i)
        w(i) = w_W(i) + RKb(4)*tmpw(i)        
#ifndef isoT
        roE(i) = roE_W(i) + RKb(4)*tmproE(i)
#endif
#ifdef ms 
        Y0(i) = Y0_W(i) + RKb(4)*tmpY0(i)
#endif        
        
        !! Final error accumulators
        e_acc_lnro(i) = e_acc_lnro(i) + RKbmbh(4)*tmplnro(i)       
        e_acc_u(i) = e_acc_u(i) + RKbmbh(4)*tmpu(i)
        e_acc_v(i) = e_acc_v(i) + RKbmbh(4)*tmpv(i) 
        e_acc_w(i) = e_acc_w(i) + RKbmbh(4)*tmpw(i)         
#ifndef isoT
        e_acc_E(i) = e_acc_E(i) + RKbmbh(4)*tmproE(i)   
#endif
#ifdef ms 
        e_acc_Y0(i) = e_acc_Y0(i) + RKbmbh(4)*tmpY0(i)
#endif
        
        !! Calculating L_infinity of error norms
        enrm_ro = abs(e_acc_lnro(i))!/(exp(lnro(i))+1.0d-9)      
        enrm_u = abs(e_acc_u(i))/(abs(u(i))+1.0d-9)  !! divide by zero mollification...
        enrm_v = abs(e_acc_v(i))/(abs(v(i))+1.0d-9)
        enrm_w = abs(e_acc_w(i))/(abs(w(i))+1.0d-9)        
#ifndef isoT
        enrm_E = abs(e_acc_E(i))/(abs(roE(i))+1.0d-9)
#endif
#ifdef ms 
        enrm_Y0 = abs(e_acc_Y0(i))/(abs(Y0(i))+1.0d-9)
#endif
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(u_W,v_W,lnro_W,roE_W,tmpu,tmpv,tmplnro,tmproE,Y0_W,tmpY0,w_W,tmpw) 
     deallocate(e_acc_lnro,e_acc_u,e_acc_v,e_acc_E,e_acc_Y0,e_acc_w)
     
     !! Finalise L_infinity error norms: find max and ensure it's >0     
     emax_np1 = max(max(max(enrm_ro,enrm_E),max(max(enrm_u,enrm_v),enrm_w)),1.0d-16)
     !! TBC: additional dependence on Y0 error norm.

     !! Set the new time   
     time = time0 + dt
     
     !! Apply BCs and update halos
     call reapply_mirror_bcs
     call halo_exchanges_all
          
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
     call reapply_mirror_bcs
     call halo_exchanges_all

     return
  end subroutine step_rk3_4S_2R_EE  
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     use thermodynamics
     integer(ikind) :: i
     real(rkind) :: cmax,umag,smin,dt_local
    
     !! Find maximum velocity magnitude
     umax = 1.0d-8;cmax = 1.0d-8
     !$OMP PARALLEL DO PRIVATE(umag) REDUCTION(max:umax,cmax)
     do i=1,npfb
        umag = u(i)*u(i) + v(i)*v(i) + w(i)*w(i)
        umax = umag
#ifndef isoT        
        cmax = (roE(i)/exp(lnro(i))-half*umag)
#endif        
     end do
     !$OMP END PARALLEL DO
     umax = sqrt(umax)
#ifndef isoT
     cmax = sqrt(cmax*gammagasm1*gammagas)     
#else
     cmax = sqrt(csq)
#endif          
     
     !! Find smallest node spacing
     smin = minval(s(1:npfb))

     !! Set time step
     dt = min(0.3*smin*smin/visc0,1.0*smin/(cmax+umax))
   
#ifdef mp     
     !! Find global time-step
     dt_local = dt
     call MPI_ALLREDUCE(dt_local,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror) 
     if(iproc.eq.0) then
        write(192,*) time,dt
        flush(192)
     end if
#else
     write(192,*) time,dt
     flush(192)         
#endif     

     return
  end subroutine set_tstep
!! ------------------------------------------------------------------------------------------------  
  subroutine set_tstep_PID
     integer(ikind) :: i
     real(rkind) :: dtfactor
     real(rkind) :: facA,facB,facC
     real(rkind) :: kappa,alph,beta,gamm,eps
     real(rkind) :: tratio_min,tratio_max
     real(rkind) :: umag2,smin
     real(rkind) :: dtmax,c,dt_local
     
     !! Use CFL to set an upper limit to dt
     smin = minval(s(1:npfb))     
     dtmax = 1.0d10
     !$omp parallel do private(umag2,c) reduction(min:dtmax)
     do i=1,npfb
        umag2 = u(i)*u(i) + v(i)*v(i) + w(i)*w(i)
#ifndef isoT        
        c = sqrt(gammagas*gammagasm1*(roE(i)/exp(lnro(i))-half*umag2))
#else
        c = sqrt(csq)        
#endif    
        dtmax = smin/(c+sqrt(umag2))  !! using smallest h and biggest speeds (not local)
     end do
     !$omp end parallel do

     dtmax = 0.8d0*dtmax
!     dtmax = min(dtmax,0.1*smin*smin/visc0) !! Viscous constraint...
     
     !! PID parameters...
     kappa=0.9
     alph=0.7/2.0;beta=0.4/2.0;gamm=0.1/2.0
     !0.7,0.4,0.1
     !0.49/p,0.34/p,0.1/p !! SENGA2
     eps = 1.0d-3
   
     !! P, I and D factors..  Calculation done in log space...
     facA = alph*log(eps/emax_np1)
     facB = beta*log(emax_n/eps)
     facC = gamm*log(eps/emax_nm1)
         
     !! Combined factor
     dtfactor = kappa*exp(facA+facB+facC)

     !! Limiting change in dt
     tratio_min = 1.0d-2
     tratio_max = 1.01d0
     if(dtfactor.lt.tratio_min) dtfactor = tratio_min
     if(dtfactor.gt.tratio_max) dtfactor = tratio_max     
     
     !! Set time new step
     dt = dt*dtfactor
     
     !! Impose upper limit
     if(dt.gt.dtmax) dt = dtmax
     
#ifdef mp     
     !! Find global time-step
     dt_local = dt
     call MPI_ALLREDUCE(dt_local,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror) 
     if(iproc.eq.0) then
        write(192,*) time,dt
        flush(192)
     end if
#else
     write(192,*) time,dt
     flush(192)         
#endif 
 
     return
  end subroutine set_tstep_PID
!! ------------------------------------------------------------------------------------------------  
end module step
