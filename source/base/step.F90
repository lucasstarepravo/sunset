module step
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains time-stepping routines, and a routine to set the time-step.
  !! Time-stepping routines start with "step_" and perform one time step. They are
  !! called from the main loop, and they themselves call routines in the rhs module.

  use kind_parameters
  use common_parameter
  use common_vars
  use mirror_boundaries
  use characteristic_boundaries
  use rhs
  use mpi_transfers
#ifdef mp  
  use mpi
#endif  
  implicit none

  private
  public set_tstep,step_rk3_4S_2R,step_rk3_4S_2R_EE, &
         set_tstep_PID
  
  !! Error norms for RK3(2)4S[2R+]C scheme
  real(rkind) :: enrm_ro,enrm_u,enrm_v,enrm_E,enrm_w
  real(rkind),dimension(nspec_max) :: enrm_Yspec

contains
!! ------------------------------------------------------------------------------------------------  
  subroutine step_rk3_4S_2R
     !! 3rd order 4step 2 register Runge Kutta
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is lnro_reg1,u_reg1,v_reg1 etc
     !! Register 2 is lnro,u,v etc
     !! Register 3 is rhs_lnro, rhs_u, rhs_v etc
     use derivatives
     integer(ikind) :: i,k,ispec
     real(rkind) :: time0
     real(rkind),dimension(:),allocatable :: u_reg1,v_reg1,w_reg1,lnro_reg1,roE_reg1
     real(rkind),dimension(:,:),allocatable :: Yspec_reg1
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb
         
     !! Set RKa,RKb with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)

     allocate(u_reg1(npfb),v_reg1(npfb),lnro_reg1(npfb),roE_reg1(npfb),w_reg1(npfb))
     allocate(rhs_u(npfb),rhs_v(npfb),rhs_lnro(npfb),rhs_roE(npfb),rhs_w(npfb))
     allocate(Yspec_reg1(npfb,nspec),rhs_Yspec(npfb,nspec))

     !! Store primary variables in register 1 (w-register)
     !$omp parallel do private(ispec)
     do i=1,npfb
        lnro_reg1(i)=lnro(i);u_reg1(i)=u(i);v_reg1(i)=v(i);w_reg1(i)=w(i);roE_reg1(i)=roE(i)
        do ispec=1,nspec
           Yspec_reg1(i,ispec) = Yspec(i,ispec)
        end do
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time
     iRKstep=0

     do k=1,3
        iRKstep = iRKstep + 1
        !! Calculate the RHS
        call calc_all_rhs
     
        !! Set w_i and new u,v
        !$omp parallel do private(ispec)
        do i=1,npfb
           !! Store next U in register 2
           lnro(i) = lnro_reg1(i) + RKa(k)*rhs_lnro(i)
           u(i) = u_reg1(i) + RKa(k)*rhs_u(i)
           v(i) = v_reg1(i) + RKa(k)*rhs_v(i)
           w(i) = w_reg1(i) + RKa(k)*rhs_w(i)           
#ifndef isoT           
           roE(i) = roE_reg1(i) + RKa(k)*rhs_roE(i)
#endif
#ifdef ms           
           do ispec=1,nspec
              Yspec(i,ispec) = Yspec_reg1(i,ispec) + RKa(k)*rhs_Yspec(i,ispec)
           end do
#endif

           !! Store next S in register 1
           lnro_reg1(i) = lnro_reg1(i) + RKb(k)*rhs_lnro(i)
           u_reg1(i) = u_reg1(i) + RKb(k)*rhs_u(i)
           v_reg1(i) = v_reg1(i) + RKb(k)*rhs_v(i) 
           w_reg1(i) = w_reg1(i) + RKb(k)*rhs_w(i)            
#ifndef isoT
           roE_reg1(i) = roE_reg1(i) + RKb(k)*rhs_roE(i)
#endif
#ifdef ms           
           do ispec=1,nspec
              Yspec_reg1(i,ispec) = Yspec_reg1(i,ispec) + RKb(k)*rhs_Yspec(i,ispec)
           end do
#endif
        end do
        !$omp end parallel do
              
        !! Apply BCs and update halos
        call apply_time_dependent_bounds        
        call reapply_mirror_bcs
        call halo_exchanges_all
        
        !! Velocity divergence
        call calc_divergence(u,v,w,divvel(1:npfb))
        call reapply_mirror_bcs_divvel_only
        call halo_exchange_divvel
        
     end do
     
     !! Final substep: returns solution straight to lnro,u,v,E (register 2)
     !! and doesn't update S
     iRKstep = iRKstep + 1
     call calc_all_rhs  
     
     !$omp parallel do private(ispec)
     do i=1,npfb
        !! Final values of prim vars
        lnro(i) = lnro_reg1(i) + RKb(4)*rhs_lnro(i)
        u(i) = u_reg1(i) + RKb(4)*rhs_u(i)
        v(i) = v_reg1(i) + RKb(4)*rhs_v(i)
        w(i) = w_reg1(i) + RKb(4)*rhs_w(i)        
#ifndef isoT
        roE(i) = roE_reg1(i) + RKb(4)*rhs_roE(i)
#endif
#ifdef ms 
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec_reg1(i,ispec) + RKb(4)*rhs_Yspec(i,ispec)
        end do
#endif        
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(u_reg1,v_reg1,w_reg1,lnro_reg1,roE_reg1,Yspec_reg1)
     deallocate(rhs_u,rhs_v,rhs_w,rhs_lnro,rhs_roE,rhs_Yspec)     

     !! Set the new time   
     time = time0 + dt
     
     !! Apply BCs and update halos
     call apply_time_dependent_bounds     
     call reapply_mirror_bcs
     call halo_exchanges_all
              
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
     call apply_time_dependent_bounds
     call reapply_mirror_bcs
     call halo_exchanges_all

     !! Velocity divergence
     call calc_divergence(u,v,w,divvel(1:npfb))
     call reapply_mirror_bcs_divvel_only
     call halo_exchange_divvel     


     return
  end subroutine step_rk3_4S_2R
!! ------------------------------------------------------------------------------------------------
  subroutine step_rk3_4S_2R_EE
     use derivatives  
     !! 3rd order 4step 2 register Runge Kutta, with embedded 2nd order
     !! scheme for error estimation.     
     !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
     !! 'U' and 'S' in comments relate to p31 of Senga2 user guide
     !! Implemented over three registers, because speed is more
     !! important than memory at present.    
     
     !! Register 1 is lnro_reg1,u_reg1,v_reg1 etc
     !! Register 2 is lnro,u,v etc...
     !! Register 3 is rhs_lnro,rhs_u,rhs_v (only used for RHS)
     !! Register 4 is e_acc_lnro,e_acc_u,e_acc_v - error accumulator
     integer(ikind) :: i,k,ispec
     real(rkind) :: time0,emax_Y
     real(rkind),dimension(:),allocatable :: u_reg1,v_reg1,w_reg1,lnro_reg1,roE_reg1
     real(rkind),dimension(:),allocatable :: e_acc_lnro,e_acc_u,e_acc_v,e_acc_E,e_acc_w
     real(rkind),dimension(:,:),allocatable :: Yspec_reg1,e_acc_Yspec
     real(rkind),dimension(3) :: RKa
     real(rkind),dimension(4) :: RKb,RKbmbh
     
     !! Push the max error storage back one
     emax_nm1 = emax_n;emax_n=emax_np1
     
     !! Set RKa,RKb,RKbmbh with dt (avoids multiplying by dt on per-node basis)
     RKa(:) = dt*rk3_4s_2r_a(:)
     RKb(:) = dt*rk3_4s_2r_b(:)
     RKbmbh(:) = dt*rk3_4s_2r_bmbh(:)

     allocate(u_reg1(npfb),v_reg1(npfb),lnro_reg1(npfb),roE_reg1(npfb),w_reg1(npfb))
     allocate(rhs_u(npfb),rhs_v(npfb),rhs_lnro(npfb),rhs_roE(npfb),rhs_w(npfb))
     allocate(e_acc_lnro(npfb),e_acc_u(npfb),e_acc_v(npfb),e_acc_E(npfb),e_acc_w(npfb))
     allocate(Yspec_reg1(npfb,nspec),rhs_Yspec(npfb,nspec),e_acc_Yspec(npfb,nspec))
     e_acc_lnro=zero;e_acc_u=zero;e_acc_v=zero;e_acc_E=zero;e_acc_Yspec=zero;e_acc_w=zero

     !! Store prim vars in register 1 (w-register)
     !$omp parallel do private(ispec)
     do i=1,npfb
        lnro_reg1(i)=lnro(i);u_reg1(i)=u(i);v_reg1(i)=v(i);w_reg1(i)=w(i);roE_reg1(i)=roE(i)
        do ispec=1,nspec
           Yspec_reg1(i,ispec)=Yspec(i,ispec)
        end do
     end do
     !$omp end parallel do             

     !! Temporary storage of time
     time0=time
     iRKstep = 0

     do k=1,3
        iRKstep = iRKstep + 1
        !! Calculate the RHS        
        call calc_all_rhs      

        !! Set the intermediate time
!        time = time0 + RK_c(k)*dt        

        !! Set w_i and new u,v
        !$omp parallel do private(ispec)
        do i=1,npfb
           !! Store next U in register 2
           lnro(i) = lnro_reg1(i) + RKa(k)*rhs_lnro(i)
           u(i) = u_reg1(i) + RKa(k)*rhs_u(i)
           v(i) = v_reg1(i) + RKa(k)*rhs_v(i)
           w(i) = w_reg1(i) + RKa(k)*rhs_w(i)           
#ifndef isoT           
           roE(i) = roE_reg1(i) + RKa(k)*rhs_roE(i)
#endif
#ifdef ms           
           do ispec=1,nspec
              Yspec(i,ispec) = Yspec_reg1(i,ispec) + RKa(k)*rhs_Yspec(i,ispec)
           end do
#endif

           !! Store next S in register 1
           lnro_reg1(i) = lnro_reg1(i) + RKb(k)*rhs_lnro(i)
           u_reg1(i) = u_reg1(i) + RKb(k)*rhs_u(i)
           v_reg1(i) = v_reg1(i) + RKb(k)*rhs_v(i) 
           w_reg1(i) = w_reg1(i) + RKb(k)*rhs_w(i)            
#ifndef isoT
           roE_reg1(i) = roE_reg1(i) + RKb(k)*rhs_roE(i)
#endif
#ifdef ms           
           do ispec=1,nspec
              Yspec_reg1(i,ispec) = Yspec_reg1(i,ispec) + RKb(k)*rhs_Yspec(i,ispec)
           end do
#endif
           
           !! Error accumulation
           e_acc_lnro(i) = e_acc_lnro(i) + RKbmbh(k)*rhs_lnro(i)       
           e_acc_u(i) = e_acc_u(i) + RKbmbh(k)*rhs_u(i)
           e_acc_v(i) = e_acc_v(i) + RKbmbh(k)*rhs_v(i)  
           e_acc_w(i) = e_acc_w(i) + RKbmbh(k)*rhs_w(i)             
#ifndef isoT
           e_acc_E(i) = e_acc_E(i) + RKbmbh(k)*rhs_roE(i)                    
#endif
#ifdef ms           
           do ispec=1,nspec
              e_acc_Yspec(i,ispec) = e_acc_Yspec(i,ispec) + RKbmbh(k)*rhs_Yspec(i,ispec)
           end do
#endif
        end do
        !$omp end parallel do
       
        !! Apply BCs and update halos
        call apply_time_dependent_bounds        
        call reapply_mirror_bcs
        call halo_exchanges_all
        
        !! Velocity divergence
        call calc_divergence(u,v,w,divvel(1:npfb))
        call reapply_mirror_bcs_divvel_only
        call halo_exchange_divvel        
        
     end do
     
     !! Final substep: returns solution straight to lnro,u,v,E (register 2)
     !! and doesn't update S
     iRKstep = iRKstep + 1
     call calc_all_rhs    
     
     enrm_ro=zero;enrm_u=zero;enrm_v=zero;enrm_E=zero;enrm_Yspec=zero;enrm_w=zero
     !$omp parallel do private(ispec) reduction(max:enrm_ro,enrm_u,enrm_v,enrm_E,enrm_Yspec,enrm_w)
     do i=1,npfb
        !! Final values of prim vars
        lnro(i) = lnro_reg1(i) + RKb(4)*rhs_lnro(i)
        u(i) = u_reg1(i) + RKb(4)*rhs_u(i)
        v(i) = v_reg1(i) + RKb(4)*rhs_v(i)
        w(i) = w_reg1(i) + RKb(4)*rhs_w(i)        
#ifndef isoT
        roE(i) = roE_reg1(i) + RKb(4)*rhs_roE(i)
#endif
#ifdef ms 
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec_reg1(i,ispec) + RKb(4)*rhs_Yspec(i,ispec)
        end do
#endif        
        
        !! Final error accumulators
        e_acc_lnro(i) = e_acc_lnro(i) + RKbmbh(4)*rhs_lnro(i)       
        e_acc_u(i) = e_acc_u(i) + RKbmbh(4)*rhs_u(i)
        e_acc_v(i) = e_acc_v(i) + RKbmbh(4)*rhs_v(i) 
        e_acc_w(i) = e_acc_w(i) + RKbmbh(4)*rhs_w(i)         
#ifndef isoT
        e_acc_E(i) = e_acc_E(i) + RKbmbh(4)*rhs_roE(i)   
#endif
#ifdef ms 
        do ispec=1,nspec
           e_acc_Yspec(i,ispec) = e_acc_Yspec(i,ispec) + RKbmbh(4)*rhs_Yspec(i,ispec)
        end do
#endif
        
        !! Calculating L_infinity norm of errors. Variables "eX_norm" provide a value for divide-by-
        !! zero mollification, and this value is scaled according to the expected magnitude of the
        !! property. (e.g. eroE_norm is much bigger than eY_norm). They are set in common parameters.
        !! N.B. for initialising simulations with velocity discontinuities, it's helpful to relax
        !! the constraint on the velocity a bit by increasing eu_norm, ev_norm & ew_norm.
        enrm_ro = max(enrm_ro, &
                     abs(e_acc_lnro(i))/(lnro(i)+elnro_norm))      
        enrm_u = max(enrm_u, &
                     abs(e_acc_u(i))/(abs(u(i)) + eu_norm))
        enrm_v = max(enrm_v, &
                     abs(e_acc_v(i))/(abs(v(i)) + ev_norm))  
        enrm_w = max(enrm_w, &
                     abs(e_acc_w(i))/(abs(w(i)) + ew_norm))      
#ifndef isoT
        enrm_E = max(enrm_E, &
                     abs(e_acc_E(i))/(abs(roE(i)) + eroE_norm))
#endif
#ifdef ms 
        do ispec=1,nspec
           enrm_Yspec(ispec) = max(enrm_Yspec(ispec), &
                                   abs(e_acc_Yspec(i,ispec))/(abs(Yspec(i,ispec)) + eY_norm))
        end do
#endif
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(u_reg1,v_reg1,w_reg1,lnro_reg1,roE_reg1,Yspec_reg1)
     deallocate(rhs_u,rhs_v,rhs_w,rhs_lnro,rhs_roE,rhs_Yspec) 
     deallocate(e_acc_lnro,e_acc_u,e_acc_v,e_acc_E,e_acc_Yspec,e_acc_w)
     
     !! Finalise L_infinity error norms: find max and ensure it's >0     
#ifdef ms
     emax_Y = maxval(enrm_Yspec(1:nspec))     
#else
     emax_Y = zero
#endif         
     emax_np1 = max( &
                    max( &
                        max(enrm_ro,enrm_E), &
                        max( &
                            max(enrm_u,enrm_v),&
                            max(enrm_w,emax_Y) &
                            ) &
                        ), &
                    1.0d-16)

     !! Set the new time   
     time = time0 + dt
     
     !! Apply BCs and update halos
     call apply_time_dependent_bounds     
     call reapply_mirror_bcs
     call halo_exchanges_all
          
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
     call apply_time_dependent_bounds
     call reapply_mirror_bcs
     call halo_exchanges_all
     
     !! Velocity divergence
     call calc_divergence(u,v,w,divvel(1:npfb))
     call reapply_mirror_bcs_divvel_only
     call halo_exchange_divvel        

     return
  end subroutine step_rk3_4S_2R_EE  
!! ------------------------------------------------------------------------------------------------
  subroutine set_tstep
     use thermodynamics
     integer(ikind) :: i
     real(rkind) :: umag,dt_local
     real(rkind) :: dt_visc,dt_therm,dt_spec  
     real(rkind) :: c,uplusc,dt_cfl_local,dt_parabolic_local
         
     call evaluate_temperature_and_pressure
     call evaluate_transport_properties   
     
     !! Find minimum values for cfl, visc, thermal diff terms
     dt_cfl = 1.0d10;dt_visc = 1.0d10;dt_therm=1.0d10;dt_spec=1.0d10
     cmax = zero;umax = zero
!     !$omp parallel do private(c,uplusc) reduction(min:dt_cfl,dt_visc,dt_therm,dt_spec) &
!     !$omp reduction(max:cmax,umax)
     do i=1,npfb
        !! Sound speed 
#ifndef isoT        
        c = evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i))    
#else
        c = sqrt(csq)
#endif                
 
        !! Max velocity and sound speed
        umax = sqrt(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        cmax = max(c,cmax)
        
        !! Max speed of information propagation
        uplusc = umax + c
        
        !! Acoustic:: s/(u+c)
        !! Slightly reduce on outflows for stability
        if(node_type(i).eq.2) then
           dt_cfl = min(dt_cfl,0.8d0*s(i)/uplusc)
        else
           dt_cfl = min(dt_cfl,s(i)/uplusc)
        endif

        !! Viscous:: s*s*ro/visc
        dt_visc = min(dt_visc,s(i)*s(i)*exp(lnro(i))/visc(i))
        
#ifndef isoT        
        !! Thermal:: s*s*ro*cp/lambda_th
        dt_therm = min(dt_therm,s(i)*s(i)*exp(lnro(i))*cp(i)/lambda_th(i))
#endif      

#ifdef ms        
        !! Molecular diffusivity::  s*s*ro/Mdiff
        dt_spec = min(dt_spec,s(i)*s(i)*exp(lnro(i))/maxval(Mdiff(i,1:nspec)))
#endif   
         
     end do
!     !$omp end parallel do

     !! Scale by characteristic lengths and coefficients
     dt_cfl = one*dt_cfl*L_char
     dt_visc = 0.3d0*dt_visc*L_char*L_char
     dt_therm = 0.3d0*dt_therm*L_char*L_char
     dt_spec = 1.5d0*dt_spec*L_char*L_char
                           
     !! Find smallest node spacing
     smin = minval(s(1:npfb))*L_char

     !! Find most restrictive parabolic constraint
     dt_parabolic = min(dt_visc,min(dt_therm,dt_spec)) 
     
     !! Set dt if not reacting. If reacting, it will be set later by PID.
#ifndef react
     dt = min(dt_parabolic,dt_cfl)
#endif     
     
#ifdef mp     
     !! Global cfl-based time-step and parabolic parts based time-step
     dt_cfl_local = dt_cfl;dt_parabolic_local=dt_parabolic
     call MPI_ALLREDUCE(dt_cfl_local,dt_cfl,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(dt_parabolic_local,dt_parabolic,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)       
                 
     !! Output time-step (only if not reacting)
#ifndef react
     !! Global time step
     dt_local = dt
     call MPI_ALLREDUCE(dt_local,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     if(iproc.eq.0) then
        write(192,*) time,dt,one
        flush(192)
     end if
#endif     
#else
#ifndef react
     write(192,*) time,dt,one
     flush(192)         
#endif     
#endif     

     return
  end subroutine set_tstep
!! ------------------------------------------------------------------------------------------------  
  subroutine set_tstep_PID
     !! Adapt the time-step using the PID controller based on intergation errors.
     !! This routine is generally only called for reacting flows, and *presumes* that 
     !! the CFL-type time-step constraints have been calculated already.
     integer(ikind) :: i
     real(rkind) :: dtfactor,emax_local
     real(rkind) :: facA,facB,facC
     real(rkind) :: tratio_min,tratio_max
     real(rkind) :: umag2,umag
     real(rkind) :: c,dt_local,dt_max
        
                  
#ifdef mp     
     !! Parallel transfer to obtain the global maximum error          
     emax_local = emax_np1
     call MPI_ALLREDUCE(emax_local,emax_np1,1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror) 
#endif     
   
     !! P, I and D factors..  Calculation done in log space...
     facA = pid_a*log(pid_tol/emax_np1)
     facB = pid_b*log(emax_n/pid_tol)
     facC = pid_c*log(pid_tol/emax_nm1)
         
     !! Combined factor
     dtfactor = pid_k*exp(facA+facB+facC)
     
     !! Suppress big changes in time step (especially increases). N.B. this significantly reduces
     !! the stiffness of the PID system, and seems faster so far.
     dtfactor = one + one*atan((dtfactor-one)/one)
       
     !! Set time new step
     dt = dt*dtfactor
     
     !! Impose upper limit based on CFL-type constraints
     dt_max = min(dt_cfl,dt_parabolic)
     if(dt.gt.dt_max) dt = dt_max

#ifdef mp     
     !! Find global time-step
     dt_local = dt
     call MPI_ALLREDUCE(dt_local,dt,1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror) 
     if(iproc.eq.0) then
        write(192,*) time,dt,dt/dt_cfl
        flush(192)
     end if
#else
     write(192,*) time,dt,dt/dt_cfl
     flush(192)         
#endif 
 

     return
  end subroutine set_tstep_PID
!! ------------------------------------------------------------------------------------------------  
end module step
