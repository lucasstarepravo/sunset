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
  public set_tstep,step_rk3_4S_2R,step_rk3_4S_2R_EE, &
         set_tstep_PID
  
  !! Error norms for RK3(2)4S[2R+]C scheme
  real(rkind) :: enrm_ro,enrm_u,enrm_v,enrm_E,enrm_w
  real(rkind),dimension(nspec) :: enrm_Yspec

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

     do k=1,3
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
        call reapply_mirror_bcs
        call halo_exchanges_all
        
        !! Velocity divergence
        call calc_divergence(u,v,w,divvel(1:npfb))
        call reapply_mirror_bcs_divvel_only
        call halo_exchange_divvel
        
     end do
     
     !! Final substep: returns solution straight to lnro,u,v,E (register 2)
     !! and doesn't update S
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
     call reapply_mirror_bcs
     call halo_exchanges_all
              
     !! Filter the solution 
     call filter_variables

     !! Apply BCs and update halos
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
     real(rkind) :: time0
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

     do k=1,3
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
        call reapply_mirror_bcs
        call halo_exchanges_all
        
        !! Velocity divergence
        call calc_divergence(u,v,w,divvel(1:npfb))
        call reapply_mirror_bcs_divvel_only
        call halo_exchange_divvel        
        
     end do
     
     !! Final substep: returns solution straight to lnro,u,v,E (register 2)
     !! and doesn't update S
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
        
        !! Calculating L_infinity of error norms
        enrm_ro = abs(e_acc_lnro(i))!/(exp(lnro(i))+1.0d-9)      
        enrm_u = abs(e_acc_u(i))/(abs(u(i))+1.0d-9)  !! divide by zero mollification...
        enrm_v = abs(e_acc_v(i))/(abs(v(i))+1.0d-9)
        enrm_w = abs(e_acc_w(i))/(abs(w(i))+1.0d-9)        
#ifndef isoT
        enrm_E = abs(e_acc_E(i))/(abs(roE(i))+1.0d-9)
#endif
#ifdef ms 
        do ispec=1,nspec
           enrm_Yspec(ispec) = abs(e_acc_Yspec(i,ispec))/(abs(Yspec(i,ispec))+1.0d-9)
        end do
#endif
     end do
     !$omp end parallel do  

     !! Deallocation
     deallocate(u_reg1,v_reg1,w_reg1,lnro_reg1,roE_reg1,Yspec_reg1)
     deallocate(rhs_u,rhs_v,rhs_w,rhs_lnro,rhs_roE,rhs_Yspec) 
     deallocate(e_acc_lnro,e_acc_u,e_acc_v,e_acc_E,e_acc_Yspec,e_acc_w)
     
     !! Finalise L_infinity error norms: find max and ensure it's >0     
     emax_np1 = max(max(max(enrm_ro,enrm_E),max(max(enrm_u,enrm_v),enrm_w)),1.0d-16)
     !! TBC: additional dependence on Yspec error norm.

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
     dt = min(0.3d0*smin*smin/visc0,one*smin/(cmax+umax))
   
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
