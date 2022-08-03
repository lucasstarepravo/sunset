module output
  !! This module contains routines to read in node distributions, modify them with shifting 
  !! pre-process to obtain boundary normal vectors, and write field data out to file.
  use kind_parameters
  use common_parameter
  use common_2d
  use omp_lib
  use neighbours
#ifdef mp
  use mpi
  use mpi_transfers
#endif    
  implicit none
  real(rkind),dimension(:,:),allocatable :: gradu,gradv,gradw !! For vorticity and drag coefficient calculations

contains
!! ------------------------------------------------------------------------------------------------
  subroutine output_to_screen
     !! This routine writes information about the simulation to screen in a fairly easy to read
     !! format. For this to work (nicely) terminal window should be 24 lines tall.
     integer(ikind) :: scr_freq=100
     real(rkind),dimension(:),allocatable :: maxphi,minphi
     integer(ikind) :: n_threads_global
     real(rkind) :: t_per_dt_global,t_last_x_global,t_run_global
     real(rkind) :: cput,cput_global
     
     allocate(maxphi(5+nspec),minphi(5+nspec))
     
    
     ts_end=omp_get_wtime()
     t_run = t_run + ts_end - ts_start
     t_per_dt = t_run/dble(itime)
     t_last_X = t_last_X + ts_end - ts_start  
     !! Output cpu-time to file.
#ifdef mp
     cput = ts_end-ts_start
     call MPI_ALLREDUCE(cput,cput_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     if(iproc.eq.0) then 
        write(191,*) itime,cput_global
        flush(191)
     end if
#else
     write(191,*) itime,ts_end-ts_start
     flush(191)  
#endif  
  
     ! Some to screen
     if(mod(itime,scr_freq).eq.0)then 
  
#ifdef mp
        !! Multi-processor
        call reduce_for_screen_output(maxphi,minphi)
     
        call MPI_ALLREDUCE(n_threads,n_threads_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)          
        call MPI_ALLREDUCE(t_per_dt,t_per_dt_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
        call MPI_ALLREDUCE(t_last_X,t_last_X_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)         
        call MPI_ALLREDUCE(t_run,t_run_global,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)              
        t_last_X_global=t_last_x_global/dble(nprocs)
        t_per_dt_global =t_per_dt_global/dble(nprocs)
        t_run_global = t_run_global/dble(nprocs)
     
        !! Profiling bits
        call MPI_ALLREDUCE(segment_time_local,segment_time_global,10,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)  

        if(iproc.eq.0) then
 
           write(6,*)"itime,time,dt=", itime,time,dt
           write(6,*) "npfb,np",npfb_global,np_global
           write(6,*) "Max |u|,|v|:",max(maxphi(1),abs(minphi(1))),max(maxphi(2),abs(minphi(2)))
           write(6,*) "Max |w|:",max(maxphi(3),abs(minphi(3)))
           write(6,*) "max/min ro:",maxphi(4),minphi(4)
           write(6,*) "max/min ro*E:",maxphi(5),minphi(5)
           write(6,*) "max/min Y0:",maxphi(6),minphi(6)
    
           write(6,*) "There are",n_threads_global,"threads spread over ",nprocs,"MPI tasks"
           write(6,*) "Wall clock run time=",t_run_global
           write(6,*) "run-time/dt=",t_per_dt_global,"Moving avg=",t_last_X_global/dble(scr_freq)

           !! Profiling
           write(6,*) "  "
           write(6,*) "Profiling:::"
           write(6,*) "MPI transfers    :",segment_time_global(1)/sum(segment_time_global(1:8))
           write(6,*) "Mirror boundaries:",segment_time_global(2)/sum(segment_time_global(1:8))
           write(6,*) "Filtering        :",segment_time_global(3)/sum(segment_time_global(1:8))
           write(6,*) "RHS boundaries   :",segment_time_global(4)/sum(segment_time_global(1:8))
           write(6,*) "RHS energy       :",segment_time_global(5)/sum(segment_time_global(1:8))
           write(6,*) "RHS velocity     :",segment_time_global(6)/sum(segment_time_global(1:8))
           write(6,*) "RHS Density      :",segment_time_global(7)/sum(segment_time_global(1:8))
           write(6,*) "RHS Species      :",segment_time_global(8)/sum(segment_time_global(1:8))           
           write(6,'(/,/,A)') "  "
                    
           
           
!        write(6,'(/,/,/,/,/,/,/,/,/,/,/,/,/,A)') "  "
        end if

        t_last_X = 0.0d0
#else
        !! Single processor
        write(6,*)"itime,time,dt=", itime,time,dt
        write(6,*) "np,npfb",np,npfb
        write(6,*) "Max |u|,|v|:",max(maxval(u(1:npfb)),abs(minval(u(1:npfb)))),max(maxval(v(1:npfb)),abs(minval(v(1:npfb))))
        write(6,*) "Max |w|:",max(maxval(w(1:npfb)),abs(minval(w(1:npfb))))
        write(6,*) "max/min ro:",exp(maxval(lnro(1:npfb))),exp(minval(lnro(1:npfb)))   
        write(6,*) "Number of threads=",n_threads,"Run time=",t_run
        write(6,*) "run-time/dt=",t_per_dt,"Moving avg=",t_last_X/dble(scr_freq)
        t_last_X = 0.0d0
        
        !! Profiling
        write(6,*) "  "
        write(6,*) "Profiling:::"
        write(6,*) "MPI transfers    :",segment_time_local(1)/sum(segment_time_local(1:8))
        write(6,*) "Mirror boundaries:",segment_time_local(2)/sum(segment_time_local(1:8))
        write(6,*) "Filtering        :",segment_time_local(3)/sum(segment_time_local(1:8))
        write(6,*) "RHS boundaries   :",segment_time_local(4)/sum(segment_time_local(1:8))
        write(6,*) "RHS energy       :",segment_time_local(5)/sum(segment_time_local(1:8))
        write(6,*) "RHS velocity     :",segment_time_local(6)/sum(segment_time_local(1:8))
        write(6,*) "RHS Density      :",segment_time_local(7)/sum(segment_time_local(1:8))
        write(6,*) "RHS Species      :",segment_time_local(8)/sum(segment_time_local(1:8))        
        write(6,'(/,/,/,/,/,A)') "  "                  
#endif          
     
     end if
     
   
     ts_start=omp_get_wtime()
     return
  end subroutine output_to_screen    
!! ------------------------------------------------------------------------------------------------
  subroutine output_layer(n_out)
     !! Little subroutine to write out field variables. 
     !! Can be converted to vtk and read into paraview.
     use derivatives
     integer(ikind),intent(in) :: n_out
     integer(ikind) :: i,j,k,np_out_local,dimsout,nspec_out
     character(70) :: fname
     real(rkind),dimension(:),allocatable :: vort,testout
     real(rkind),dimension(:,:),allocatable :: testoutv
     real(rkind) :: tmpT,xn,yn,tmpro

     !! Only output from the first sheet
!     if(iprocZ.eq.0)then     
if(.true.)then     

        !! Calculate the vorticity 
        allocate(gradu(npfb,dims),gradv(npfb,dims),gradw(npfb,dims));gradw=zero
        allocate(vort(npfb))

        call calc_gradient(u,gradu)
        call calc_gradient(v,gradv)
#ifdef dim3
        call calc_gradient(w,gradw)
#endif     
        !$omp parallel do 
        do i=1,npfb
#ifdef dim3
           vort(i) = (gradw(i,2) - gradv(i,3))**two &
                   + (gradu(i,3) - gradw(i,1))**two &
                   + (gradv(i,1) - gradu(i,2))**two  !! Vorticity magnitude if 3D
           vort(i) = sqrt(vort(i))
#else     
           vort(i) = gradv(i,1) - gradu(i,2)
#endif        
        end do
        !$omp end parallel do

!! TEMPORARY OUTPUT DENSITY GRADIENT AS VORT.
!     call calc_gradient(lnro,gradu)
!     !$omp parallel do private(tmpT)
!     do i=1,npfb
!        tmpT = gradu(i,1)**two + gradu(i,2)**two
!        vort(i) = exp(lnro(i))*tmpT
!     end do
!     !$omp end parallel do
!     gradu=zero
!     gradv=zero
!! END TEMPORARY     
     
    
        if(nb.ne.0)then
           do j=1,nb  !! On WALL boundaries, vorticity requires some rotation...
              i=boundary_list(j)
              if(node_type(i).eq.0)then
                 xn = rnorm(i,1);yn = rnorm(i,2)
#ifdef dim3           
                 vort(i) = (gradw(i,2) - (xn*gradv(i,3)-yn*gradu(i,3)))**two &
                         + ((xn*gradu(i,3)+yn*gradv(i,3)) - gradw(i,1))**two &
                         + (xn*gradv(i,1)-yn*gradu(i,1)-xn*gradu(i,2)-yn*gradv(i,2))**two
                 vort(i) = sqrt(vort(i))                      
#else
                 vort(i) = xn*(gradv(i,1)-gradu(i,2)) - yn*(gradu(i,1) + gradv(i,2))
#endif
              end if
           end do
        end if     
        deallocate(gradu,gradv,gradw)     
          
        !! set the name of the file...
        !! first number is processor number, second is dump number (allowed up to 9999 processors)
#ifdef mp
        k=10000+iproc
#else
        k=10000
#endif     
        if( n_out .lt. 10 ) then 
           write(fname,'(A17,I5,A1,I1)') './data_out/layer_',k,'_',n_out        
        else if( n_out .lt. 100 ) then 
           write(fname,'(A17,I5,A1,I2)') './data_out/layer_',k,'_',n_out        
        else if( n_out .lt. 1000 ) then
           write(fname,'(A17,I5,A1,I3)') './data_out/layer_',k,'_',n_out        
        else
           write(fname,'(A17,I5,A1,I4)') './data_out/layer_',k,'_',n_out        
        end if 
     
        !! Local number nodes output
        np_out_local = npfb_layer 
        
        !! N species output
        nspec_out = nspec


        !! Write the main dump files
        open(unit = 20,file=fname)  
        write(20,*) np_out_local
        do i=1,np_out_local
           tmpro = exp(lnro(i))
#ifdef isoT
           tmpT = csq*(tmpro-one) !! if isoT, output |u| in E and pressure in tmpT. N.B. p_out=p-csq*rho0
           roE(i) = tmpro*sqrt(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))             
#else
           tmpT=(roE(i)/tmpro-0.5*(u(i)*u(i)+v(i)*v(i) + w(i)*w(i)))*gammagasm1/Rs0 !! The temperature  
#endif        
#ifdef dim3
           write(20,*) rp(i,1),rp(i,2),rp(i,3),s(i),node_type(i),tmpro, &
                       u(i),v(i),w(i),vort(i),roE(i)/tmpro,tmpT,Yspec(i,1)
!if(i.le.npfb) then        
!           write(20,*) rp(i,1)+0.5*iprocX,rp(i,2)+2.0*iprocY,rp(i,3)+0.5*iprocZ,s(i),node_type(i) &
!                       ,tmpro,u(i),v(i),w(i),dble(ij_count(i)) &
!                       ,roE(i)/tmpro,zero,divvel(i)  !! Diagnostic/debugging output
!else if(i.le.np_nohalo) then        
!           write(20,*) rp(i,1)+0.5*iprocX,rp(i,2)+2.0*iprocY,rp(i,3)+0.5*iprocZ,s(i),node_type(i) &
!                       ,tmpro,u(i),v(i),w(i),zero &
!                       ,roE(i)/tmpro,one,divvel(i)  !! Diagnostic/debugging output
!else 
!           write(20,*) rp(i,1)+0.5*iprocX,rp(i,2)+2.0*iprocY,rp(i,3)+0.5*iprocZ,s(i),node_type(i) &
!                       ,tmpro,u(i),v(i),w(i),zero &
!                       ,roE(i)/tmpro,two,divvel(i)  !! Diagnostic/debugging output
!end if             

#else
           tmpT=divvel(i)
           write(20,*) rp(i,1),rp(i,2),s(i),node_type(i),tmpro,u(i),v(i),vort(i), &
                       roE(i)/tmpro,tmpT,Yspec(i,1)
!if(i.le.npfb) then        
!           write(20,*) rp(i,1)+0.4*iprocX,rp(i,2)+4.4*iprocY,h(i),node_type(i),tmpro,u(i),v(i),dble(ij_count(i)) &
!                  ,roE(i)/tmpro,zero,Yspec(i,1)  !! Diagnostic/debugging output
!else if(i.le.np_nohalo) then        
!           write(20,*) rp(i,1)+0.4*iprocX,rp(i,2)+4.4*iprocY,h(i),node_type(i),tmpro,u(i),v(i),zero &
!                  ,roE(i)/tmpro,one,Yspec(i,1)  !! Diagnostic/debugging output
!else 
!           write(20,*) rp(i,1)+0.4*iprocX,rp(i,2)+4.4*iprocY,h(i),node_type(i),tmpro,u(i),v(i),zero &
!                  ,roE(i)/tmpro,two,Yspec(i,1)  !! Diagnostic/debugging output               
!end if             


#endif
        end do

        flush(20)
        close(20)
        if(allocated(vort)) deallocate(vort)     
     end if

     !! Write the time,dump number and # nodes to file
#ifdef dim3
     dimsout = 3
#else
     dimsout = 2
#endif    
#ifdef mp  
!     call MPI_ALLREDUCE(np_out_local,np_global,1,MPI_INT,MPI_SUM,MPI_COMM_WORLD,ierror)     
     if(iproc.eq.0) then 
        write(21,*) time,dimsout,npfb_global,n_out,nprocsX*nprocsY*nprocsZ,nspec_out
        flush(21)
     end if
#else
     write(21,*) time,dimsout,npfb,n_out,1,nspec_out
     flush(21)
#endif     


     return
  end subroutine output_layer
!! ------------------------------------------------------------------------------------------------
  subroutine statistics
     !! This is a temporary subroutine to calculate the total mass and energy in the domain
     integer(ikind) :: i
     real(rkind) :: tot_mass,tot_vol,tmpro,dVi,tot_mass_tmp,tot_vol_tmp,tot_roE,tot_roE_tmp
     
     !! Evaluate the mean velocity and adjust pressure gradient if required
     call velocity_control
     
     !! Check conservation of mass and energy
     call mass_and_energy_check
     
     !! Check mean internal energy and use PID controller and heat-sink if required
     call int_energy_control

     !! Error evaluation for Taylor Green vortices?
     call error_TG

     !! Calculate the lift and drag on all solid obstacles
!     call liftdrag

     !! Calculate how well balanced the MPI decomposition is
!     call check_load_balance

     !! Check the conservation of the species equations
     call species_check

     return
  end subroutine statistics  
!! ------------------------------------------------------------------------------------------------   
  subroutine check_load_balance  
     real(rkind),dimension(:),allocatable :: prof_tmp,prof_tmp_local,load_n_n
     integer(ikind),dimension(:),allocatable :: sum_n_n,sum_n_n_local
     integer(ikind) :: sum_sum_n_n
#ifdef mp     
     !! Calculate relative amounts of work being done by each processor
     allocate(prof_tmp_local(nprocs),prof_tmp(nprocs));prof_tmp_local=zero
     prof_tmp_local(iproc+1) = sum(segment_time_local(5:7))
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
     use derivatives
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
           sigma(1,1) = visc0*(fourthirds*gradu0(1)-twothirds*gradv0(2)) - csq*(exp(lnro(i))-one)
           sigma(1,2) = visc0*(gradu0(2)+gradv0(1))
           sigma(2,1) = visc0*(gradu0(2)+gradv0(1))           
           sigma(2,2) = visc0*(fourthirds*gradv0(2)-twothirds*gradu0(1)) - csq*(exp(lnro(i))-one)
          
           Fn(:) = matmul(sigma,rnorm(i,:))   !! Force on surface (sigma.n)
                     
           force(:) = force(:) + Fn(:)*s(i)  !! Integrate over surface... s(i) is node spacing
        end if
     end do
     !$omp end parallel do
     deallocate(gradu,gradv,gradw)

#ifdef mp
     force_tmp = force
     call MPI_ALLREDUCE(force_tmp,force,dims,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)     
     if(iproc.eq.0)then
        write(194,*) time,force
        flush(194)
     end if        
#else
     write(194,*) time,force
     flush(194)
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
        tmpro = exp(lnro(i))
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
        write(195,*) time,tot_vel,tot_u,driving_force(1)
        flush(195)
     end if
#else
     write(195,*) time,tot_vel,tot_u,driving_force(1)    
     flush(195)
#endif       
     
     
     return
  end subroutine velocity_control   
!! ------------------------------------------------------------------------------------------------
  subroutine mass_and_energy_check
     !! This subroutine calculates the mean density and total energy in the domain.
     integer(ikind) :: i
     real(rkind) :: tot_mass,tot_vol,tmpro,dVi,tot_mass_tmp,tot_vol_tmp,tot_roE,tot_roE_tmp
    
   
     tot_mass = zero;tot_vol = zero;tot_roE = zero
     !$omp parallel do private(tmpro,dVi) reduction(+:tot_mass,tot_vol,tot_roE)
     do i=1,npfb
        tmpro = exp(lnro(i))
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
     
#ifdef mp
     tot_mass_tmp = tot_mass;tot_vol_tmp = tot_vol;tot_roE_tmp = tot_roE
     call MPI_ALLREDUCE(tot_mass_tmp,tot_mass,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_vol_tmp,tot_vol,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_roE_tmp,tot_roE,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
#endif     
     
     !! Relative difference in mass from I.C.s
     tot_mass = (tot_mass-rho_char*tot_vol)/(rho_char*tot_vol)  
         
#ifdef mp
     if(iproc.eq.0)then
        write(193,*) time,tot_mass
        flush(193)
        write(197,*) time,tot_roE
        flush(197)
     end if
#else
     write(193,*) time,tot_mass
     flush(193)
     write(197,*) time,tot_roE
     flush(197)
#endif

     return
  end subroutine mass_and_energy_check 
!! ------------------------------------------------------------------------------------------------   
  subroutine species_check
     !! This subroutine calculates total quantity of each species in the domain
     integer(ikind) :: i,ispec
     real(rkind) :: tot_vol,tmpro,dVi,tot_vol_tmp,tmpY
     real(rkind),dimension(:),allocatable :: tot_Yspec,tot_Yspec_tmp
#ifdef ms     
     
     !! Allocate and zero accumulators
     allocate(tot_Yspec(nspec),tot_Yspec_tmp(nspec))   
     tot_Yspec = zero;tot_vol = zero
     
     !$omp parallel do private(tmpro,dVi,ispec,tmpY) reduction(+:tot_Yspec,tot_vol)
     do i=1,npfb
        tmpro = exp(lnro(i))
        dVi = s(i)*s(i) !! assume square nodes for now...
#ifdef dim3
        dVi = dVi*dz
#endif        

        do ispec=1,nspec
           tmpY = Yspec(i,ispec)
           tot_Yspec(ispec) = tot_Yspec(ispec) + dVi*tmpY*tmpro
        end do
!        tot_vol = tot_vol + dVi
     end do
     !$omp end parallel do
     
#ifdef mp
     tot_Yspec_tmp = tot_Yspec
     tot_vol_tmp = tot_vol;
     call MPI_ALLREDUCE(tot_Yspec_tmp,tot_Yspec,nspec,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
!    call MPI_ALLREDUCE(tot_vol_tmp,tot_vol,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     if(iproc.eq.0)then
        write(199,*) time,tot_Yspec(:)
        flush(199)
     end if
#else
     write(199,*) time,tot_Yspec(:)
     flush(199)
#endif

#endif

     return
  end subroutine species_check    
!! ------------------------------------------------------------------------------------------------
  subroutine int_energy_control
     !! Output the mean internal energy over the domain
     !! also set heat_sink_mag via a P.I.D. controller if required.
     integer(ikind) :: i
     real(rkind) :: tot_int_e,tot_vol,tmpro,dVi,tmpvel,tmp_int_e
     real(rkind) :: tot_int_e_tmp,tot_vol_tmp
       
     tot_int_e = zero
     tot_vol = zero
     !$omp parallel do private(tmpro,tmpvel,dVi,tmp_int_e) reduction(+:tot_int_e,tot_vol)
     do i=1,npfb
        tmpro = exp(lnro(i))
        dVi = s(i)*s(i) !! assume square nodes for now...
#ifdef dim3
        dVi = dVi*dz
#endif               
        tmpvel = (u(i)*u(i)+v(i)*v(i)+w(i)*w(i))
        tmp_int_e = roE(i)/tmpro - half*tmpvel  !! internal energy
        tot_int_e = tot_int_e + tmp_int_e*dVi
        tot_vol = tot_vol + dVi
     end do
     !$omp end parallel do

#ifdef mp
     tot_int_e_tmp = tot_int_e;tot_vol_tmp = tot_vol
     call MPI_ALLREDUCE(tot_int_e_tmp,tot_int_e,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(tot_vol_tmp,tot_vol,1,MPI_DOUBLE_PRECISION,MPI_SUM,MPI_COMM_WORLD,ierror)
#endif         

     !! Normalise over volume
     tot_int_e = (tot_int_e/tot_vol)

     !! Set the initial (and henceforwards target) value:
     if(itime.eq.1) mean_int_energy0 = tot_int_e
        
#ifdef mp
     if(iproc.eq.0)then
        write(198,*) time,tot_int_e
        flush(198)
     end if
#else
     write(198,*) time,tot_int_e
     flush(198)
#endif       
     
     
     return
  end subroutine int_energy_control 
!! ------------------------------------------------------------------------------------------------    
  subroutine poiseuille_l2norm
     !! Output the L2norm of the velocity compared with Poiseuille flow
     !! Umax=1, unit width domain...
     !! N.B. This routine needs updating for multi-processor simulations
     integer(ikind) :: i,j,jj
     real(rkind) :: y,uexact,local_error,sum_e,sum_exact,tot_vol,dVi
     real(rkind) :: N,X1,X2,y1
       
     sum_e = zero
     sum_exact = zero
     tot_vol = zero
     !$omp parallel do private(y,uexact,local_error,dVi,j,jj,N,X1,X2,y1) &
     !$omp reduction(+:sum_e,sum_exact,tot_vol)
     do i=1,npfb
        y = rp(i,2)
        uexact = (half-y)*(half+y)*4.0d0
        y1 = y+half


        jj = 30!floor(10.0d0/sqrt(time+0.1d0)); % vary the degree of expansion to avoid NaNs...
        do j=1,jj
           N=(2*j-1)*pi
           X1 = sin(N*y1)/N**3.0d0
           X2 = exp(-N*N*visc0*time)
           uexact = uexact - 32.0d0*X1*X2
        end do       
        
        dVi = h(i)*h(i) !! assume square nodes for now...
        local_error = u(i)-uexact
        sum_e = sum_e + local_error**two
        sum_exact = sum_exact + uexact**two
        tot_vol = tot_vol + dVi
     end do
     !$omp end parallel do
     sum_e = dsqrt(sum_e/sum_exact)

     write(196,*) time,sum_e
     flush(196)  
     return
  end subroutine poiseuille_l2norm  
!! ------------------------------------------------------------------------------------------------
  subroutine error_TG
    !! N.B. This routine needs updating for multiprocessor simulations
    use kind_parameters
    use common_parameter
    use common_2d
    implicit none
    integer(ikind) :: i
    real(rkind) :: u_exact,v_exact,x,y,expo
    real(rkind) :: U_ex,U_num,error_sum,L2error,Ex_sum,error_sum_local,ex_sum_local
  
    error_sum = zero;Ex_sum =zero
    !$omp parallel do private(x,y,u_exact,v_exact,U_ex,U_num) reduction(+:Ex_sum,error_sum)
    do i=1,npfb
       x=rp(i,1)
       y=rp(i,2)
       expo = exp(-8.0d0*pi*pi*time/Re)
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
    write(196,*) time,L2error
    flush(196) 
#endif    
  
  
    return
  end subroutine error_TG   
!! ------------------------------------------------------------------------------------------------ 
end module output
