module output
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to write main output files (field data), perform reduction 
  !! calculations to obtain and output global statistics, and write data to standard-out.
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
     
     allocate(maxphi(6+nspec),minphi(6+nspec))
     
    
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
           write(6,*) "npfb,np",npfb_global,np_global,"n_out real",time/dt_out
           write(6,*) "Max |u|,|v|:",max(maxphi(1),abs(minphi(1))),max(maxphi(2),abs(minphi(2)))
           write(6,*) "Max |w|:",max(maxphi(3),abs(minphi(3)))
           write(6,*) "max/min ro:",maxphi(4),minphi(4)
           write(6,*) "max/min ro*E:",maxphi(5),minphi(5)
           write(6,*) "max/min T:",maxphi(6),minphi(6)   
           write(6,*) "There are",n_threads_global,"threads spread over ",nprocs,"MPI tasks"
           write(6,*) "Wall clock run time=",t_run_global
           write(6,*) "run-time/itime=",t_per_dt_global,"Moving avg=",t_last_X_global/dble(scr_freq)

           !! Profiling
           write(6,*) "  "
           write(6,*) "Profiling:::"
           write(6,*) "MPI transfers             :",segment_time_global(1)/sum(segment_time_global(1:8))
           write(6,*) "Mirror boundaries         :",segment_time_global(2)/sum(segment_time_global(1:8))
           write(6,*) "Filtering                 :",segment_time_global(3)/sum(segment_time_global(1:8))
           write(6,*) "Gradients                 :",segment_time_global(4)/sum(segment_time_global(1:8))
           write(6,*) "Laplacians                :",segment_time_global(5)/sum(segment_time_global(1:8))
           write(6,*) "Chemistry                 :",segment_time_global(6)/sum(segment_time_global(1:8))
           write(6,*) "Boundary 2nd derivatives  :",segment_time_global(7)/sum(segment_time_global(1:8))
           write(6,*) "RHS building              :",segment_time_global(8)/sum(segment_time_global(1:8))           
           write(6,'(/,/,A)') "  "
                                         
        end if

        t_last_X = 0.0d0
#else
        !! Single processor
        write(6,*)"itime,time,dt=", itime,time,dt
        write(6,*) "np,npfb",np,npfb,"n_out real",time/dt_out
        write(6,*) "Max |u|,|v|:",max(maxval(u(1:npfb)),abs(minval(u(1:npfb)))),max(maxval(v(1:npfb)),abs(minval(v(1:npfb))))
        write(6,*) "Max |w|:",max(maxval(w(1:npfb)),abs(minval(w(1:npfb))))
        write(6,*) "Max/min roE:",maxval(roE(1:npfb)),minval(roE(1:npfb))
        write(6,*) "Max/min T:",maxval(T(1:npfb)),minval(T(1:npfb))        
        write(6,*) "max/min ro:",exp(maxval(lnro(1:npfb))),exp(minval(lnro(1:npfb)))   
        write(6,*) "Number of threads=",n_threads,"Run time=",t_run
        write(6,*) "run-time/itime=",t_per_dt,"Moving avg=",t_last_X/dble(scr_freq)
        t_last_X = 0.0d0
        
        !! Profiling
        write(6,*) "  "
        write(6,*) "Profiling:::"
        write(6,*) "MPI transfers             :",segment_time_local(1)/sum(segment_time_local(1:8))
        write(6,*) "Mirror boundaries         :",segment_time_local(2)/sum(segment_time_local(1:8))
        write(6,*) "Filtering                 :",segment_time_local(3)/sum(segment_time_local(1:8))
        write(6,*) "Gradients                 :",segment_time_local(4)/sum(segment_time_local(1:8))
        write(6,*) "Laplacians                :",segment_time_local(5)/sum(segment_time_local(1:8))
        write(6,*) "Chemistry                 :",segment_time_local(6)/sum(segment_time_local(1:8))
        write(6,*) "Boundary 2nd derivatives  :",segment_time_local(7)/sum(segment_time_local(1:8))
        write(6,*) "RHS building              :",segment_time_local(8)/sum(segment_time_local(1:8))        
        write(6,'(/,/,/,A)') "  "                  
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
     real(rkind) :: tmpT,xn,yn,tmpro,tmpVort

     !! Only output from the first sheet
     if(iprocZ.eq.0)then     
!if(.true.)then     

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
        
        !! Number of species out
        nspec_out = 1
#ifdef output_composition
        nspec_out = nspec 
#endif                
        
        !! Write the main dump files
        open(unit = 20,file=fname)  
        write(20,*) np_out_local
        do i=1,np_out_local
           tmpro = exp(lnro(i))
#ifndef isoT
           tmpT = T(i)
#else
           tmpT = zero
#endif           
 
           !! Pass something to tmpVort (we use vorticity to output other things sometimes during
           !! debugging...)
           tmpVort = alpha_out(i)!vort(i)

#ifdef dim3
           write(20,*) rp(i,1),rp(i,2),rp(i,3),s(i),h(i),node_type(i),tmpro, &
                       u(i),v(i),w(i),tmpVort(i),tmpT,Yspec(i,1:nspec_out)        
#else
           write(20,*) rp(i,1),rp(i,2),s(i),h(i),node_type(i),tmpro,u(i),v(i),tmpVort, &
                       tmpT,Yspec(i,1:nspec_out)
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
  subroutine output_laminar_flame_structure(n_out)
     !! Output data for a laminar 1D flame profile.
     integer(ikind),intent(in) :: n_out !! Number of output file
#ifdef output_composition     
     integer(ikind) :: i,j,k
     real(rkind) :: x,y
     character(70) :: fname
   
     !! set the name of the file...
     !! first number is processor number, second is dump number (allowed up to 9999 processors)
     k = 10000        
#ifdef mp
     k = k + iproc
#endif                
     if( n_out .lt. 10 ) then 
        write(fname,'(A16,I5,A1,I1)') './data_out/flame',k,'_',n_out        
     else if( n_out .lt. 100 ) then 
        write(fname,'(A16,I5,A1,I2)') './data_out/flame',k,'_',n_out          
     else if( n_out .lt. 1000 ) then
        write(fname,'(A16,I5,A1,I3)') './data_out/flame',k,'_',n_out         
     else
        write(fname,'(A16,I5,A1,I4)') './data_out/flame',k,'_',n_out           
     end if   
        
     !! Write the main dump files
     open(unit = 20,file=fname)  
        
     do i=1,npfb
        x=rp(i,1);y=rp(i,2)
        if(abs(y).le.s(i)) then  !! For nodes within a node-spacing of y=0
           write(20,*) x,y,u(i),v(i),w(i),exp(lnro(i)),roE(i),T(i),p(i),Yspec(i,1:nspec)
!write(511,*) x,y,u(i),v(i),w(i),exp(lnro(i)),roE(i),T(i),p(i),Yspec(i,1:nspec)              
        end if
     end do
#endif
    return
  end subroutine output_laminar_flame_structure 
!! ------------------------------------------------------------------------------------------------  
end module output
