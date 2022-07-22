program sunset
  !! This is the main program of the sunset code.
  use kind_parameters
  use common_parameter
  use common_2d
  use setup
  use output 
  use neighbours
  use labf
  use fd
  use step
#ifdef mp  
  use mpi
#endif  
  implicit none

  integer(ikind) :: n_out,m_out

#ifdef mp  
  call MPI_INIT(ierror)
#endif  

  !! Initial conditions
  call initial_setup  
  call setup_domain

  !! Build the neighbour lists
  call find_neighbours

  !! Adapt the stencils by reducing h (only if not restarting)
#ifndef restart  
  call adapt_stencils
#endif  

  !! Calculate all the interparticle weights and any moments we might need
  call calc_labf_weights
  if(nb.ne.0) call calc_boundary_weights
  call calc_labf_sums
#ifdef dim3
  call calc_fd_weights
#endif    

  !! Calculate the filter coefficients
  call filter_coefficients   
  
  !! Create initial fields for primary variables
  call initial_solution

  !! Initialise the time-step (to something small to be safe...)
  call set_tstep;dt=0.0001*dt  

  !! Initialise time profiling and output counter...
  n_out = 0;ts_start = omp_get_wtime()
  m_out = 0
        
  !! MAIN TIME LOOP ---------------------------------------------------
  do while (time.le.time_end)
    
     !! Output, conditionally: at start, subsequently every dt_out
     if(itime.eq.0.or.time.gt.n_out*dt_out) then 
!     if(itime.eq.0.or.mod(itime,1).eq.0)then
        n_out = n_out + 1
        call output_layer(n_out)
     end if        
    
     !! Set the time step
     call set_tstep     
!     call set_tstep_PID  

     !! Perform one time step
!     call step_rk4
!     call step_rk3_2N
     call step_rk3_4S_2R
!     call step_rk3_4S_2R_EE     

     !! Calculate time-profiling and output to screen
     itime = itime + 1
     call output_to_screen

     !! Call routines to evaluate global statistics and adjust forcing terms if desired
     call statistics
     
  end do
  !! END MAIN TIME LOOP -----------------------------------------------
  
  !! Deallocate particle properties and neighbour lists
  call deallocate_weights
#ifdef mp    
!  call MPI_FINALIZE(ierror)
  call MPI_Abort(MPI_COMM_WORLD, n_out, ierror)       
#else
  stop
#endif  
end program sunset
!! ------------------------------------------------------------------------------------------------
subroutine deallocate_weights
  use kind_parameters
  use common_parameter
  use common_2d
  deallocate(rp,u,v,w,lnro,roE,Y0,s)
  deallocate(ij_count,ij_link)
  deallocate(irelation,vrelation)
  deallocate(node_type,internal_list)
  if(allocated(boundary_list)) deallocate(boundary_list)
  if(allocated(ij_w_grad)) then
     deallocate(ij_w_grad,ij_wb_grad2,ij_w_hyp)
     deallocate(ij_w_grad_sum,ij_wb_grad2_sum,ij_w_hyp_sum)     
  end if
  if(allocated(ij_link_fd)) then
     deallocate(ij_link_fd)
     deallocate(ij_fd_grad,ij_fd_grad2,ij_fd_hyp)
     deallocate(zlayer_index_global)
  end if
  return
end subroutine deallocate_weights
