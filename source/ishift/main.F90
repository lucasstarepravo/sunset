program main
  use kind_parameters
  use common_parameter
  use common_2d
  use inputoutput 
  use neighbours
  use omp_lib
  implicit none

  integer(ikind) :: n,i
    
  write(6,*) "Enter number of processors in X direction"
  read(5,*) nprocsX
  if(nprocsX.lt.1) then
     write(6,*) "A positive number of processors in X please. Stopping"
     stop
  end if
  write(6,*) "Enter number of processors in Y direction"
  read(5,*) nprocsY
  if(nprocsY.lt.1) then
     write(6,*) "A positive number of processors in Y please. Stopping"
     stop
  end if  
  
  xbcond = 0;ybcond=0
  write(6,*) "Enter X boundary condition: periodic(1),symmetric(2)"
  read(5,*) xbcond
  write(6,*) "Enter Y boundary condition: periodic(1),symmetric(2)"
  read(5,*) ybcond
   
  write(6,*) "User specified processor grid:",nprocsX,nprocsY
  write(6,*) "Total # of processors:",nprocsX*nprocsY
  nprocs = nprocsX*nprocsY

  !! Initial conditions
  call initial_setup  

  !! Main bulk
  call setup_domain
  
  !! Re-arranging for output
  call remove_fd_nodes
  call rearrange_nodes
  
  !! Save nodes (to IPART in this directory
  call output_newnodes

  stop
end program main
!! ------------------------------------------------------------------------------------------------
subroutine initial_setup  
  use kind_parameters
  use common_parameter
  use common_2d
  use omp_lib
  
  !! Particles per smoothing length and supportsize/h
  hovs = 2.5   !! 2.2 normally, 2.5 can be safer...
  ss = 2.0
  nplink = 2*4*ceiling(ss*hovs)**2  !! # square stencil w/ side length 2*dr*hovs*ss
  
end subroutine initial_setup
!! ------------------------------------------------------------------------------------------------

