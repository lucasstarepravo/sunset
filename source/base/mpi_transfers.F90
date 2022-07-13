module mpi_transfers
  !! This module contains routines to construct halos and do parallel data transfers
  !! ----------------------------------------------------------------------------------------------
  !! 2D DOMAIN DECOMPOSITION in x-y. Each processor has up to 8 neighbour processors
  !!
  !!     SEND ORDER    
  !!
  !!  ---------------------------   
  !!    6         |   | 2+npy+1       
  !!  --------------------------- 
  !!    5         |   | 2+npy+1     
  !!  ---------------------------  
  !!    4         | 1 | 2+npy+1      
  !!  ---------------------------   
  !!    3         |   | 2+npy+1       
  !!  ---------------------------   
  !!    2+npy     | 2 | 2+2*npy       
  !!  ---------------------------   
  !!    2+npy-1   |   | 2+2*npy-1      
  !!  ---------------------------   
  !!    2+npy-2   |   | 2+2*npy-2       
  !!  ---------------------------    
  !!
  !!   RECEIVE ORDER
  !!
  !!  ---------------------------   
  !!    2+2*npy-3 |   | 2+npy-3       
  !!  --------------------------- 
  !!    2+2*npy-2 |   | 2+npy-2     
  !!  ---------------------------  
  !!    2+2*npy-1 | 2 | 2+npy-1      
  !!  ---------------------------   
  !!    2+2*npy   |   | 2+npy       
  !!  ---------------------------   
  !!    2+npy+1   | 1 | 3       
  !!  ---------------------------   
  !!    2+npy+2   |   | 4      
  !!  ---------------------------   
  !!    2+npy+3   |   | 5       
  !!  ---------------------------      
  !!
  !!
  !! Transfers are SEND-RECEIVE for ODD processors, and RECEIVE-SEND for EVEN processors in EACH
  !! direction. Hence, if any periodic boundaries are required, need either 1 or EVEN number of 
  !! processors in periodic direction.
  !!
  !! UP-DOWN transfers are ALWAYS cyclic (in new framework...)
  !!  
  !! ----------------------------------------------------------------------------------------------
  use kind_parameters
  use common_parameter
  use common_2d
  use omp_lib
#ifdef mp
  use mpi
#endif 
  implicit none

contains
!! ------------------------------------------------------------------------------------------------  
  subroutine halo_exchanges_all  
     !! If using mpi, this calls routines to transfer all properties between halos. If not using
     !! mpi, it does nothing
     segment_tstart = omp_get_wtime()

#ifdef mp 
     !! u-velocity
     call halo_exchange(u)
    
     !! v-velocity
     call halo_exchange(v)

     !! w-velocity
#ifdef dim3
     call halo_exchange(w)
#endif     

     !! density logarithm
     call halo_exchange(lnro)

#ifndef isoT
     !! Energy
     call halo_exchange(roE)
#endif
#ifdef ms     
     !! Species
     call halo_exchange(Y0)
#endif     

#endif     

     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(1) = segment_time_local(1) + segment_tend - segment_tstart
     return
  end subroutine halo_exchanges_all
!! ------------------------------------------------------------------------------------------------
  subroutine halo_exchange_divvel  
     !! If using mpi, this calls routines to transfer divvel between halos. If not using
     !! mpi, it does nothing
     segment_tstart = omp_get_wtime()

#ifdef mp 
     !! Velocity divergence
     call halo_exchange(divvel)
#endif     


     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(1) = segment_time_local(1) + segment_tend - segment_tstart
     return
  end subroutine halo_exchange_divvel
#ifdef mp  
!! ------------------------------------------------------------------------------------------------
  subroutine processor_mapping
     integer(ikind) :: ii,i,j
     integer(ikind),dimension(:),allocatable :: iproc_thiscolumn
  
     !! Number of processors in X and Y decomposition - check match with schedule     
     if(nprocsX*nprocsY.ne.nprocs) then
        write(6,*) "ERROR: nprocs doesn't match the decomposition schedule from ishift. STOPPING"
        call MPI_Abort(MPI_COMM_WORLD, ii, ierror)
     end if
     
     !! Allocate space for the send and recieve lists
     allocate(iproc_S_LR(2*nprocsY),iproc_R_LR(2*nprocsY),iproc_S_UD(2),iproc_R_UD(2))
     iproc_S_LR=-1;iproc_R_LR=-1       !! Send-leftright,receive-leftright
     iproc_S_UD=-1;iproc_R_UD=-1       !! Send-updown, receive-updown
 
     !! Indices of this processor in X,Y grid 
     iprocX=iproc/nprocsY   
     iprocY=mod(iproc,nprocsY)
     
     !! Build a list of the indices of the processors in this column
     allocate(iproc_thiscolumn(nprocsY))
     iproc_thiscolumn(1) = iproc
     if(nprocsY.ne.1) then
        j=1
        do i=iprocY+1,nprocsY-1
           j = j+1
           iproc_thiscolumn(j) = iproc_thiscolumn(1) + j - 1
        end do
        do i=1,iprocY
           j = j+1
           iproc_thiscolumn(j) = iproc_thiscolumn(1) + j - 1 - nprocsY
        end do
     end if
     
     !! UP-DOWN SEND LIST
     if(nprocsY.ne.1) then
        !! Up, cyclic
        if(iprocY.ne.nprocsY-1) then
           iproc_S_UD(1) = iproc + 1
        else
           iproc_S_UD(1) = iproc + 1 - nprocsY
        end if
        !! Down, cyclic
        if(iprocY.ne.0)then
           iproc_S_UD(2) = iproc - 1
        else
           iproc_S_UD(2) = iproc - 1 + nprocsY
        end if
     end if
     
     !! LEFT-RIGHT RECEIVE LIST
     if(nprocsX.ne.1) then
        !! Left
        if(iprocX.ne.0) then
           iproc_S_LR(1:nprocsY) = iproc_thiscolumn(:) - nprocsY
        end if
        !! Right
        if(iprocX.ne.nprocsX-1) then
           iproc_S_LR(nprocsY+1:2*nprocsY) = iproc_thiscolumn(:) + nprocsY
        end if
        !! Left & periodic
        if(xbcond.eq.1.and.iprocX.eq.0) then
           iproc_S_LR(1:nprocsY) = iproc_thiscolumn(:) - nprocsY + nprocsX*nprocsY           
        end if              
        !! Right & periodic
        if(xbcond.eq.1.and.iprocX.eq.nprocsX-1)then
           iproc_S_LR(nprocsY+1:2*nprocsY) = iproc_thiscolumn(:) + nprocsY - nprocsX*nprocsY
        end if                         
     end if
     
     
     !! Rebuild thiscolumn list (reversed for RECEIVE)
     iproc_thiscolumn(1) = iproc
     if(nprocsY.ne.1) then
        j=1
        do i=1,iprocY
           j = j+1
           iproc_thiscolumn(j) = iproc_thiscolumn(1) - j + 1
        end do
        do i=iprocY+1,nprocsY-1
           j = j+1
           iproc_thiscolumn(j) = iproc_thiscolumn(1) - j + 1 + nprocsY
        end do
     end if   
 
   
     !! UP-DOWN RECEIVE LIST
     if(nprocsY.ne.1) then
        !! Down, cyclic
        if(iprocY.ne.0) then
           iproc_R_UD(1) = iproc - 1
        else
           iproc_R_UD(1) = iproc - 1 + nprocsY
        end if
        !! Up, cyclic
        if(iprocY.ne.nprocsY-1)then
           iproc_R_UD(2) = iproc + 1
        else
           iproc_R_UD(2) = iproc + 1 - nprocsY
        end if
     end if
     
     !! LEFT-RIGHT RECEIVE LIST
     if(nprocsX.ne.1) then
        !! Right
        if(iprocX.ne.nprocsX-1) then
           iproc_R_LR(1:nprocsY) = iproc_thiscolumn(:) + nprocsY
        end if
        !! Left
        if(iprocX.ne.0) then
           iproc_R_LR(nprocsY+1:2*nprocsY) = iproc_thiscolumn(:) - nprocsY
        end if
        !! right & periodic
        if(xbcond.eq.1.and.iprocX.eq.nprocsX-1) then
           iproc_R_LR(1:nprocsY) = iproc_thiscolumn(:) + nprocsY - nprocsX*nprocsY           
        end if              
        !! Left & periodic
        if(xbcond.eq.1.and.iprocX.eq.0)then
           iproc_R_LR(nprocsY+1:2*nprocsY) = iproc_thiscolumn(:) - nprocsY + nprocsX*nprocsY
        end if                                      
     end if   


     write(6,*) "new lists-S-LR",iproc,iproc_S_LR(:)
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)          
     write(6,*) "new lists-R-LR",iproc,iproc_R_LR(:)     
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)          
     write(6,*) "new lists-S-UD",iproc,iproc_S_UD(:)
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)          
     write(6,*) "new lists-R-UD",iproc,iproc_R_UD(:)     
     call MPI_BARRIER( MPI_COMM_WORLD, ierror)     
       
     
     return
  end subroutine processor_mapping
!! ------------------------------------------------------------------------------------------------
  subroutine build_halos
     integer(ikind) :: i,nstart,nend,maxhalo,suminhalo,is,ie,k
     integer(ikind) :: maxhalo_UD,maxhalo_LR
     real(rkind) :: x,y,ss_local,halo_fac
     integer(ikind),dimension(:,:),allocatable :: halo_lists_tmp
     integer(ikind),dimension(:,:),allocatable :: halo_lists_tmp_LR,halo_lists_tmp_UD
     !! Routine builds a list of nodes which are to be exported as halos for adjacent processors

     !! How much extra?
     halo_fac = 2.0d0
          
     !! Space for halos in temporary arrays
     allocate(halo_lists_tmp_LR(np,2*nprocsY),halo_lists_tmp_UD(np,2))
     
     !! Prepare space for halo sizes
     allocate(nhalo_LR(2*nprocsY),nhalo_UD(2));nhalo_LR=0;nhalo_UD=0
     allocate(inhalo_LR(2*nprocsY),inhalo_UD(2));inhalo_LR=0;inhalo_UD=0
     
     !! Build halo UP
     if(iproc_S_UD(1).ge.0) then
        do i=1,np
           y=rp(i,2);ss_local = halo_fac*ss*h(i)      
           if(abs(y-YU_thisproc).le.ss_local) then     
              nhalo_UD(1) = nhalo_UD(1) + 1
              halo_lists_tmp_UD(nhalo_UD(1),1) = i              
           end if
        end do
     end if

     !! Build halo DOWN     
     if(iproc_S_UD(2).ge.0)then            
        do i=1,np
           y=rp(i,2);ss_local = halo_fac*ss*h(i)
           if(abs(y-YD_thisproc).le.ss_local) then
              nhalo_UD(2) = nhalo_UD(2) + 1
              halo_lists_tmp_UD(nhalo_UD(2),2) = i              
           end if
        end do           
     end if     

     !! Build halos LEFT
     do k=1,nprocsY
        if(iproc_S_LR(k).ge.0) then           
           do i=1,np
              x=rp(i,1);ss_local = halo_fac*ss*h(i)      
              if(abs(x-XL_thisproc).le.ss_local) then     
                 nhalo_LR(k) = nhalo_LR(k) + 1
                 halo_lists_tmp_LR(nhalo_LR(k),k) = i
              end if
           end do
        end if
     end do

     !! Build halos RIGHT
     do k=nprocsY+1,2*nprocsY
        if(iproc_S_LR(k).ge.0) then           
           do i=1,np
              x=rp(i,1);ss_local = halo_fac*ss*h(i)      
              if(abs(x-XR_thisproc).le.ss_local) then     
                 nhalo_LR(k) = nhalo_LR(k) + 1
                 halo_lists_tmp_LR(nhalo_LR(k),k) = i
              end if
           end do
        end if
     end do

     !! Store halos in final array
     maxhalo_LR = maxval(nhalo_LR(1:2*nprocsY))
     maxhalo_UD = maxval(nhalo_UD(1:2))
     allocate(halo_lists_LR(maxhalo_LR,2*nprocsY),halo_lists_UD(maxhalo_UD,2))
     halo_lists_LR=0;halo_lists_UD=0

     do k=1,2
        if(iproc_S_UD(k).ne.-1)then
           halo_lists_UD(1:nhalo_UD(k),k) = halo_lists_tmp_UD(1:nhalo_UD(k),k)
        end if
     end do
     do k=1,2*nprocsY
        if(iproc_S_LR(k).ne.-1)then
           halo_lists_LR(1:nhalo_LR(k),k) = halo_lists_tmp_LR(1:nhalo_LR(k),k)
        end if
     end do         
     deallocate(halo_lists_tmp_UD,halo_lists_tmp_LR)

     !! Exchange halo sizes and node positions between processors   
     call send_halo_sizes_nodes
     
    
     suminhalo = sum(inhalo_LR(1:2*nprocsY)) + sum(inhalo_UD(1:2))         
     np = np_nohalo + suminhalo     
            
     return
  end subroutine build_halos  
!! ------------------------------------------------------------------------------------------------  
  subroutine halo_exchange(phi)    
     !! This routine does halo exchanges for phi= u,v,lnro,E,Y0
     real(rkind),dimension(:),intent(inout) :: phi     
     real(rkind),dimension(:),allocatable :: halo_phi
     integer(ikind) :: i,j,k,tag
     logical :: xodd,yodd

     xodd = .true.     
     if(mod(iprocX,2).eq.0) then
        xodd = .false.
     end if
     yodd = .true.     
     if(mod(iprocY,2).eq.0) then
        yodd = .false.
     end if
     
     
     !! Up and down first
     do k=1,2
        if(yodd) then !! ODD Y: SEND FIRST, RECEIVE SECOND
           if(iproc_S_UD(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_UD(k)))
              do i=1,nhalo_UD(k)
                 j=halo_lists_UD(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_UD(k),MPI_DOUBLE_PRECISION,iproc_S_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_UD(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_UD(k) + 100*k
              call MPI_RECV(phi(nrecstart_UDLR(k):nrecstart_UDLR(k)+inhalo_UD(k)-1),&
                            inhalo_UD(k),MPI_DOUBLE_PRECISION,iproc_R_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN Y: RECEIVE FIRST, SEND SECOND
           if(iproc_R_UD(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_UD(k) + 100*k
              call MPI_RECV(phi(nrecstart_UDLR(k):nrecstart_UDLR(k)+inhalo_UD(k)-1), &
                            inhalo_UD(k),MPI_DOUBLE_PRECISION,iproc_R_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_UD(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_UD(k)))
              do i=1,nhalo_UD(k)
                 j=halo_lists_UD(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_UD(k),MPI_DOUBLE_PRECISION,iproc_S_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do
     
     !! All left and right
     do k=1,2*nprocsY
        if(xodd) then !! ODD X: SEND FIRST, RECEIVE SECOND
           if(iproc_S_LR(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_LR(k)))
              do i=1,nhalo_LR(k)
                 j=halo_lists_LR(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_LR(k),MPI_DOUBLE_PRECISION,iproc_S_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_LR(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_LR(k) + 100*k
              call MPI_RECV(phi(nrecstart_UDLR(2+k):nrecstart_UDLR(2+k)+inhalo_LR(k)-1),&
                            inhalo_LR(k),MPI_DOUBLE_PRECISION,iproc_R_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN X: RECEIVE FIRST, SEND SECOND
           if(iproc_R_LR(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_LR(k) + 100*k
              call MPI_RECV(phi(nrecstart_UDLR(2+k):nrecstart_UDLR(2+k)+inhalo_LR(k)-1), &
                            inhalo_LR(k),MPI_DOUBLE_PRECISION,iproc_R_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_LR(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_LR(k)))
              do i=1,nhalo_LR(k)
                 j=halo_lists_LR(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_LR(k),MPI_DOUBLE_PRECISION,iproc_S_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do                                    
     
     return
  end subroutine halo_exchange
!! ------------------------------------------------------------------------------------------------  
  subroutine halo_exchange_int(phi) 
     !! This routine does halo exchanges for integers (only used at startup)   
     integer(ikind),dimension(:),intent(inout) :: phi     
     integer(ikind),dimension(:),allocatable :: halo_phi
     integer(ikind) :: i,j,k,tag
     logical :: xodd,yodd

     xodd = .true.     
     if(mod(iprocX,2).eq.0) then
        xodd = .false.
     end if
     yodd = .true.     
     if(mod(iprocY,2).eq.0) then
        yodd = .false.
     end if

     !! Up and down first
     do k=1,2
        if(yodd) then !! ODD Y: SEND FIRST, RECEIVE SECOND
           if(iproc_S_UD(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_UD(k)))
              do i=1,nhalo_UD(k)
                 j=halo_lists_UD(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_UD(k),MPI_INT,iproc_S_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_UD(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_UD(k) + 100*k
              call MPI_RECV(phi(nrecstart_UDLR(k):nrecstart_UDLR(k)+inhalo_UD(k)-1),&
                            inhalo_UD(k),MPI_INT,iproc_R_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN Y: RECEIVE FIRST, SEND SECOND
           if(iproc_R_UD(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_UD(k) + 100*k
              call MPI_RECV(phi(nrecstart_UDLR(k):nrecstart_UDLR(k)+inhalo_UD(k)-1), &
                            inhalo_UD(k),MPI_INT,iproc_R_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_UD(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_UD(k)))
              do i=1,nhalo_UD(k)
                 j=halo_lists_UD(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_UD(k),MPI_INT,iproc_S_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do
     
     !! All left and right
     do k=1,2*nprocsY
        if(xodd) then !! ODD X: SEND FIRST, RECEIVE SECOND
           if(iproc_S_LR(k).ge.0)then !! Send neighbour exists
              !! Put data in an array for sending
              allocate(halo_phi(nhalo_LR(k)))
              do i=1,nhalo_LR(k)
                 j=halo_lists_LR(i,k)
                 halo_phi(i) = phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_LR(k),MPI_INT,iproc_S_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
           if(iproc_R_LR(k).ge.0) then !! Receive neighbour exists
              tag = iproc_R_LR(k) + 100*k
              call MPI_RECV(phi(nrecstart_UDLR(2+k):nrecstart_UDLR(2+k)+inhalo_LR(k)-1),&
                            inhalo_LR(k),MPI_INT,iproc_R_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
        else !! EVEN X: RECEIVE FIRST, SEND SECOND
           if(iproc_R_LR(k).ge.0)then !! Receive neighbour exists
              !! Receive the data
              tag = iproc_R_LR(k) + 100*k
              call MPI_RECV(phi(nrecstart_UDLR(2+k):nrecstart_UDLR(2+k)+inhalo_LR(k)-1), &
                            inhalo_LR(k),MPI_INT,iproc_R_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
           end if
           if(iproc_S_LR(k).ge.0)then !! Send neighbour exists
              allocate(halo_phi(nhalo_LR(k)))
              do i=1,nhalo_LR(k)
                 j=halo_lists_LR(i,k)
                 halo_phi(i)=phi(j)
              end do
              
              !! Send the data
              tag = iproc + 100*k
              call MPI_SEND(halo_phi,nhalo_LR(k),MPI_INT,iproc_S_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_phi)
           end if
        end if
     end do     
                       
     return
  end subroutine halo_exchange_int
!! ------------------------------------------------------------------------------------------------  
  subroutine reduce_for_screen_output(maxphi,minphi)  
     !! Find the maximum and minimum values for output to screen
     real(rkind),dimension(:),intent(out) :: maxphi,minphi
     
     
     
     call MPI_ALLREDUCE(maxval(u(1:np)),maxphi(1),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(v(1:np)),maxphi(2),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(maxval(w(1:np)),maxphi(3),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(lnro(1:np)),maxphi(4),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(roE(1:np)),maxphi(5),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(maxval(Y0(1:np)),maxphi(6),1,MPI_DOUBLE_PRECISION,MPI_MAX,MPI_COMM_WORLD,ierror)                    
  
     call MPI_ALLREDUCE(minval(u(1:np)),minphi(1),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(v(1:np)),minphi(2),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)               
     call MPI_ALLREDUCE(minval(w(1:np)),minphi(3),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(lnro(1:np)),minphi(4),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(roE(1:np)),minphi(5),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)
     call MPI_ALLREDUCE(minval(Y0(1:np)),minphi(6),1,MPI_DOUBLE_PRECISION,MPI_MIN,MPI_COMM_WORLD,ierror)                    
  
     !! Density from density logarithm
     maxphi(4) = exp(maxphi(4));minphi(4)=exp(minphi(4))
  
     return
  end subroutine reduce_for_screen_output
!! ------------------------------------------------------------------------------------------------  
  subroutine refine_halos
     !! This routine uses the rough halos (based on distance from XL,XR,YU,YD _thisproc), and
     !! checks to see whether a node has a neighbour in that halo. If it does, it is marked for a
     !! refined halo. Reduces halo size to (nearly) minimum.
     !! Assumes that find_neighbours has already been called...
     integer(ikind),dimension(:),allocatable :: halo_essential_all
     integer(ikind),dimension(:),allocatable :: halo_essential_from,halo_essential,halo_list_tmp
     integer(ikind) i,j,k,jstart,maxhalo,suminhalo,nhalo_new,tag
     logical :: xodd,yodd

     !! Identify required halo nodes
     allocate(halo_essential_all(np));halo_essential_all = 0
     !$omp parallel do private(j,k)
     do i=1,npfb
        do k=1,ij_count(i)
           j=ij_link(k,i)
           if(j.gt.np_nohalo) halo_essential_all(j) = 1
        end do
     end do
     !$omp end parallel do

     xodd = .true.     
     if(mod(iprocX,2).eq.0) then
        xodd = .false.
     end if
     yodd = .true.     
     if(mod(iprocY,2).eq.0) then
        yodd = .false.
     end if          
     
     !! Up and down first
     do k=1,2
        if(yodd) then !! ODD Y: SEND FIRST, RECEIVE SECOND
           if(iproc_R_UD(k).ge.0)then !! This neighbour sent something
              !! Build halo_essential_from ready for transfer
              allocate(halo_essential_from(inhalo_UD(k)))
              halo_essential_from = halo_essential_all(nrecstart_UDLR(k):nrecstart_UDLR(k)+inhalo_UD(k)-1)

              !! Send the data back to originator processor
              tag = 100*k+iproc
              call MPI_SEND(halo_essential_from,inhalo_UD(k),MPI_INT,iproc_R_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_essential_from)
           end if
           if(iproc_S_UD(k).ge.0) then !! Receive neighbour exists
              allocate(halo_essential(nhalo_UD(k)))
              tag = 100*k+iproc_S_UD(k)
              call MPI_RECV(halo_essential,nhalo_UD(k),MPI_INT,iproc_S_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
                            
              !! Count the new halo size
              nhalo_new = 0
              do i=1,nhalo_UD(k)
                 if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
              end do
              write(6,*) iproc,"Reduced halo_UD",k," from",nhalo_UD(k),"to",nhalo_new              
              
              !! Build the new halo lists...
              allocate(halo_list_tmp(nhalo_new))
              j=0
              do i=1,nhalo_UD(k)
                 if(halo_essential(i).eq.1)then
                    j=j+1
                    halo_list_tmp(j) = halo_lists_UD(i,k)
                 end if
              end do
              
              !! Move new halo list to old halo list
              halo_lists_UD(:,k)=0;halo_lists_UD(1:nhalo_new,k) = halo_list_tmp;nhalo_UD(k) = nhalo_new
              deallocate(halo_list_tmp)
              deallocate(halo_essential)                                         
           end if
        else !! EVEN Y: RECEIVE FIRST, SEND SECOND
           if(iproc_S_UD(k).ge.0) then !! Receive neighbour exists
              allocate(halo_essential(nhalo_UD(k)))
              tag = 100*k+iproc_S_UD(k)
              call MPI_RECV(halo_essential,nhalo_UD(k),MPI_INT,iproc_S_UD(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)

              !! Count the new halo size
              nhalo_new = 0
              do i=1,nhalo_UD(k)
                 if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
              end do
              write(6,*) iproc,"Reduced halo_UD",k," from",nhalo_UD(k),"to",nhalo_new              
              
              !! Build the new halo lists...
              allocate(halo_list_tmp(nhalo_new))
              j=0
              do i=1,nhalo_UD(k)
                 if(halo_essential(i).eq.1)then
                    j=j+1
                    halo_list_tmp(j) = halo_lists_UD(i,k)
                 end if
              end do
              
              !! Move new halo list to old halo list
              halo_lists_UD(:,k)=0;halo_lists_UD(1:nhalo_new,k) = halo_list_tmp;nhalo_UD(k) = nhalo_new
              deallocate(halo_list_tmp)
              deallocate(halo_essential)                                         
           end if
           if(iproc_R_UD(k).ge.0)then !! This neighbour sent something
              !! Build halo_essential_from ready for transfer
              allocate(halo_essential_from(inhalo_UD(k)))
              halo_essential_from = halo_essential_all(nrecstart_UDLR(k):nrecstart_UDLR(k)+inhalo_UD(k)-1)

              !! Send the data back to originator processor
              tag = 100*k+iproc
              call MPI_SEND(halo_essential_from,inhalo_UD(k),MPI_INT,iproc_R_UD(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_essential_from)
           end if
        end if
     end do
     
     !! Then left and right
     do k=1,2*nprocsY
        if(xodd) then !! ODD X: SEND FIRST, RECEIVE SECOND
           if(iproc_R_LR(k).ge.0)then !! This neighbour sent something
              !! Build halo_essential_from ready for transfer
              allocate(halo_essential_from(inhalo_LR(k)))
              halo_essential_from = halo_essential_all(nrecstart_UDLR(k+2):nrecstart_UDLR(k+2)+inhalo_LR(k)-1)

              !! Send the data back to originator processor
              tag = 100*k+iproc
              call MPI_SEND(halo_essential_from,inhalo_LR(k),MPI_INT,iproc_R_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_essential_from)
           end if
           if(iproc_S_LR(k).ge.0) then !! Receive neighbour exists
              allocate(halo_essential(nhalo_LR(k)))
              tag = 100*k+iproc_S_LR(k)
              call MPI_RECV(halo_essential,nhalo_LR(k),MPI_INT,iproc_S_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)
                            
              !! Count the new halo size
              nhalo_new = 0
              do i=1,nhalo_LR(k)
                 if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
              end do
              write(6,*) iproc,"Reduced halo_LR",k," from",nhalo_LR(k),"to",nhalo_new              
              
              !! Build the new halo lists...
              allocate(halo_list_tmp(nhalo_new))
              j=0
              do i=1,nhalo_LR(k)
                 if(halo_essential(i).eq.1)then
                    j=j+1
                    halo_list_tmp(j) = halo_lists_LR(i,k)
                 end if
              end do
              
              !! Move new halo list to old halo list
              halo_lists_LR(:,k)=0;halo_lists_LR(1:nhalo_new,k) = halo_list_tmp;nhalo_LR(k) = nhalo_new
              deallocate(halo_list_tmp)
              deallocate(halo_essential)                                         
           end if
        else !! EVEN X: RECEIVE FIRST, SEND SECOND
           if(iproc_S_LR(k).ge.0) then !! Receive neighbour exists
              allocate(halo_essential(nhalo_LR(k)))
              tag = 100*k+iproc_S_LR(k)
              call MPI_RECV(halo_essential,nhalo_LR(k),MPI_INT,iproc_S_LR(k),tag,MPI_COMM_WORLD, &
                            MPI_STATUS_IGNORE,ierror)

              !! Count the new halo size
              nhalo_new = 0
              do i=1,nhalo_LR(k)
                 if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
              end do
              write(6,*) iproc,"Reduced halo_LR",k," from",nhalo_LR(k),"to",nhalo_new              
              
              !! Build the new halo lists...
              allocate(halo_list_tmp(nhalo_new))
              j=0
              do i=1,nhalo_LR(k)
                 if(halo_essential(i).eq.1)then
                    j=j+1
                    halo_list_tmp(j) = halo_lists_LR(i,k)
                 end if
              end do
              
              !! Move new halo list to old halo list
              halo_lists_LR(:,k)=0;halo_lists_LR(1:nhalo_new,k) = halo_list_tmp;nhalo_LR(k) = nhalo_new
              deallocate(halo_list_tmp)
              deallocate(halo_essential)                                         
           end if
           if(iproc_R_LR(k).ge.0)then !! This neighbour sent something
              !! Build halo_essential_from ready for transfer
              allocate(halo_essential_from(inhalo_LR(k)))
              halo_essential_from = halo_essential_all(nrecstart_UDLR(k+2):nrecstart_UDLR(k+2)+inhalo_LR(k)-1)

              !! Send the data back to originator processor
              tag = 100*k+iproc
              call MPI_SEND(halo_essential_from,inhalo_LR(k),MPI_INT,iproc_R_LR(k), &
                            tag,MPI_COMM_WORLD,ierror)
              deallocate(halo_essential_from)
           end if
        end if
     end do     
     
         
     !! Re-set np to without halos
     np = np_nohalo
     deallocate(nrecstart_UDLR)
        
     !! Exchange halo sizes and node lists between processors
     call send_halo_sizes_nodes 
    
     suminhalo = sum(inhalo_LR(1:2*nprocsY)) + sum(inhalo_UD(1:2))         
     np = np_nohalo + suminhalo     
        
     write(6,*) "New halos built",iproc,npfb,np_nohalo,np
  
     return
  end subroutine refine_halos
!! ------------------------------------------------------------------------------------------------  
  subroutine send_halo_sizes_nodes
     !! This subroutine exchanges the sizes of halos, and the lists of halo node positions between
     !! processors
     integer(ikind) :: i,k,nrec_end,is,ie
   
   
     !! UP-DOWN
     do k=1,2
       
        !! Check whether k-th SEND exists
        if(iproc_S_UD(k).ge.0) then
           !! Send the halo size
           call MPI_SEND(nhalo_UD(k),1,MPI_INT,iproc_S_UD(k),1000*k+iproc,MPI_COMM_WORLD,ierror)
        end if
        !! Check whether the k-th RECEIVE neighbour exists
        if(iproc_R_UD(k).ge.0)then
           !! Recieve the size
           call MPI_RECV(inhalo_UD(k),1,MPI_INT,iproc_R_UD(k),1000*k+iproc_R_UD(k), & 
                         MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)        
        end if         
     
     end do        

     write(6,*) iproc,"UD-out",nhalo_UD(:)
     write(6,*) iproc,"UD-in",inhalo_UD(:)
    

     !! LEFT-RIGHT
     do k=1,2*nprocsY
       
        !! Check whether k-th SEND exists
        if(iproc_S_LR(k).ge.0) then
           !! Send the halo size
           call MPI_SEND(nhalo_LR(k),1,MPI_INT,iproc_S_LR(k),1000*k+iproc,MPI_COMM_WORLD,ierror)
        end if
        !! Check whether the k-th RECEIVE neighbour exists
        if(iproc_R_LR(k).ge.0)then
           !! Recieve the size
           call MPI_RECV(inhalo_LR(k),1,MPI_INT,iproc_R_LR(k),1000*k+iproc_R_LR(k), & 
                         MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)        
        end if         
     
     end do        
     
            
     write(6,*) iproc,"LR-out",nhalo_LR(:)
     write(6,*) iproc,"LR-in",inhalo_LR(:)
     
     
     !! Mark out any send & receive processors if halo has zero size
     !! Also build the starting index for recieves
     allocate(nrecstart_UDLR(2+2*nprocsY))
     nrec_end = np_nohalo
     do k=1,2
        if(nhalo_UD(k).eq.0) iproc_S_UD(k)=-1
        if(inhalo_UD(k).eq.0) iproc_R_UD(k)=-1
        nrecstart_UDLR(k) = nrec_end + 1
        nrec_end = nrec_end + inhalo_UD(k)
     end do
     do k=1,2*nprocsY
        if(nhalo_LR(k).eq.0) iproc_S_LR(k)=-1
        if(inhalo_LR(k).eq.0) iproc_R_LR(k)=-1
        nrecstart_UDLR(2+k) = nrec_end + 1
        nrec_end = nrec_end + inhalo_LR(k)
     end do

     !! Loop over all dimensions and exchange node positions
#ifdef dim3
     do i=1,3
#else
     do i=1,2
#endif     
        call halo_exchange(rp(:,i))
     end do

     !! Adjust positions for X-periodicity
     if(xbcond.eq.1.and.nprocsX.gt.1) then
        if(iprocX.eq.0) then
           do k=2+nprocsY+1,2+2*nprocsY
              is = nrecstart_UDLR(k);ie = is + inhalo_LR(k-2)-1
              rp(is:ie,1) = rp(is:ie,1) - (xmax-xmin)
           end do
        end if           
        if(iprocX.eq.nprocsX-1) then
           do k=2+1,2+nprocsY        
              is = nrecstart_UDLR(k);ie = is + inhalo_LR(k-2)-1
              rp(is:ie,1) = rp(is:ie,1) + (xmax-xmin)
           end do
        end if   
     end if    

    
     return
  end subroutine send_halo_sizes_nodes
!! ------------------------------------------------------------------------------------------------  
#endif  
end module mpi_transfers

