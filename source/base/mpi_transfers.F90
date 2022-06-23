module mpi_transfers
  !! This module contains routines to construct halos and do parallel data transfers
  !! ----------------------------------------------------------------------------------------------
  !! 2D DOMAIN DECOMPOSITION in x-y. Each processor has up to 8 neighbour processors
  !!
  !!     SEND ORDER       RECEIVE ORDER
  !!  ---------------   ---------------
  !!    8  | 7 | 6        4  | 3 | 2
  !!  ---------------   ---------------
  !!    1  |   | 5        5  |   | 1 
  !!  ---------------   ---------------
  !!    2  | 3 | 4        6  | 7 | 8  
  !!  ---------------   ---------------
  !!
  !! Transfers are SEND-RECEIVE for ODD processors, and RECEIVE-SEND for EVEN processors in EACH
  !! direction. Hence, if any periodic boundaries are required, need either 1 or EVEN number of 
  !! processors in periodic direction.
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
#ifdef mp
!! ------------------------------------------------------------------------------------------------
  subroutine processor_mapping
     integer(ikind) :: ii
  
     !! Number of processors in X and Y decomposition - check match with schedule     
     if(nprocsX*nprocsY.ne.nprocs) then
        write(6,*) "ERROR: nprocs doesn't match the decomposition schedule from ishift. STOPPING"
        call MPI_Abort(MPI_COMM_WORLD, ii, ierror)
     end if
 
     !! Indices of this processor in X,Y grid 
     iprocX=iproc/nprocsY   
     iprocY=mod(iproc,nprocsY)
         
     !! Lists for 2 dimensional decomposition - SEND LIST
     iprocNS = -1
     if(nprocsY.ne.1)then
        if(iprocY.ne.0) iprocNS(3) = iproc - 1 !D
        if(iprocY.ne.nprocsY-1) iprocNS(7) = iproc + 1 !U
     end if
     !! Y-periodicity
     if(ybcond.eq.1.and.nprocsY.gt.1) then
        if(iprocY.eq.0) iprocNS(3) = iproc - 1 + nprocsY  !D
        if(iprocY.eq.nprocsY-1) iprocNS(7) = iproc + 1 - nprocsY   !U
     end if
     if(iprocX.ne.0) then
        iprocNS(1) = iproc-nprocsY !L
        if(iprocNS(3).ne.-1) iprocNS(2) = iprocNS(3) - nprocsY  !LD
        if(iprocNS(7).ne.-1) iprocNS(8) = iprocNS(7) - nprocsY  !LU
     end if
     if(iprocX.ne.nprocsX-1) then
        iprocNS(5) = iproc+nprocsY !R
        if(iprocNS(3).ne.-1) iprocNS(4) = iprocNS(3) + nprocsY  !RD
        if(iprocNS(7).ne.-1) iprocNS(6) = iprocNS(7) + nprocsY  !RU
     end if
     !! X-periodicity
     if(xbcond.eq.1.and.nprocsX.gt.1) then
        if(iprocX.eq.0) then
           iprocNS(1) = iproc + nprocs - nprocsY !L
           if(iprocNS(3).ne.-1) iprocNS(2) = iprocNS(3) + nprocs - nprocsY !LD
           if(iprocNS(7).ne.-1) iprocNS(8) = iprocNS(7) + nprocs - nprocsY !LU
        end if
        if(iprocX.eq.nprocsX-1) then
           iprocNS(5) = iproc - nprocs + nprocsY
           if(iprocNS(3).ne.-1) iprocNS(4) = iprocNS(3) - nprocs + nprocsY !RD
           if(iprocNS(7).ne.-1) iprocNS(6) = iprocNS(7) - nprocs + nprocsY !RU
        end if     
     end if



     !! Lists for 2 dimensional decomposition - RECEIVE LIST
     iprocNR = -1
     if(nprocsY.ne.1)then
        if(iprocY.ne.0) iprocNR(7) = iproc - 1 !D
        if(iprocY.ne.nprocsY-1) iprocNR(3) = iproc + 1 !U
     end if
     !! Y-periodicity
     if(ybcond.eq.1.and.nprocsY.gt.1) then
        if(iprocY.eq.0) iprocNR(7) = iproc - 1 + nprocsY    !D
        if(iprocY.eq.nprocsY-1) iprocNR(3) = iproc + 1 - nprocsY         !U
     end if         
     if(iprocX.ne.0) then
        iprocNR(5) = iproc-nprocsY
        if(iprocNR(7).ne.-1) iprocNR(6) = iprocNR(7) - nprocsY  !LD
        if(iprocNR(3).ne.-1) iprocNR(4) = iprocNR(3) - nprocsY  !LU
     end if
     if(iprocX.ne.nprocsX-1) then
        iprocNR(1) = iproc+nprocsY
        if(iprocNR(7).ne.-1) iprocNR(8) = iprocNR(7) + nprocsY  !RD
        if(iprocNR(3).ne.-1) iprocNR(2) = iprocNR(3) + nprocsY  !RU
     end if 
     !! X-periodicity
     if(xbcond.eq.1.and.nprocsX.gt.1) then
        if(iprocX.eq.0) then
           iprocNR(5) = iproc + nprocs - nprocsY !L
           if(iprocNR(7).ne.-1) iprocNR(6) = iprocNR(7) + nprocs - nprocsY !LD
           if(iprocNR(3).ne.-1) iprocNR(4) = iprocNR(3) + nprocs - nprocsY !LU
        end if
        if(iprocX.eq.nprocsX-1) then
           iprocNR(1) = iproc - nprocs + nprocsY
           if(iprocNR(7).ne.-1) iprocNR(8) = iprocNR(7) - nprocs + nprocsY !RD
           if(iprocNR(3).ne.-1) iprocNR(2) = iprocNR(3) - nprocs + nprocsY !RU
        end if     
     end if
     
     write(6,*) "process",iproc,"NS:",iprocNS(:)
     write(6,*) "process",iproc,"NR:",iprocNR(:)     

     
     return
  end subroutine processor_mapping
!! ------------------------------------------------------------------------------------------------
  subroutine build_halos
     integer(ikind) :: i,nstart,nend,maxhalo,suminhalo,is,ie,k
     real(rkind) :: x,y,ss_local,halo_fac
     integer(ikind),dimension(:,:),allocatable :: halo_lists_tmp
     !! Routine builds a list of nodes which are to be exported as halos for adjacent processors
          
     !! How much extra?
     halo_fac = 6.0d0
      
     !! Allocate a temporary halo list
     allocate(halo_lists_tmp(np,8))
     nhalo=0
     
     !! Build halo LEFT    
     if(iprocNS(1).ge.0) then           
        do i=1,np
           x=rp(i,1);ss_local = halo_fac*ss*h(i)      
           if(abs(x-XL_thisproc).le.ss_local) then     
              nhalo(1) = nhalo(1) + 1
              halo_lists_tmp(nhalo(1),1) = i
           end if
        end do
     end if
              
     !! Build halo RIGHT
     if(iprocNS(5).ge.0) then  
        do i=1,np
           x=rp(i,1);ss_local = halo_fac*ss*h(i)
           if(abs(x-XR_thisproc).le.ss_local) then
              nhalo(5) = nhalo(5) + 1
              halo_lists_tmp(nhalo(5),5) = i              
           end if
        end do
     end if
     
     !! Build halo UP
     if(iprocNS(7).ge.0) then
        do i=1,np
           y=rp(i,2);ss_local = halo_fac*ss*h(i)      
           if(abs(y-YU_thisproc).le.ss_local) then     
              nhalo(7) = nhalo(7) + 1
              halo_lists_tmp(nhalo(7),7) = i              
           end if
        end do
     end if

     !! Build halo DOWN     
     if(iprocNS(3).ge.0)then            
        do i=1,np
           y=rp(i,2);ss_local = halo_fac*ss*h(i)
           if(abs(y-YD_thisproc).le.ss_local) then
              nhalo(3) = nhalo(3) + 1
              halo_lists_tmp(nhalo(3),3) = i              
           end if
        end do           
     end if
     
     !! Build halo LEFT-DOWN 
     if(iprocNS(2).ge.0) then           
        do i=1,np
           x=rp(i,1);y=rp(i,2);ss_local = halo_fac*ss*h(i)     
!           if(abs(x-XL_thisproc).le.ss_local.and.abs(y-YD_thisproc).le.ss_local) then     
           if(abs(x-XL_thisproc).le.ss_local)then
              nhalo(2) = nhalo(2) + 1
              halo_lists_tmp(nhalo(2),2) = i
           end if
        end do
     end if     

     !! Build halo RIGHT-DOWN 
     if(iprocNS(4).ge.0) then           
        do i=1,np
           x=rp(i,1);y=rp(i,2);ss_local = halo_fac*ss*h(i)      
!           if(abs(x-XR_thisproc).le.ss_local.and.abs(y-YD_thisproc).le.ss_local) then     
           if(abs(x-XR_thisproc).le.ss_local)then
              nhalo(4) = nhalo(4) + 1
              halo_lists_tmp(nhalo(4),4) = i
           end if
        end do
     end if     
     
     !! Build halo RIGHT-UP
     if(iprocNS(6).ge.0) then           
        do i=1,np
           x=rp(i,1);y=rp(i,2);ss_local = halo_fac*ss*h(i)    
!           if(abs(x-XR_thisproc).le.ss_local.and.abs(y-YU_thisproc).le.ss_local) then     
           if(abs(x-XR_thisproc).le.ss_local)then
              nhalo(6) = nhalo(6) + 1
              halo_lists_tmp(nhalo(6),6) = i
           end if
        end do
     end if          
     
     !! Build halo LEFT-UP 
     if(iprocNS(8).ge.0) then           
        do i=1,np
           x=rp(i,1);y=rp(i,2);ss_local = halo_fac*ss*h(i)     
!           if(abs(x-XL_thisproc).le.ss_local.and.abs(y-YU_thisproc).le.ss_local) then     
           if(abs(x-XL_thisproc).le.ss_local)then           
              nhalo(8) = nhalo(8) + 1
              halo_lists_tmp(nhalo(8),8) = i
           end if
        end do
     end if          
     
     !! Store halos in final array
     maxhalo = maxval(nhalo(1:8))
     allocate(halo_lists(maxhalo,8));halo_lists=0
    
     do k=1,8
        if(iprocNS(k).ne.-1)then
           halo_lists(1:nhalo(k),k) = halo_lists_tmp(1:nhalo(k),k)
        end if
     end do
     
     deallocate(halo_lists_tmp)
      
             
     !! Exchange halo sizes and node positions between processors   
     call send_halo_sizes_nodes
     
    

     !! Update size and node-types
     suminhalo = sum(inhalo(1:8))
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

     !! Loop over each possible neighbour
     do k =1,8

        if(k.ne.3.and.k.ne.7) then !! Transfers to left and right plus diagonals (1,2,4,5,6,8)
     
           !! XODD: 
           if(xodd) then  !! SEND LEFT, RECEIVE FROM RIGHT   

              if(iprocNS(k).ge.0)then !! Send neighbour exists
                 allocate(halo_phi(nhalo(k)))
                 do i=1,nhalo(k)
                    j=halo_lists(i,k)
                   halo_phi(i) = phi(j)
                 end do

                 !! Send the data
                 tag = iproc + 100*k
                 call MPI_SEND(halo_phi,nhalo(k),MPI_DOUBLE_PRECISION,iprocNS(k),tag,MPI_COMM_WORLD,ierror)
                 deallocate(halo_phi) 
              end if   
              if(iprocNR(k).ge.0)then !! Receive neighbour exists
                 !! Receive the data
                 tag = iprocNR(k) + 100*k
                 call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_DOUBLE_PRECISION &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              end if  
    
           !! XEVEN:
           else     !! RECEIVE FROM RIGHT, SEND LEFT
              if(iprocNR(k).ge.0)then !! Receive neighbour exists
                 !! Receive the data
                 tag = iprocNR(k) + 100*k
                 call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_DOUBLE_PRECISION &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              end if  
              if(iprocNS(k).ge.0)then !! Send neighbour exists
                 allocate(halo_phi(nhalo(k)))
                 do i=1,nhalo(k)
                    j=halo_lists(i,k)
                    halo_phi(i) = phi(j)
                 end do

                 !! Send the data
                 tag = iproc + 100*k
                 call MPI_SEND(halo_phi,nhalo(k),MPI_DOUBLE_PRECISION,iprocNS(k),tag,MPI_COMM_WORLD,ierror) 
                 deallocate(halo_phi)            
              end if               
           end if
        
        else  !! Transfers up and down (3 and 7)
           if(yodd) then !! SEND DOWN, RECEIVE FROM UP
              if(iprocNS(k).ge.0)then !! Send neighbour exists
                 allocate(halo_phi(nhalo(k)))
                 do i=1,nhalo(k)
                    j=halo_lists(i,k)
                   halo_phi(i) = phi(j)
                 end do

                 !! Send the data
                 tag = iproc + 100*k
                 call MPI_SEND(halo_phi,nhalo(k),MPI_DOUBLE_PRECISION,iprocNS(k),tag,MPI_COMM_WORLD,ierror)
                 deallocate(halo_phi) 
              end if   
              if(iprocNR(k).ge.0)then !! Receive neighbour exists
                 !! Receive the data
                 tag = iprocNR(k) + 100*k
                 call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_DOUBLE_PRECISION &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              end if           
                     
           else          !! RECEIVE FROM UP, SEND DOWN
              if(iprocNR(k).ge.0)then !! Receive neighbour exists
                 !! Receive the data
                 tag = iprocNR(k) + 100*k
                 call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_DOUBLE_PRECISION &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              end if  
              if(iprocNS(k).ge.0)then !! Send neighbour exists
                 allocate(halo_phi(nhalo(k)))
                 do i=1,nhalo(k)
                    j=halo_lists(i,k)
                    halo_phi(i) = phi(j)
                 end do

                 !! Send the data
                 tag = iproc + 100*k
                 call MPI_SEND(halo_phi,nhalo(k),MPI_DOUBLE_PRECISION,iprocNS(k),tag,MPI_COMM_WORLD,ierror) 
                 deallocate(halo_phi)            
              end if                                 
           end if
        end if
     end do
   
     return
  end subroutine halo_exchange
!! ------------------------------------------------------------------------------------------------  
  subroutine halo_exchange_backup(phi)    
     !! This routine does halo exchanges for phi= u,v,lnro,E,Y0
     real(rkind),dimension(:),intent(inout) :: phi     
     real(rkind),dimension(:),allocatable :: halo_phi
     integer(ikind) :: i,j,nrec_start,nrec_end,k,tag

     !! Loop over each possible neighbour processor
     do k=1,8
    
        !! If the k-th send neighbour exists
        if(iprocNS(k).ge.0) then  

           !! Allocate space
           allocate(halo_phi(nhalo(k)))
!           !$omp parallel do private(j)
           do i=1,nhalo(k)
              j=halo_lists(i,k)
              halo_phi(i) = phi(j)
           end do
!           !$omp end parallel do

           !! Send the data
           tag = iproc + 100*k
           call MPI_SEND(halo_phi,nhalo(k),MPI_DOUBLE_PRECISION,iprocNS(k),tag,MPI_COMM_WORLD,ierror)
           
           !! Clear space
           deallocate(halo_phi)
        end if

        
        !! if the k-th receive neighbour exists
        if(iprocNR(k).ge.0)then
                
           !! Receive the data
           tag = iprocNR(k) + 100*k
           call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_DOUBLE_PRECISION &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
                
        end if
                        
     end do
     

     return
  end subroutine halo_exchange_backup 
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

     !! Loop over each possible neighbour
     do k =1,8

        if(k.ne.3.and.k.ne.7) then !! Transfers to left and right plus diagonals (1,2,4,5,6,8)
     
           !! XODD: 
           if(xodd) then  !! SEND LEFT, RECEIVE FROM RIGHT   

              if(iprocNS(k).ge.0)then !! Send neighbour exists
                 allocate(halo_phi(nhalo(k)))
                 do i=1,nhalo(k)
                    j=halo_lists(i,k)
                   halo_phi(i) = phi(j)
                 end do

                 !! Send the data
                 tag = iproc + 100*k
                 call MPI_SEND(halo_phi,nhalo(k),MPI_INT,iprocNS(k),tag,MPI_COMM_WORLD,ierror)
                 deallocate(halo_phi) 
              end if   
              if(iprocNR(k).ge.0)then !! Receive neighbour exists
                 !! Receive the data
                 tag = iprocNR(k) + 100*k
                 call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_INT &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              end if  
    
           !! XEVEN:
           else     !! RECEIVE FROM RIGHT, SEND LEFT
              if(iprocNR(k).ge.0)then !! Receive neighbour exists
                 !! Receive the data
                 tag = iprocNR(k) + 100*k
                 call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_INT &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              end if  
              if(iprocNS(k).ge.0)then !! Send neighbour exists
                 allocate(halo_phi(nhalo(k)))
                 do i=1,nhalo(k)
                    j=halo_lists(i,k)
                    halo_phi(i) = phi(j)
                 end do

                 !! Send the data
                 tag = iproc + 100*k
                 call MPI_SEND(halo_phi,nhalo(k),MPI_INT,iprocNS(k),tag,MPI_COMM_WORLD,ierror) 
                 deallocate(halo_phi)            
              end if               
           end if
        
        else  !! Transfers up and down (3 and 7)
           if(yodd) then !! SEND DOWN, RECEIVE FROM UP
              if(iprocNS(k).ge.0)then !! Send neighbour exists
                 allocate(halo_phi(nhalo(k)))
                 do i=1,nhalo(k)
                    j=halo_lists(i,k)
                   halo_phi(i) = phi(j)
                 end do

                 !! Send the data
                 tag = iproc + 100*k
                 call MPI_SEND(halo_phi,nhalo(k),MPI_INT,iprocNS(k),tag,MPI_COMM_WORLD,ierror)
                 deallocate(halo_phi) 
              end if   
              if(iprocNR(k).ge.0)then !! Receive neighbour exists
                 !! Receive the data
                 tag = iprocNR(k) + 100*k
                 call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_INT &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              end if           
                     
           else          !! RECEIVE FROM UP, SEND DOWN
              if(iprocNR(k).ge.0)then !! Receive neighbour exists
                 !! Receive the data
                 tag = iprocNR(k) + 100*k
                 call MPI_RECV(phi(nrecstart(k):nrecstart(k)+inhalo(k)-1),inhalo(k),MPI_INT &
                         ,iprocNR(k),tag,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)
              end if  
              if(iprocNS(k).ge.0)then !! Send neighbour exists
                 allocate(halo_phi(nhalo(k)))
                 do i=1,nhalo(k)
                    j=halo_lists(i,k)
                    halo_phi(i) = phi(j)
                 end do

                 !! Send the data
                 tag = iproc + 100*k
                 call MPI_SEND(halo_phi,nhalo(k),MPI_INT,iprocNS(k),tag,MPI_COMM_WORLD,ierror) 
                 deallocate(halo_phi)            
              end if                                 
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
     integer(ikind) i,j,k,jstart,jend,maxhalo,suminhalo,nhalo_new
     logical :: xodd,yodd

     !! Identify required halo nodes
     allocate(halo_essential_all(np));halo_essential_all =0
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

     jend=np_nohalo
     !! Loop over all RECEIVE neighbours
     do k=1,8     
        if(k.ne.3.and.k.ne.7) then
           if(xodd) then !! SEND LEFT, RECEIVE FROM RIGHT        
              !! Check whether this receive neighbour sent anything
              if(iprocNR(k).ge.0)then        
                 !! Build halo_essential_from ready for transfer
                 allocate(halo_essential_from(inhalo(k)))
                 halo_essential_from = halo_essential_all(nrecstart(k):nrecstart(k)+inhalo(k)-1)
                 !! Send back to originator processor
                 call MPI_SEND(halo_essential_from,inhalo(k),MPI_INT,iprocNR(k),100*k,MPI_COMM_WORLD,ierror)          
                 deallocate(halo_essential_from)
              end if        
              !! Check whether this was an originator processor 
              if(iprocNS(k).ge.0)then

                 !! Receive the essential list
                 allocate(halo_essential(nhalo(k)))
                 call MPI_RECV(halo_essential,nhalo(k),MPI_INT,iprocNS(k),100*k,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)

                 !! Count the new halo size
                 nhalo_new = 0
                 do i=1,nhalo(k)
                    if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
                 end do
           
                 write(6,*) iproc,"Reduced halo",k," from",nhalo(k),"to",nhalo_new
                 !! Build the new halo lists...
                 allocate(halo_list_tmp(nhalo_new))
                 j=0
                 do i=1,nhalo(k)
                    if(halo_essential(i).eq.1)then
                       j=j+1
                       halo_list_tmp(j) = halo_lists(i,k)
                    end if
                 end do
           
                 !! Move new halo list to old halo list
                 halo_lists(:,k)=0;halo_lists(1:nhalo_new,k) = halo_list_tmp;nhalo(k) = nhalo_new           
                 deallocate(halo_list_tmp,halo_essential)                   
              end if
          else  !! RECEIVE FROM RIGHT, SEND LEFT
              !! Check whether this was an originator processor 
              if(iprocNS(k).ge.0)then

                 !! Receive the essential list
                 allocate(halo_essential(nhalo(k)))
                 call MPI_RECV(halo_essential,nhalo(k),MPI_INT,iprocNS(k),100*k,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)

                 !! Count the new halo size
                 nhalo_new = 0
                 do i=1,nhalo(k)
                    if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
                 end do
           
                 write(6,*) iproc,"Reduced halo",k," from",nhalo(k),"to",nhalo_new
                 !! Build the new halo lists...
                 allocate(halo_list_tmp(nhalo_new))
                 j=0
                 do i=1,nhalo(k)
                    if(halo_essential(i).eq.1)then
                       j=j+1
                       halo_list_tmp(j) = halo_lists(i,k)
                    end if
                 end do
           
                 !! Move new halo list to old halo list
                 halo_lists(:,k)=0;halo_lists(1:nhalo_new,k) = halo_list_tmp;nhalo(k) = nhalo_new           
                 deallocate(halo_list_tmp,halo_essential)                   
              end if          
              !! Check whether this receive neighbour sent anything
              if(iprocNR(k).ge.0)then        
                 !! Build halo_essential_from ready for transfer
                 allocate(halo_essential_from(inhalo(k)))
                 halo_essential_from = halo_essential_all(nrecstart(k):nrecstart(k)+inhalo(k)-1)
                 !! Send back to originator processor
                 call MPI_SEND(halo_essential_from,inhalo(k),MPI_INT,iprocNR(k),100*k,MPI_COMM_WORLD,ierror)          
                 deallocate(halo_essential_from)
              end if                                                    
          end if
       else ! Transfers up and down (3 and 7)
           if(yodd) then !! SEND DOWN, RECEIVE FROM UP
              !! Check whether this receive neighbour sent anything
              if(iprocNR(k).ge.0)then        
                 !! Build halo_essential_from ready for transfer
                 allocate(halo_essential_from(inhalo(k)))
                 halo_essential_from = halo_essential_all(nrecstart(k):nrecstart(k)+inhalo(k)-1)
                 !! Send back to originator processor
                 call MPI_SEND(halo_essential_from,inhalo(k),MPI_INT,iprocNR(k),100*k,MPI_COMM_WORLD,ierror)          
                 deallocate(halo_essential_from)
              end if        
              !! Check whether this was an originator processor 
              if(iprocNS(k).ge.0)then

                 !! Receive the essential list
                 allocate(halo_essential(nhalo(k)))
                 call MPI_RECV(halo_essential,nhalo(k),MPI_INT,iprocNS(k),100*k,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)

                 !! Count the new halo size
                 nhalo_new = 0
                 do i=1,nhalo(k)
                    if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
                 end do
           
                 write(6,*) iproc,"Reduced halo",k," from",nhalo(k),"to",nhalo_new
                 !! Build the new halo lists...
                 allocate(halo_list_tmp(nhalo_new))
                 j=0
                 do i=1,nhalo(k)
                    if(halo_essential(i).eq.1)then
                       j=j+1
                       halo_list_tmp(j) = halo_lists(i,k)
                    end if
                 end do
           
                 !! Move new halo list to old halo list
                 halo_lists(:,k)=0;halo_lists(1:nhalo_new,k) = halo_list_tmp;nhalo(k) = nhalo_new           
                 deallocate(halo_list_tmp,halo_essential)                   
              end if
          else  !! RECEIVE FROM UP, SEND DOWN
              !! Check whether this was an originator processor 
              if(iprocNS(k).ge.0)then

                 !! Receive the essential list
                 allocate(halo_essential(nhalo(k)))
                 call MPI_RECV(halo_essential,nhalo(k),MPI_INT,iprocNS(k),100*k,MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)

                 !! Count the new halo size
                 nhalo_new = 0
                 do i=1,nhalo(k)
                    if(halo_essential(i).eq.1) nhalo_new = nhalo_new + 1
                 end do
           
                 write(6,*) iproc,"Reduced halo",k," from",nhalo(k),"to",nhalo_new
                 !! Build the new halo lists...
                 allocate(halo_list_tmp(nhalo_new))
                 j=0
                 do i=1,nhalo(k)
                    if(halo_essential(i).eq.1)then
                       j=j+1
                       halo_list_tmp(j) = halo_lists(i,k)
                    end if
                 end do
           
                 !! Move new halo list to old halo list
                 halo_lists(:,k)=0;halo_lists(1:nhalo_new,k) = halo_list_tmp;nhalo(k) = nhalo_new           
                 deallocate(halo_list_tmp,halo_essential)                   
              end if          
              !! Check whether this receive neighbour sent anything
              if(iprocNR(k).ge.0)then        
                 !! Build halo_essential_from ready for transfer
                 allocate(halo_essential_from(inhalo(k)))
                 halo_essential_from = halo_essential_all(nrecstart(k):nrecstart(k)+inhalo(k)-1)
                 !! Send back to originator processor
                 call MPI_SEND(halo_essential_from,inhalo(k),MPI_INT,iprocNR(k),100*k,MPI_COMM_WORLD,ierror)          
                 deallocate(halo_essential_from)
              end if                                                    
          end if             
       end if  
     end do
        
     !! Re-set np to without halos
     np = np_nohalo
        
     !! Exchange halo sizes and node lists between processors
     call send_halo_sizes_nodes
     
     !! Update size and node-types
     suminhalo = sum(inhalo(1:8))
     np = np_nohalo + suminhalo

        
     write(6,*) "New halos built",iproc,npfb,np_nohalo,np
  
     return
  end subroutine refine_halos
!! ------------------------------------------------------------------------------------------------  
  subroutine send_halo_sizes_nodes
     !! This subroutine exchanges the sizes of halos, and the lists of halo node positions between
     !! processors
     integer(ikind) :: i,k,nrec_end,is,ie
   
     inhalo = 0       
     !! Loop over all send neighbour processors
     do k=1,8

        !! Check whether the k-th SEND neighbour exists
        if(iprocNS(k).ge.0)then
           !! Send the size
           call MPI_SEND(nhalo(k),1,MPI_INT,iprocNS(k),1000*k+iproc,MPI_COMM_WORLD,ierror)                      
        end if
        
        !! Check whether the k-th RECEIVE neighbour exists
        if(iprocNR(k).ge.0)then
           !! Recieve the size
           call MPI_RECV(inhalo(k),1,MPI_INT,iprocNR(k),1000*k+iprocNR(k),MPI_COMM_WORLD,MPI_STATUS_IGNORE,ierror)        
        end if         
     
     end do
           
     write(6,*) "exchanged nhalo sizes",iproc,nhalo     
     write(6,*) "exchanged inhalo sizes",iproc,inhalo  
     
     !! If the k-th halo has size zero, mark out the k-th send/receive processors
     !! Also build the starting index for recieves
     nrec_end = np_nohalo
     do k=1,8
        if(nhalo(k).eq.0) iprocNS(k)=-1
        if(inhalo(k).eq.0) iprocNR(k)=-1
        nrecstart(k) = nrec_end + 1;
        nrec_end = nrec_end + inhalo(k)
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
           is = nrecstart(5);ie = is + inhalo(5)-1
           rp(is:ie,1) = rp(is:ie,1) - (xmax-xmin)

           if(iprocNR(4).ne.-1) then
              is = nrecstart(4);ie = is + inhalo(4)-1
              rp(is:ie,1) = rp(is:ie,1) - (xmax-xmin)
           end if
           if(iprocNR(6).ne.-1) then
              is = nrecstart(6);ie = is + inhalo(6)-1
              rp(is:ie,1) = rp(is:ie,1) - (xmax-xmin)
           end if
        end if
        if(iprocX.eq.nprocsX-1) then
           is = nrecstart(1);ie = is + inhalo(1)-1
           rp(is:ie,1) = rp(is:ie,1) + (xmax-xmin)

           if(iprocNR(2).ne.-1) then
              is = nrecstart(2);ie = is + inhalo(2)-1
              rp(is:ie,1) = rp(is:ie,1) + (xmax-xmin)
           end if
           if(iprocNR(8).ne.-1) then
              is = nrecstart(8);ie = is + inhalo(8)-1
              rp(is:ie,1) = rp(is:ie,1) + (xmax-xmin)
           end if
        end if    
     end if     
     !! Adjust positions for Y-periodicity
     if(ybcond.eq.1.and.nprocsY.gt.1) then
        if(iprocY.eq.0) then
           is = nrecstart(7);ie = is + inhalo(7)-1
           rp(is:ie,2) = rp(is:ie,2) - (ymax-ymin)
           
           if(iprocNR(6).ne.-1) then
              is = nrecstart(6);ie = is + inhalo(6)-1
              rp(is:ie,2) = rp(is:ie,2) - (ymax-ymin)
           end if
           if(iprocNR(8).ne.-1) then
              is = nrecstart(8);ie = is + inhalo(8)-1
              rp(is:ie,2) = rp(is:ie,2) - (ymax-ymin)           
           end if
        end if
        if(iprocY.eq.nprocsY-1) then
           is = nrecstart(3);ie = is + inhalo(3)-1
           rp(is:ie,2) = rp(is:ie,2) + (ymax-ymin)
           
           if(iprocNR(2).ne.-1) then
              is = nrecstart(2);ie = is + inhalo(2)-1
              rp(is:ie,2) = rp(is:ie,2) + (ymax-ymin)           
           end if
           if(iprocNR(4).ne.-1) then
              is = nrecstart(4);ie = is + inhalo(4)-1
              rp(is:ie,2) = rp(is:ie,2) + (ymax-ymin)                      
           end if
        end if    
     end if   
               
     return
  end subroutine send_halo_sizes_nodes
!! ------------------------------------------------------------------------------------------------  
#endif  
end module mpi_transfers
