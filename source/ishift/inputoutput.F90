module inputoutput
  !! This module contains routines to read in node distributions, modify them with shifting 
  !! pre-process to obtain boundary normal vectors, and write field data out to file.
  use kind_parameters
  use common_parameter
  use common_2d
  use omp_lib
  use neighbours
  implicit none

  real(rkind),dimension(:),allocatable :: x,y,xn,yn,ds  !! Properties for post-shift tidying
  integer(ikind),dimension(:),allocatable :: nt  

  !! Start and end indices for domain decomposition
  integer(ikind),dimension(:),allocatable :: nstart,nend

contains
!! ------------------------------------------------------------------------------------------------  
  subroutine iteratively_shift(kk)
     use boundaries     
     !! Subroutine to create a nice shifted distribution...
     !! N.B. Because the shifting is only based on rij<2s (i.e. half the stencil), we can get away
     !! with not updating neighbour lists or re-doing boundaries within each iteration. Makes life
     !! easier for the multi-process version.
     integer(ikind),intent(in) :: kk
     integer(ikind) :: ll,i,j,k
     real(rkind),dimension(dims) :: rij,gradw,dr_tmp
     real(rkind),dimension(:,:),allocatable :: dr
     real(rkind) :: rad,tmp,qq
     real(rkind) :: qkd_mag,ns,drmag
    
     allocate(dr(npfb,dims))
 
     if(allocated(irelation)) deallocate(irelation,vrelation)
     call create_mirror_particles
     if(allocated(ij_count)) deallocate(ij_count,ij_link)
     call find_neighbours


     !! Low order shifting loop
     do ll=1,kk

!        if(allocated(irelation)) deallocate(irelation,vrelation)
!        call create_mirror_particles
!        if(allocated(ij_count)) deallocate(ij_count,ij_link)
!        call find_neighbours


      
        !! Find shifting vector...
        !$OMP PARALLEL DO PRIVATE(k,j,rij,rad,qq,gradw,dr_tmp,qkd_mag)
        do i=1,npfb
           if(node_type(i).eq.999) then
              qkd_mag = 1.0d-1*h(i)
              dr_tmp = zero
              do k=1,ij_count(i)
                 j=ij_link(i,k)
                 rij = rp(i,:)-rp(j,:);rad=sqrt(dot_product(rij,rij));qq=hovs*rad/h(i)

                 gradw(:) = qkd_mag*((half*qq - one)**1.0)*rij(:)/max(rad,epsilon(rad))
                 if(qq.gt.2.0) gradw(:) = zero            
                 dr_tmp = dr_tmp + gradw(:)
              end do
              dr(i,:) = dr_tmp
              rad = sqrt(dot_product(dr_tmp,dr_tmp))

              if(rad.gt.0.1d0*h(i))then
                 dr(i,:) = dr_tmp*0.1d0*h(i)/rad
              end if
           end if

        end do
        !$OMP END PARALLEL DO
        
        
        !! Move particles...
        !$OMP PARALLEL DO
        do i=1,npfb
           if(node_type(i).eq.999) then
              rp(i,:) = rp(i,:) - dr(i,:)
           end if
        end do
        !$OMP END PARALLEL DO

        !! Move mirrors        
        !$omp parallel do private(j,k) 
        do i=npfb+1,np
           j=irelation(i)
           k=vrelation(i)
           if(k.eq.1) then
              rp(i,:) = rp(i,:) - dr(j,:)
           end if
           if(k.eq.2) then
              rp(i,1) = rp(i,1) + dr(j,1);rp(i,2) = rp(i,2) - dr(j,2)
           end if
           if(k.eq.3) then
              rp(i,1) = rp(i,1) - dr(j,1);rp(i,2) = rp(i,2) + dr(j,2)
           end if
           if(k.eq.4) then
              rp(i,:) = rp(i,:) + dr(j,:)
           end if
        end do
        !$omp end parallel do

write(6,*) "Shifting iteration",ll,"of ",kk

     end do
     deallocate(ij_link,ij_count)     


     deallocate(dr)
     return
  end subroutine iteratively_shift
!! ------------------------------------------------------------------------------------------------
  subroutine setup_domain
     !! Reads in boundary patches
     use boundaries
     integer(ikind) i,j,ii,jj,npfb_tmp,k
     real(rkind) :: ns,dummy,prox,rad,radmin,dx,dy,smag
     real(rkind),dimension(dims) :: rij
     real(rkind),dimension(:,:),allocatable :: tmp_vec

     !! STEP 1: Load IPART (some params, plus list of nodes + boundary normals)
     open(13,file='../gen/IPART')
     read(13,*) nb,npfb,dummy      !! dummy is largest s(i) in domain...
     read(13,*) xmin,xmax,ymin,ymax
     read(13,*) xbcond,ybcond
     !! Calculate some useful constants
     smax = dummy;h0 = hovs*dummy;sup_size = ss*h0;h2=h0*h0;h3=h2*h0
        
    
     allocate(rp(4*npfb,dims),rnorm(4*npfb,dims),h(4*npfb),s(4*npfb));rp=0.0d0;rnorm=0.0d0
     allocate(node_type(4*npfb));node_type=0
     allocate(fd_parent(2*npfb));fd_parent=0
          
     !! Load all nodes. Build FD stencils near boundaries on the fly.
     npfb_tmp = npfb
     nb = 0;ii = 0
    
     do i=1,npfb_tmp
        ii = ii + 1
        read(13,*) rp(ii,1:2),jj,rnorm(ii,1:2),dummy
        h(ii) = dummy*hovs
        s(ii) = dummy
        node_type(ii) = jj
        if(jj.ge.0.and.jj.le.2) then !! If it is a boundary node
           k = ii !! k is the index of the parent node
           nb = nb + 1
           do j=1,4  !! Make 4 additional nodes
              ii = ii + 1
              rp(ii,:) = rp(k,:) + rnorm(k,:)*dble(j)*s(k)   !! Moving along an FD stencil
              rnorm(ii,:)=rnorm(k,:)          !! Copy normals
              h(ii)=h(k);s(ii)=s(k)          !! length-scales
              node_type(ii) = -j           !! and node type
              fd_parent(ii) = k            !! and lineage
              npfb = npfb + 1           
           end do
        end if
     end do

     write(6,*) nb,npfb     
     close(13)     

     !! Randomly perturb the nodes
     smag = 0.5d0
     do i=1,npfb
        if(node_type(i).ge.3)then  !! only perturb interior nodes
        dx = rand();dx = smag*(dx - 0.5d0)*2.0d0*s(i)
        dy = rand();dy = smag*(dy - 0.5d0)*2.0d0*s(i)        
        rp(i,:) = rp(i,:) + (/dx,dy/)
        end if
!        write(32,*) rp(i,:)
     end do

     write(6,*) "Before iterative shifting:",nb,npfb,np

     call iteratively_shift(10)
     call iteratively_shift(10)
     call iteratively_shift(10)
     call iteratively_shift(10)
     call iteratively_shift(10)               
     
     !! Output post-shift distribution
!     do i=1,npfb
!        write(31,*) rp(i,:)
!     end do
     write(6,*) "After iterative shifting:",nb,npfb,np


                    
     return
  end subroutine setup_domain
!! ------------------------------------------------------------------------------------------------  
  subroutine output_newnodes
     integer(ikind) :: i,n,j
     
     n= npfb - 4*nb
     open(212,file='../../IPART')
     write(212,*) nb,n*nprocsZ,smax
     write(212,*) xmin,xmax,ymin,ymax
     write(212,*) xbcond,ybcond
     write(212,*) nprocsX,nprocsY,nprocsZ
  
     !! Indices of each column 
     write(212,*) nprocs*nprocsZ
     do j=1,nprocsZ
        do i=1,nprocs
           write(212,*) (j-1)*n + nstart(i),(j-1)*n + nend(i)
        end do
     end do
     deallocate(nstart,nend)
     
     do j=1,nprocsZ
        do i=1,n
           write(212,*) x(i),y(i),nt(i),xn(i),yn(i),ds(i)
        end do
     end do
     
     deallocate(x,y,xn,yn,ds,nt)
 
     close(212)
     
     write(6,*) "Output written, ending"

     return
  end subroutine output_newnodes
!! ------------------------------------------------------------------------------------------------
  subroutine remove_fd_nodes
     integer(ikind) :: i,n,j
     
     !! Number of nodes without FD stencils
     n = npfb - 4*nb
  
     allocate(x(n),y(n),xn(n),yn(n),ds(n),nt(n))
     
     j=0
     do i=1,npfb      !! Loop over all nodes
        if(node_type(i).ge.0) then !! Exclude FD stencil
           j=j+1

           !! Populate new arrays
           x(j) = rp(i,1);y(j)=rp(i,2)
           xn(j) = rnorm(i,1);yn(j) = rnorm(i,2)
           ds(j) = s(i)
           nt(j) = node_type(i)
        end if
     end do
     
     !! Output a sanity check
     write(6,*) "Expected",n,"found",j
        
     return
  end subroutine remove_fd_nodes  
!! ------------------------------------------------------------------------------------------------
  subroutine rearrange_nodes
     integer(ikind) :: kk,nband,nptmp,nl_ini,nl_end,ii,nblock,ll,nl_end_column
     integer(ikind) :: nshift
   
     !! First sort nodes in X
     write(6,*) "About to quicksort"
     call quicksort(x,1,npfb-4*nb)
     write(6,*) "Quicksorted nodes ordered increasing x"
     
     
    
     !! How many particles (total) need sorting
     nptmp = npfb-4*nb
     
     !! allocation of index limits
     allocate(nstart(nprocs),nend(nprocs))
     
     !! Loop over all bands along X
     do kk=1,nprocsX

           !! Approximate particles per band:
           nband = ceiling(dble(nptmp/nprocsX))
        
           nl_ini = (kk-1)*nband + 1
           nl_end = nl_ini - 1 + nband
     
           !! Add a few to final process if npfb_global isn't divisible by nprocs
           if(kk.eq.nprocsX) then 
           write(6,*) "hi",kk,nprocsX
              ii = nl_end - nptmp   !! How many too many (or too few if negative)
              nband = nband - ii;nl_end = nl_ini - 1 + nband
           end if
   
    
        !! If more than 1 processor in Y direction, for each band of X, sort in Y
        if(nprocsY.ne.1)then
        
    
           !! Sort the range nl_ini:nl_end by y-position            
           write(6,*) "sorting in y for processor band ",kk,nl_ini,nl_end,nband
           call quicksorty(y,nl_ini,nl_end)
           write(6,*) "sorted in y for processor band ",kk   
           
           !! Shuffle the blocks to create cyclical structure in y
           nblock = ceiling(dble(nband/nprocsY))      
           nshift = floor(0.5*nblock) 
           call shift_indices(nl_ini,nl_end,nshift)
                       
                  
        
           !! Calculate the block sizes
           nl_end_column = nl_end
           nl_end = nl_ini - 1
           do ll=1,nprocsY
              if(ll.ne.nprocsY) then
                 nl_ini = nl_end + 1
                 nl_end = nl_ini + nblock              
              else
                 nl_ini = nl_end + 1
                 nl_end = nl_end_column              
              end if
                                            
              !! Store the start and end indices of the block
              nstart((kk-1)*nprocsY+ll) = nl_ini
              nend((kk-1)*nprocsY+ll) = nl_end
           
           end do
        
        
        
        else
           !! Store the start and end indices of the Xbands (each column)
           nstart(kk) = nl_ini
           nend(kk) = nl_end              
                   
                            
        end if
               
     end do
     
      
     return
  end subroutine rearrange_nodes
!! ------------------------------------------------------------------------------------------------ 
  subroutine shift_indices(istart,iend,nswap)
     integer(ikind), intent(in) :: istart,iend,nswap 
     real(rkind),dimension(:),allocatable :: x_tmp,y_tmp,xn_tmp,yn_tmp,ds_tmp
     integer(ikind),dimension(:),allocatable :: nt_tmp
     integer(ikind) :: band_size,shift_size,i_old,i_new,i
     
     !! Sizes
     band_size = 1 + iend - istart
     shift_size = band_size - nswap
     write(6,*) "shift size,bandsize",shift_size,band_size
     
     !! Make some space
     allocate(x_tmp(nswap),y_tmp(nswap),xn_tmp(nswap),yn_tmp(nswap))
     allocate(ds_tmp(nswap),nt_tmp(nswap))
     
     write(6,*) "shift indices", istart,iend,nswap
     !! Temporary store of the final nswap elements
     x_tmp(1:nswap) = x(iend-nswap+1:iend)
     y_tmp(1:nswap) = y(iend-nswap+1:iend)
     xn_tmp(1:nswap) = xn(iend-nswap+1:iend)
     yn_tmp(1:nswap) = yn(iend-nswap+1:iend)
     ds_tmp(1:nswap) = ds(iend-nswap+1:iend)
     nt_tmp(1:nswap) = nt(iend-nswap+1:iend)       
     
     
     !! Shift
     do i=1,shift_size
        !! original index
        i_old = iend - nswap + 1 - i
        !! new index
        i_new = i_old + nswap
        
        x(i_new) = x(i_old)
        y(i_new) = y(i_old)
        xn(i_new) = xn(i_old)
        yn(i_new) = yn(i_old)
        ds(i_new) = ds(i_old)
        nt(i_new) = nt(i_old)                                                
     end do
     
     !! Copy temp back to start of band
     x(istart:istart+nswap-1) = x_tmp(1:nswap)
     y(istart:istart+nswap-1) = y_tmp(1:nswap)
     xn(istart:istart+nswap-1) = xn_tmp(1:nswap)
     yn(istart:istart+nswap-1) = yn_tmp(1:nswap)
     ds(istart:istart+nswap-1) = ds_tmp(1:nswap)
     nt(istart:istart+nswap-1) = nt_tmp(1:nswap)                         
     
     deallocate(x_tmp,y_tmp,xn_tmp,yn_tmp,ds_tmp,nt_tmp)
     
     return
  end subroutine shift_indices
!! ------------------------------------------------------------------------------------------------ 
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
recursive subroutine quicksort(a, first, last)
  implicit none
  double precision  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     call swap_nodes(i,j)
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksort(a, first, i-1)
  if (j+1 < last)  call quicksort(a, j+1, last)
end subroutine quicksort
!! ------------------------------------------------------------------------------------------------
  subroutine swap_nodes(i,j)
     integer :: i,j
     double precision :: tmp
     integer :: itmp

     !! x is already swapped by sub-routine quicksort     
     tmp = y(j);y(j)=y(i);y(i)=tmp
     tmp = xn(j);xn(j)=xn(i);xn(i)=tmp
     tmp = yn(j);yn(j)=yn(i);yn(i)=tmp
     tmp = ds(j);ds(j)=ds(i);ds(i)=tmp
     itmp = nt(j);nt(j)=nt(i);nt(i)=itmp                    
     return
  end subroutine swap_nodes
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
recursive subroutine quicksorty(a, first, last)
  implicit none
  double precision  a(*), x, t
  integer first, last
  integer i, j

  x = a( (first+last) / 2 )
  i = first
  j = last
  do
     do while (a(i) < x)
        i=i+1
     end do
     do while (x < a(j))
        j=j-1
     end do
     if (i >= j) exit
     t = a(i);  a(i) = a(j);  a(j) = t
     call swap_nodesy(i,j)
     i=i+1
     j=j-1
  end do
  if (first < i-1) call quicksorty(a, first, i-1)
  if (j+1 < last)  call quicksorty(a, j+1, last)
end subroutine quicksorty
!! ------------------------------------------------------------------------------------------------
  subroutine swap_nodesy(i,j)
     integer :: i,j
     double precision :: tmp
     integer :: itmp

     !! y is already swapped by sub-routine quicksort     
     tmp = x(j);x(j)=x(i);x(i)=tmp
     tmp = xn(j);xn(j)=xn(i);xn(i)=tmp
     tmp = yn(j);yn(j)=yn(i);yn(i)=tmp
     tmp = ds(j);ds(j)=ds(i);ds(i)=tmp
     itmp = nt(j);nt(j)=nt(i);nt(i)=itmp                        
     return
  end subroutine swap_nodesy
!! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
end module inputoutput
