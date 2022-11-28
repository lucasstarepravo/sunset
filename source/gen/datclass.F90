program datgen
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! DATGEN program to generate node sets for input
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  use kind_parameters
  use common_parameter
  use global_variables 
  implicit none

  !! xbcond/ybcond = 0,1,2 for none, periodic or symmetric BCs respectively.
  !! btype = 0,1,2,3 for wall, inflow, outflow or periodic/symmetric respectively

  real(rkind) :: x,y

  integer ipart,itest
  integer i,j,icha,nn,ve_model,ii,jj
  double precision h0,r_mag,yl,D_cyl,S_cyl
  double precision areafluid,ns,varresratio
  integer xbcond,ybcond
  
  double precision :: a0,a1,a2,a3,a4,a5 !! NACA coefficients
  double precision :: dydx,temp,tmp2,tmp,xtec,rtec,atec,thtec
  real(rkind),dimension(2) :: tmpN,rn
  integer(ikind) :: iround,nround,nsearch,ipartrow,ipartrowm1,nnear,npdps,idown,iup
  logical :: keepgoing,keepgoing2,skip
  real(rkind),dimension(:,:),allocatable :: pdp
  real(rkind) :: dxmin,dxmax,dist2bound,dxtmp,minpdp
  real(rkind),dimension(:),allocatable :: tmpX
  real(rkind),dimension(:),allocatable :: pdp_x,pdp_y,pdp_dist2
  logical,dimension(:),allocatable :: pdp_active
  integer(ikind),dimension(:),allocatable :: pdp_freeindices
  real(rkind),dimension(5) :: pdp_x_new,pdp_y_new
  real(rkind) :: thup,thdown,th_increment,sumdG,sumG
  integer(ikind) :: npdps_new,inew,block_left,block_right,block_new,block_delete,block_end
  real(rkind) :: b0,b1,b2,b3 !! Coefficients in dx functions

  write(*,*) 'Cases: '
  write(*,*) '  case 1:  box for Rayleigh-Taylor'
  write(*,*) '  case 2:  unit torus with 4 cylinders'
  write(*,*) '  case 3:  unit torus'
  write(*,*) '  case 4:  porous blobby'
  write(*,*) '  case 5:  NACA 0012'  
  write(*,*) '  case 6:  Channel flows with obstacle'    
  write(*,*) '  case 7:  Something periodic with obstacle'      
  write(*,*) '  '
  write(*,*) 'Input test case number: '
  read(*,*) itest


  select case (itest) 
!! ------------------------------------------------------------------------------------------------
  case(1)   
      !! EMPTY
!! ------------------------------------------------------------------------------------------------
  case(2) !! GRID!
     yl=0.0125d0!0.0125d0  ! channel width
     xl=1.0d0 ! channel length
     dx0=xl/1000.0       !15
     xbcond=0;ybcond=1     
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 1/)  
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 0
     call make_boundary_edge_vectors
     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     !! Uniform resolution
     dxmax = dx0  
     dxmin = dx0
     dxb=dx0
     dx_in=dx0;dx_out=dx0
     call make_boundary_particles
     call make_boundary_blobs               
     ipart = nb   
     
     dx = dx0
     x = xb_min + 5.0*dx
     do while(x.lt.xb_max-4.5*dx)
        y = yb_min + 0.5*dx
        do while(y.lt.yb_max)
              ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;dxp(ipart) = dx
           y = y + dx
        enddo
        x = x + dx
     enddo     
                        
     npfb = ipart
     dx0 = maxval(dxp(1:npfb))
     

     write(*,*) 'nb,npfb= ', nb,npfb,nbio
       

!! ------------------------------------------------------------------------------------------------
  case(3) !! Kolmogorov flow

     yl=2.0d0*pi;h0=yl
     xl=1.0d0*yl
     dx0=yl/64.0!0.025d0
     xbcond=1;ybcond=1     

     !   JRCK boundary conditions
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 3, 3, 3/)  
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)

     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     !! outsource making wall particles... JRCK
     dxb=dx0;dx_in=dx0;dx_out=dx0
     call make_boundary_particles
     ipart = nb   

      ! -- fluid particles-- 
     dx = dx0
     x = xb_min + 0.5*dx
     do while(x.lt.xb_max)
        y = yb_min + 0.5*dx
        do while(y.lt.yb_max)
              ipart=ipart+1;xp(ipart)= x;yp(ipart)= y;dxp(ipart) = dx
           y = y + dx
        enddo
        x = x + dx
     enddo
    
     npfb = ipart
     write(*,*) 'nb,npfb= ', nb,npfb

!! ------------------------------------------------------------------------------------------------
case(4) !! A sort of porous media... for porous Rayleigh-Taylor stuff

     D_cyl = 1.0d0
     S_cyl = D_cyl*1.2d0
     h0=D_cyl/2.0d0      !cylinder radius
     yl=S_cyl*sqrt(3.0d0)*2.0d0  ! height
     xl=2.0*S_cyl ! width
     dx0=h0/30.0       !15
     xbcond=1;ybcond=1     
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 3, 3, 3/)  
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 13
     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,6),blob_rotation(nb_blobs),blob_ellipse(nb_blobs))
     b0=2.5d0*h0;b1=b0*sqrt(3.0d0)/2.0d0;b2=b0/2.0d0
     blob_centre(1,:)=(/0.d0,0.d0/); !! Central row
     blob_centre(2,:)=(/-S_cyl,0.d0/);      
     blob_centre(3,:)=(/S_cyl,0.d0/); 
     blob_centre(4,:)=(/0.d0,S_cyl*sqrt(3.0d0)/); !! Upper row
     blob_centre(5,:)=(/-S_cyl,S_cyl*sqrt(3.0d0)/);                
     blob_centre(6,:)=(/S_cyl,S_cyl*sqrt(3.0d0)/);                  
     blob_centre(7,:)=(/0.d0,-S_cyl*sqrt(3.0d0)/); !! Lower row
     blob_centre(8,:)=(/-S_cyl,-S_cyl*sqrt(3.0d0)/);                
     blob_centre(9,:)=(/S_cyl,-S_cyl*sqrt(3.0d0)/);                  

     blob_centre(10,:) = (/-0.5d0*S_cyl,S_cyl*sqrt(3.0d0)/2.0d0/) !! Row 1/2
     blob_centre(11,:) = (/0.5d0*S_cyl,S_cyl*sqrt(3.0d0)/2.0d0/) 
     blob_centre(12,:) = (/-0.5d0*S_cyl,-S_cyl*sqrt(3.0d0)/2.0d0/) !! Row -1/2
     blob_centre(13,:) = (/0.5d0*S_cyl,-S_cyl*sqrt(3.0d0)/2.0d0/) 



     do i=1,nb_blobs
        blob_coeffs(i,:)=h0*(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0/);blob_rotation(i)=-pi/9.0d0;blob_ellipse(i)=1
     end do

     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     varresratio = 3.0d0  !! Ratio for scaling near the solid objects
     dxmax = dx0  
     dxmin = dx0/varresratio
     dxb=dx0/varresratio;dx_in=4.0d0*dxmax;dx_out=dx_in !! dx for solids and in/outs...!! Ratio for scaling far field...
     call make_boundary_particles
     call make_boundary_blobs               
     ipart = nb   
     
     b0 = 3.0d0;b1 = 40.0d0;b2 = 30.0d0
         
     !! Initialise a line of potential dot points   
     nsearch = ceiling(yb_max-yb_min)/dxmin/2.0d0
     allocate(pdp_x(5*nsearch),pdp_y(5*nsearch))
     npdps = nsearch
     y=yb_min-0.5d0*dxmin
     i=0
     do while (y.lt.yb_max)
        y = y + dxmin/2.0d0;i=i+1
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_x(i) = xb_min + temp
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_y(i) = y 
     end do          
     npdps=i         
     allocate(pdp_dist2(5*nsearch));pdp_dist2=0.0d0   
        
     minpdp = minval(pdp_x(1:npdps))
     do while (minpdp.le.xb_max)  !! Keep going until all PDPs have passed out the end of the domain
        
        !! Pick the left-most pdp
        j=minloc(pdp_x(1:npdps),DIM=1)
        keepgoing = .true.
        x=pdp_x(j);y=pdp_y(j)
        pdp_dist2(1:npdps) = (x-pdp_x(1:npdps))**2.0d0 + (y-pdp_y(1:npdps))**2.0d0
        
        !! How far are we from the boundaries?
        dist2bound = xb_max-xb_min + 100.0d0
        i=1
        sumdG=0.0d0;sumG=0.0d0
        do while(i.le.nb)
           tmpN(1) = x - xp(i);tmpN(2) = y - yp(i);temp = sqrt(dot_product(tmpN,tmpN))
           if(i.gt.nbio.and.temp.le.dist2bound) dist2bound = temp
           i=i+1
           if(temp.le.4.25*dxp(i-1)) keepgoing = .false. !! Too close to bound, leave it.           
        end do        
        
        !! And what is the spacing, based on dist2bound?
        if(dist2bound.le.b0*dx0) then  !! Close - set to dxmin
           dx = dxmin                     
        else if(dist2bound.le.b1*dx0)then  !! A bit further out, smoothly vary from dxmin to dx0
           dx = 0.5d0*(dx0+dxmin) - 0.5d0*(dx0-dxmin)*cos((dist2bound-b0*dx0)*pi/((b1-b0)*dx0))  
        else if(dist2bound.le.b1*dx0+b2*dx_in)then  !! Further still: linearly vary from dx0 to dx_in
           dx = dx0 + (dx_in-dx0)*((dist2bound-b1*dx0)/(b2*dx_in))
        else     !! Far out: set to dx_in
           dx = dx_in
        end if       
      
               
        !! Check whether to place a particle here, based on some criteria
        !! --------------------------------------------------------------
           !! Are we within the object!!?!?!
           do i=1,nb_blobs
              temp = sqrt((x-blob_centre(i,1))**2. + (y-blob_centre(i,2))**2.)
              if(x-blob_centre(i,1).ge.0.0d0)then
                 tmp2 = asin((y-blob_centre(i,2))/temp)
              else
                 tmp2 = pi-asin((y-blob_centre(i,2))/temp)
              endif              
              tmp2 = tmp2 - blob_rotation(i)
              if(blob_ellipse(i).eq.0)then !! blob is blob
                 r_mag = blob_coeffs(i,1) + blob_coeffs(i,2)*sin(tmp2) + blob_coeffs(i,3)*sin(2.0*tmp2) + &
                      blob_coeffs(i,4)*sin(3.0*tmp2) + blob_coeffs(i,5)*sin(4.0*tmp2) + blob_coeffs(i,6)*sin(5.0*tmp2)
              else              !! blob is ellipse
                 r_mag = blob_coeffs(i,1)*blob_coeffs(i,2) &
                 /sqrt(blob_coeffs(i,2)*blob_coeffs(i,2)*cos(tmp2)**2. + &
                   blob_coeffs(i,1)*blob_coeffs(i,1)*sin(tmp2)**2.)                       
              end if                      
              if(temp.le.r_mag)then
                 keepgoing = .false.
              end if
           end do
           
           !! Are we too close to a boundary?
           if(x-0.5d0*dx.le.xb_min.or.x+0.5d0*dx.ge.xb_max.or.y-0.5d0*dx.le.yb_min.or.y+0.5d0*dx.ge.yb_max)then
              keepgoing = .false.
           end if            
                                  

        !! END CRITERIA
        !! --------------------------------------------------------------

        !! Place a particle here
        if(keepgoing) then
           ipart = ipart + 1
           call random_number(temp);temp = temp -0.5d0;xp(ipart) = pdp_x(j) + temp*dxmin*0.5d0
           call random_number(temp);temp = temp -0.5d0;yp(ipart) = pdp_y(j) + temp*dxmin*0.5d0
           dxp(ipart) = dx     
        end if           
        
        !! Deactive all pdps within dx of this pdp
        !! Search down
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i-1;
           if(i.eq.0) then
              idown = 0
              thdown = -0.49999999d0*pi
              keepgoing = .false.
           else if(pdp_dist2(i).ge.dx*dx) then             
              idown = i
              thdown = atan2((pdp_y(idown)-y),(pdp_x(idown)-x))
              keepgoing = .false.    
 
           end if       
        end do
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i+1
           if(i.eq.npdps+1) then
              iup = npdps+1
              thup = 0.4999999999d0*pi
              keepgoing = .false.              
           else if(pdp_dist2(i).ge.dx*dx) then 
              iup = i
              thup = atan2((pdp_y(iup)-y),(pdp_x(iup)-x))              
              keepgoing = .false.
           
           end if
        end do     
        
        !! Temporary store for the new pdps
        th_increment = (thup-thdown)/5.0d0
        inew = 0
        do i=1,5
           temp = y + dx*sin(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)    
           if(temp.ge.yb_min.and.temp.le.yb_max) then
              inew = inew + 1
              pdp_x_new(inew) = x + dx*cos(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)
              pdp_y_new(inew) = temp
           end if            
        end do
 

        block_new = idown + 1
        block_right = block_new + inew
        block_delete = iup - idown - 1
        npdps_new = npdps + inew - block_delete

        !! Shunt indices above pdp of interest
        if(block_delete.gt.inew) then !! Shunt to LEFT
           do i=block_right,npdps_new
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do
        end if
        if(block_delete.lt.inew) then !! Shunt to RIGHT
           do i=npdps_new,block_right,-1
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do

        end if
                  
        !! Insert any new pdps
        if(inew.ne.0) then          
           pdp_x(block_new:block_new+inew-1) = pdp_x_new(1:inew)
           pdp_y(block_new:block_new+inew-1) = pdp_y_new(1:inew)        
        end if       

        npdps = npdps_new
                  
                
        !! How left is the left-most PDP?                 
        minpdp = minval(pdp_x(1:npdps))
                
        write(6,*) ipart,(minval(pdp_x(1:npdps))-xb_min)/(xb_max-xb_min)
     end do  
                         
                                                  
                        
     npfb = ipart
     dx0 = maxval(dxp(1:npfb))
     

     write(*,*) 'nb,npfb= ', nb,npfb,nbio



!! ------------------------------------------------------------------------------------------------
case(5) !! Inflow/outflow tube for simple flames

     yl=0.0125d0!0.0125d0  ! channel width
     xl=1.0d0 ! channel length
     dx0=xl/1000.0       !15
     xbcond=0;ybcond=1     
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 1/)  
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 0
!     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,6),blob_rotation(nb_blobs),blob_ellipse(nb_blobs))
!     b0=2.5d0*h0;b1=b0*sqrt(3.0d0)/2.0d0;b2=b0/2.0d0
!     blob_centre(1,:)=(/0.d0,0.d0/); !! Central
!     do i=1,nb_blobs
!        blob_coeffs(i,:)=h0*(/1.0d0,0.4d0,0.0d0,0.0d0,0.0d0,0.0d0/);blob_rotation(i)=-pi/9.0d0;blob_ellipse(i)=1
!     end do

     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     varresratio = 1.0d0  !! Ratio for scaling near the solid objects
     dxmax = dx0  
     dxmin = dx0/varresratio
     dxb=dx0/varresratio;dx_in=1.0d0*dxmax;dx_out=dx0*1.0d0  !! dx for solids and in/outs..
     call make_boundary_particles
     call make_boundary_blobs               
     ipart = nb   
     
     !! Wall to first growth; region of first growth; region of 2nd growth.
     b0 = 3.0d0;b1 = 40.0d0;b2 = 50.0d0
         
     !! Initialise a line of potential dot points   
     nsearch = ceiling(yb_max-yb_min)/dxmin/2.0d0
     allocate(pdp_x(5*nsearch),pdp_y(5*nsearch))
     npdps = nsearch
     y=yb_min-0.5d0*dxmin
     i=0
     do while (y.lt.yb_max)
        y = y + dxmin/2.0d0;i=i+1
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_x(i) = xb_min + temp
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_y(i) = y 
     end do          
     npdps=i         
     allocate(pdp_dist2(5*nsearch));pdp_dist2=0.0d0   
        
     minpdp = minval(pdp_x(1:npdps))
     do while (minpdp.le.xb_max)  !! Keep going until all PDPs have passed out the end of the domain
        
        !! Pick the left-most pdp
        j=minloc(pdp_x(1:npdps),DIM=1)
        keepgoing = .true.
        x=pdp_x(j);y=pdp_y(j)
        pdp_dist2(1:npdps) = (x-pdp_x(1:npdps))**2.0d0 + (y-pdp_y(1:npdps))**2.0d0
        
        !! How far are we from the boundaries?
        dist2bound = xb_max-xb_min + 100.0d0
        i=1
        sumdG=0.0d0;sumG=0.0d0
        do while(i.le.nb)
           tmpN(1) = x - xp(i);tmpN(2) = y - yp(i);temp = sqrt(dot_product(tmpN,tmpN))
           if(i.gt.nbio.and.temp.le.dist2bound) dist2bound = temp
           i=i+1
           if(temp.le.4.25*dxp(i-1)) keepgoing = .false. !! Too close to bound, leave it.           
        end do        
        
        !! And what is the spacing, based on dist2bound?
!! Over-ride object tests
if(abs(x).le.0.075) then
   dist2bound=0.0d0
else
   dist2bound=min(abs(x-0.075),abs(x+0.075))
!   dist2bound=abs(x+0.4)
endif   
   
        if(dist2bound.le.b0*dx0) then  !! Close - set to dxmin
           dx = dxmin                     
        else if(dist2bound.le.b1*dx0)then  !! A bit further out, smoothly vary from dxmin to dx0
           dx = 0.5d0*(dx0+dxmin) - 0.5d0*(dx0-dxmin)*cos((dist2bound-b0*dx0)*pi/((b1-b0)*dx0))  
        else if(dist2bound.le.b1*dx0+b2*dx_out)then  !! Further still: linearly vary from dx0 to dx_out
           dx = dx0 + (dx_out-dx0)*((dist2bound-b1*dx0)/(b2*dx_out))
        else     !! Far out: set to dx_out
           dx = dx_out
        end if               
      
               
        !! Check whether to place a particle here, based on some criteria
        !! --------------------------------------------------------------
           !! Are we within the object!!?!?!
           do i=1,nb_blobs
              temp = sqrt((x-blob_centre(i,1))**2. + (y-blob_centre(i,2))**2.)
              if(x-blob_centre(i,1).ge.0.0d0)then
                 tmp2 = asin((y-blob_centre(i,2))/temp)
              else
                 tmp2 = pi-asin((y-blob_centre(i,2))/temp)
              endif              
              tmp2 = tmp2 - blob_rotation(i)
              if(blob_ellipse(i).eq.0)then !! blob is blob
                 r_mag = blob_coeffs(i,1) + blob_coeffs(i,2)*sin(tmp2) + blob_coeffs(i,3)*sin(2.0*tmp2) + &
                      blob_coeffs(i,4)*sin(3.0*tmp2) + blob_coeffs(i,5)*sin(4.0*tmp2) + blob_coeffs(i,6)*sin(5.0*tmp2)
              else              !! blob is ellipse
                 r_mag = blob_coeffs(i,1)*blob_coeffs(i,2) &
                 /sqrt(blob_coeffs(i,2)*blob_coeffs(i,2)*cos(tmp2)**2. + &
                   blob_coeffs(i,1)*blob_coeffs(i,1)*sin(tmp2)**2.)                       
              end if                      
              if(temp.le.r_mag)then
                 keepgoing = .false.
              end if
           end do
           
           !! Are we too close to a boundary?
           if(x-0.5d0*dx.le.xb_min.or.x+0.5d0*dx.ge.xb_max.or.y-0.5d0*dx.le.yb_min.or.y+0.5d0*dx.ge.yb_max)then
              keepgoing = .false.
           end if            
                                  

        !! END CRITERIA
        !! --------------------------------------------------------------

        !! Place a particle here
        if(keepgoing) then
           ipart = ipart + 1
           call random_number(temp);temp = temp -0.5d0;xp(ipart) = pdp_x(j) + temp*dxmin*0.5d0
           call random_number(temp);temp = temp -0.5d0;yp(ipart) = pdp_y(j) + temp*dxmin*0.5d0
           dxp(ipart) = dx     
        end if           
        
        !! Deactive all pdps within dx of this pdp
        !! Search down
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i-1;
           if(i.eq.0) then
              idown = 0
              thdown = -0.49999999d0*pi
              keepgoing = .false.
           else if(pdp_dist2(i).ge.dx*dx) then             
              idown = i
              thdown = atan2((pdp_y(idown)-y),(pdp_x(idown)-x))
              keepgoing = .false.    
 
           end if       
        end do
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i+1
           if(i.eq.npdps+1) then
              iup = npdps+1
              thup = 0.4999999999d0*pi
              keepgoing = .false.              
           else if(pdp_dist2(i).ge.dx*dx) then 
              iup = i
              thup = atan2((pdp_y(iup)-y),(pdp_x(iup)-x))              
              keepgoing = .false.
           
           end if
        end do     
        
        !! Temporary store for the new pdps
        th_increment = (thup-thdown)/5.0d0
        inew = 0
        do i=1,5
           temp = y + dx*sin(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)    
           if(temp.ge.yb_min.and.temp.le.yb_max) then
              inew = inew + 1
              pdp_x_new(inew) = x + dx*cos(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)
              pdp_y_new(inew) = temp
           end if            
        end do
 

        block_new = idown + 1
        block_right = block_new + inew
        block_delete = iup - idown - 1
        npdps_new = npdps + inew - block_delete

        !! Shunt indices above pdp of interest
        if(block_delete.gt.inew) then !! Shunt to LEFT
           do i=block_right,npdps_new
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do
        end if
        if(block_delete.lt.inew) then !! Shunt to RIGHT
           do i=npdps_new,block_right,-1
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do

        end if
                  
        !! Insert any new pdps
        if(inew.ne.0) then          
           pdp_x(block_new:block_new+inew-1) = pdp_x_new(1:inew)
           pdp_y(block_new:block_new+inew-1) = pdp_y_new(1:inew)        
        end if       

        npdps = npdps_new
                  
                
        !! How left is the left-most PDP?                 
        minpdp = minval(pdp_x(1:npdps))
                
        write(6,*) ipart,(minval(pdp_x(1:npdps))-xb_min)/(xb_max-xb_min)
     end do  
                         
                                                  
                        
     npfb = ipart
     dx0 = maxval(dxp(1:npfb))

     write(*,*) 'nb,npfb= ', nb,npfb,nbio
     
!! ------------------------------------------------------------------------------------------------
case(6) !! Channel flows, propagating front

     xl=1.0d0 ! channel length
     h0=xl/20.0d0   !cylinder radius
     yl=8.0d0*h0  ! channel width
     dx0=h0/20.0       !15
     xbcond=0;ybcond=2     
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 2, 3, 1/)  
     b_node(1,:) = (/ -0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/ 0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/ 0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/ -0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 1
     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,6),blob_rotation(nb_blobs),blob_ellipse(nb_blobs))
     b0=2.5d0*h0;b1=b0*sqrt(3.0d0)/2.0d0;b2=b0/2.0d0
     blob_centre(1,:)=(/-0.0d0*h0,0.d0/); !! Central
     do i=1,nb_blobs
        blob_coeffs(i,:)=h0*(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0/);blob_rotation(i)=-pi/9.0d0;blob_ellipse(i)=1
     end do

     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     varresratio = 2.0d0  !! Ratio for scaling near the solid objects
     dxmax = dx0  
     dxmin = dx0/varresratio
     dxb=dx0/varresratio;dx_in=2.0d0*dxmax;dx_out=2.0d0*dxmax  !! dx for solids and in/outs...!! 
     call make_boundary_particles
     call make_boundary_blobs               
     ipart = nb   
     
     !! Wall to first growth; region of first growth; region of 2nd growth.
     b0 = 3.0d0;b1 = 40.0d0;b2 = 50.0d0
         
     !! Initialise a line of potential dot points   
     nsearch = ceiling(yb_max-yb_min)/dxmin/2.0d0
     allocate(pdp_x(5*nsearch),pdp_y(5*nsearch))
     npdps = nsearch
     y=yb_min-0.5d0*dxmin
     i=0
     do while (y.lt.yb_max)
        y = y + dxmin/2.0d0;i=i+1
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_x(i) = xb_min + temp
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_y(i) = y 
     end do          
     npdps=i         
     allocate(pdp_dist2(5*nsearch));pdp_dist2=0.0d0   
        
     minpdp = minval(pdp_x(1:npdps))
     do while (minpdp.le.xb_max)  !! Keep going until all PDPs have passed out the end of the domain
        
        !! Pick the left-most pdp
        j=minloc(pdp_x(1:npdps),DIM=1)
        keepgoing = .true.
        x=pdp_x(j);y=pdp_y(j)
        pdp_dist2(1:npdps) = (x-pdp_x(1:npdps))**2.0d0 + (y-pdp_y(1:npdps))**2.0d0
        
        !! How far are we from the boundaries?
        dist2bound = xb_max-xb_min + 100.0d0
        i=1
        sumdG=0.0d0;sumG=0.0d0
        do while(i.le.nb)
           tmpN(1) = x - xp(i);tmpN(2) = y - yp(i);temp = sqrt(dot_product(tmpN,tmpN))
           if(i.gt.nbio.and.temp.le.dist2bound) dist2bound = temp
           i=i+1
           if(temp.le.4.25*dxp(i-1)) keepgoing = .false. !! Too close to bound, leave it.           
        end do     
     
!! Stretch high-res region downstream of flameholder        
if((x-blob_centre(1,1)).gt.0.0d0) then
 temp = max(exp(-(x-blob_centre(1,1))/h0),1.0d0)
 dist2bound = dist2bound*temp
 dxio = dx_out
else
 dxio = dx_in 
endif            
        
        !! And what is the spacing, based on dist2bound?
        if(dist2bound.le.b0*dx0) then  !! Close - set to dxmin
           dx = dxmin                     
        else if(dist2bound.le.b1*dx0)then  !! A bit further out, smoothly vary from dxmin to dx0
           dx = 0.5d0*(dx0+dxmin) - 0.5d0*(dx0-dxmin)*cos((dist2bound-b0*dx0)*pi/((b1-b0)*dx0))  
        else if(dist2bound.le.b1*dx0+b2*dxio)then  !! Further still: linearly vary from dx0 to dxio
           dx = dx0 + (dxio-dx0)*((dist2bound-b1*dx0)/(b2*dxio))
        else     !! Far out: set to dxio
           dx = dxio
        end if       
      
               
        !! Check whether to place a particle here, based on some criteria
        !! --------------------------------------------------------------
           !! Are we within the object!!?!?!
           do i=1,nb_blobs
              temp = sqrt((x-blob_centre(i,1))**2. + (y-blob_centre(i,2))**2.)
              if(x-blob_centre(i,1).ge.0.0d0)then
                 tmp2 = asin((y-blob_centre(i,2))/temp)
              else
                 tmp2 = pi-asin((y-blob_centre(i,2))/temp)
              endif              
              tmp2 = tmp2 - blob_rotation(i)
              if(blob_ellipse(i).eq.0)then !! blob is blob
                 r_mag = blob_coeffs(i,1) + blob_coeffs(i,2)*sin(tmp2) + blob_coeffs(i,3)*sin(2.0*tmp2) + &
                      blob_coeffs(i,4)*sin(3.0*tmp2) + blob_coeffs(i,5)*sin(4.0*tmp2) + blob_coeffs(i,6)*sin(5.0*tmp2)
              else              !! blob is ellipse
                 r_mag = blob_coeffs(i,1)*blob_coeffs(i,2) &
                 /sqrt(blob_coeffs(i,2)*blob_coeffs(i,2)*cos(tmp2)**2. + &
                   blob_coeffs(i,1)*blob_coeffs(i,1)*sin(tmp2)**2.)                       
              end if                      
              if(temp.le.r_mag)then
                 keepgoing = .false.
              end if
           end do
           
           !! Are we too close to a boundary?
           if(x-0.5d0*dx.le.xb_min.or.x+0.5d0*dx.ge.xb_max.or.y-0.5d0*dx.le.yb_min.or.y+0.5d0*dx.ge.yb_max)then
              keepgoing = .false.
           end if            
                                  

        !! END CRITERIA
        !! --------------------------------------------------------------

        !! Place a particle here
        if(keepgoing) then
           ipart = ipart + 1
           call random_number(temp);temp = temp -0.5d0;xp(ipart) = pdp_x(j) + temp*dxmin*0.5d0
           call random_number(temp);temp = temp -0.5d0;yp(ipart) = pdp_y(j) + temp*dxmin*0.5d0
           dxp(ipart) = dx     
        end if           
        
        !! Deactive all pdps within dx of this pdp
        !! Search down
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i-1;
           if(i.eq.0) then
              idown = 0
              thdown = -0.49999999d0*pi
              keepgoing = .false.
           else if(pdp_dist2(i).ge.dx*dx) then             
              idown = i
              thdown = atan2((pdp_y(idown)-y),(pdp_x(idown)-x))
              keepgoing = .false.    
 
           end if       
        end do
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i+1
           if(i.eq.npdps+1) then
              iup = npdps+1
              thup = 0.4999999999d0*pi
              keepgoing = .false.              
           else if(pdp_dist2(i).ge.dx*dx) then 
              iup = i
              thup = atan2((pdp_y(iup)-y),(pdp_x(iup)-x))              
              keepgoing = .false.
           
           end if
        end do     
        
        !! Temporary store for the new pdps
        th_increment = (thup-thdown)/5.0d0
        inew = 0
        do i=1,5
           temp = y + dx*sin(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)    
           if(temp.ge.yb_min.and.temp.le.yb_max) then
              inew = inew + 1
              pdp_x_new(inew) = x + dx*cos(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)
              pdp_y_new(inew) = temp
           end if            
        end do
 

        block_new = idown + 1
        block_right = block_new + inew
        block_delete = iup - idown - 1
        npdps_new = npdps + inew - block_delete

        !! Shunt indices above pdp of interest
        if(block_delete.gt.inew) then !! Shunt to LEFT
           do i=block_right,npdps_new
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do
        end if
        if(block_delete.lt.inew) then !! Shunt to RIGHT
           do i=npdps_new,block_right,-1
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do

        end if
                  
        !! Insert any new pdps
        if(inew.ne.0) then          
           pdp_x(block_new:block_new+inew-1) = pdp_x_new(1:inew)
           pdp_y(block_new:block_new+inew-1) = pdp_y_new(1:inew)        
        end if       

        npdps = npdps_new
                  
                
        !! How left is the left-most PDP?                 
        minpdp = minval(pdp_x(1:npdps))
                
        write(6,*) ipart,(minval(pdp_x(1:npdps))-xb_min)/(xb_max-xb_min)
     end do  
                         
                                                  
                        
     npfb = ipart
     dx0 = maxval(dxp(1:npfb))
     
     write(*,*) 'nb,npfb= ', nb,npfb,nbio

!! ------------------------------------------------------------------------------------------------
case(7) !! Something periodic

     D_cyl = 1.0d0
     S_cyl = D_cyl*1.2d0
     h0=D_cyl/2.0d0      !cylinder radius
     yl=2.0d0*S_cyl ! box height
     xl=sqrt(3.0d0)*S_cyl ! channel length
     dx0=D_cyl/50.0       !75
     xbcond=1;ybcond=1     
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 3, 3, 3/)  
     b_node(1,:) = (/-0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/-0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 7
     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,6),blob_rotation(nb_blobs),blob_ellipse(nb_blobs))
     blob_centre(1,:) = (/0.0d0,0.0d0/)   !! Row 0
     blob_centre(2,:) = (/0.0d0,-S_cyl/)
     blob_centre(3,:) = (/0.0d0,S_cyl/)
     blob_centre(4,:) = (/S_cyl*sqrt(3.0d0)/2.0d0,-0.5d0*S_cyl/) !! Row 1/2
     blob_centre(5,:) = (/S_cyl*sqrt(3.0d0)/2.0d0,0.5d0*S_cyl/)
     blob_centre(6,:) = (/-S_cyl*sqrt(3.0d0)/2.0d0,-0.5d0*S_cyl/) !! Row -1/2
     blob_centre(7,:) = (/-S_cyl*sqrt(3.0d0)/2.0d0,0.5d0*S_cyl/)
                          
     
     do i=1,nb_blobs
        blob_coeffs(i,:)=h0*(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0/);blob_rotation(i)=-pi/9.0d0;blob_ellipse(i)=1
     end do

     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     varresratio = 3.0d0  !! Ratio for scaling near the solid objects
     dxmax = dx0  
     dxmin = dx0/varresratio
     dxb=dx0/varresratio;dx_in=2.0d0*dxmax;dx_out=dx_in  !! dx for solids and in/outs...!! Ratio for scaling far field...
     call make_boundary_particles
     call make_boundary_blobs               
     ipart = nb   
     
     b0 = 2.0d0;b1 = 40.0d0;b2 = 30.0d0
         
     !! Initialise a line of potential dot points   
     nsearch = ceiling(yb_max-yb_min)/dxmin/2.0d0
     allocate(pdp_x(10*nsearch),pdp_y(10*nsearch))
     npdps = nsearch
     y=yb_min-0.5d0*dxmin
     i=0
     do while (y.lt.yb_max)
        y = y + dxmin/2.0d0;i=i+1
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_x(i) = xb_min + temp
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_y(i) = y 
     end do          
     npdps=i         
     allocate(pdp_dist2(10*nsearch));pdp_dist2=0.0d0   
        
     minpdp = minval(pdp_x(1:npdps))
     do while (minpdp.le.xb_max)  !! Keep going until all PDPs have passed out the end of the domain
        
        !! Pick the left-most pdp
        j=minloc(pdp_x(1:npdps),DIM=1)
        keepgoing = .true.
        x=pdp_x(j);y=pdp_y(j)
        pdp_dist2(1:npdps) = (x-pdp_x(1:npdps))**2.0d0 + (y-pdp_y(1:npdps))**2.0d0
        
        !! How far are we from the boundaries?
        dist2bound = xb_max-xb_min + 100.0d0
        i=1
        sumdG=0.0d0;sumG=0.0d0
        do while(i.le.nb)
           tmpN(1) = x - xp(i);tmpN(2) = y - yp(i);temp = sqrt(dot_product(tmpN,tmpN))
           if(i.gt.nbio.and.temp.le.dist2bound) dist2bound = temp
           i=i+1
           if(temp.le.4.25*dxp(i-1)) keepgoing = .false. !! Too close to bound, leave it.           
        end do        
        
        !! And what is the spacing, based on dist2bound?
        if(dist2bound.le.b0*dx0) then  !! Close - set to dxmin
           dx = dxmin                     
        else if(dist2bound.le.b1*dx0)then  !! A bit further out, smoothly vary from dxmin to dx0
           dx = 0.5d0*(dx0+dxmin) - 0.5d0*(dx0-dxmin)*cos((dist2bound-b0*dx0)*pi/((b1-b0)*dx0))  
        else if(dist2bound.le.b1*dx0+b2*dx_in)then  !! Further still: linearly vary from dx0 to dx_in
           dx = dx0 + (dx_in-dx0)*((dist2bound-b1*dx0)/(b2*dx_in))
        else     !! Far out: set to dx_in
           dx = dx_in
        end if       
      
               
        !! Check whether to place a particle here, based on some criteria
        !! --------------------------------------------------------------
           !! Are we within the object!!?!?!
           do i=1,nb_blobs
              temp = sqrt((x-blob_centre(i,1))**2. + (y-blob_centre(i,2))**2.)
              if(x-blob_centre(i,1).ge.0.0d0)then
                 tmp2 = asin((y-blob_centre(i,2))/temp)
              else
                 tmp2 = pi-asin((y-blob_centre(i,2))/temp)
              endif              
              tmp2 = tmp2 - blob_rotation(i)
              if(blob_ellipse(i).eq.0)then !! blob is blob
                 r_mag = blob_coeffs(i,1) + blob_coeffs(i,2)*sin(tmp2) + blob_coeffs(i,3)*sin(2.0*tmp2) + &
                      blob_coeffs(i,4)*sin(3.0*tmp2) + blob_coeffs(i,5)*sin(4.0*tmp2) + blob_coeffs(i,6)*sin(5.0*tmp2)
              else              !! blob is ellipse
                 r_mag = blob_coeffs(i,1)*blob_coeffs(i,2) &
                 /sqrt(blob_coeffs(i,2)*blob_coeffs(i,2)*cos(tmp2)**2. + &
                   blob_coeffs(i,1)*blob_coeffs(i,1)*sin(tmp2)**2.)                       
              end if                      
              if(temp.le.r_mag)then
                 keepgoing = .false.
              end if
           end do
           
           !! Are we too close to a boundary?
           if(x-0.5d0*dx.le.xb_min.or.x+0.5d0*dx.ge.xb_max.or.y-0.5d0*dx.le.yb_min.or.y+0.5d0*dx.ge.yb_max)then
              keepgoing = .false.
           end if            
                                  

        !! END CRITERIA
        !! --------------------------------------------------------------

        !! Place a particle here
        if(keepgoing) then
           ipart = ipart + 1
           call random_number(temp);temp = temp -0.5d0;xp(ipart) = pdp_x(j) + temp*dxmin*0.5d0
           call random_number(temp);temp = temp -0.5d0;yp(ipart) = pdp_y(j) + temp*dxmin*0.5d0
           dxp(ipart) = dx     
        end if           
        
        !! Deactive all pdps within dx of this pdp
        !! Search down
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i-1;
           if(i.eq.0) then
              idown = 0
              thdown = -0.49999999d0*pi
              keepgoing = .false.
           else if(pdp_dist2(i).ge.dx*dx) then             
              idown = i
              thdown = atan2((pdp_y(idown)-y),(pdp_x(idown)-x))
              keepgoing = .false.    
 
           end if       
        end do
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i+1
           if(i.eq.npdps+1) then
              iup = npdps+1
              thup = 0.4999999999d0*pi
              keepgoing = .false.              
           else if(pdp_dist2(i).ge.dx*dx) then 
              iup = i
              thup = atan2((pdp_y(iup)-y),(pdp_x(iup)-x))              
              keepgoing = .false.
           
           end if
        end do     
        
        !! Temporary store for the new pdps
        th_increment = (thup-thdown)/5.0d0
        inew = 0
        do i=1,5
           temp = y + dx*sin(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)    
           if(temp.ge.yb_min.and.temp.le.yb_max) then
              inew = inew + 1
              pdp_x_new(inew) = x + dx*cos(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)
              pdp_y_new(inew) = temp
           end if            
        end do
 

        block_new = idown + 1
        block_right = block_new + inew
        block_delete = iup - idown - 1
        npdps_new = npdps + inew - block_delete

        !! Shunt indices above pdp of interest
        if(block_delete.gt.inew) then !! Shunt to LEFT
           do i=block_right,npdps_new
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do
        end if
        if(block_delete.lt.inew) then !! Shunt to RIGHT
           do i=npdps_new,block_right,-1
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do

        end if
                  
        !! Insert any new pdps
        if(inew.ne.0) then          
           pdp_x(block_new:block_new+inew-1) = pdp_x_new(1:inew)
           pdp_y(block_new:block_new+inew-1) = pdp_y_new(1:inew)        
        end if       

        npdps = npdps_new
                  
                
        !! How left is the left-most PDP?                 
        minpdp = minval(pdp_x(1:npdps))
                
        write(6,*) ipart,(minval(pdp_x(1:npdps))-xb_min)/(xb_max-xb_min)
     end do  
                         
                                                  
                        
     npfb = ipart
     dx0 = maxval(dxp(1:npfb))
     

     write(*,*) 'nb,npfb= ', nb,npfb,nbio
!! ------------------------------------------------------------------------------------------------
  case(8)
  !! Loads of cylinders!!...

     D_cyl = 1.0d0
     S_cyl = D_cyl*1.2d0
     h0=D_cyl/2.0d0      !cylinder radius
     yl=4.0d0*S_cyl ! box height
     xl=3.0d0*sqrt(3.0d0)*S_cyl ! channel length
     dx0=D_cyl/100.0       !75
     xbcond=1;ybcond=1     
     
     nb_patches = 4
     allocate(b_node(nb_patches,2),b_edge(nb_patches,2))
     allocate(b_type(nb_patches))
     b_type(:) = (/ 3, 3, 3, 3/)  
     b_node(1,:) = (/-0.5d0*xl, -0.5d0*yl /)
     b_node(2,:) = (/0.5d0*xl, -0.5d0*yl /)
     b_node(3,:) = (/0.5d0*xl, 0.5d0*yl /)
     b_node(4,:) = (/-0.5d0*xl, 0.5d0*yl /)
     nb_blobs = 31
     allocate(blob_centre(nb_blobs,2),blob_coeffs(nb_blobs,6),blob_rotation(nb_blobs),blob_ellipse(nb_blobs))
     blob_centre(1,:) = (/0.0d0,0.0d0/)   !! Column 0
     blob_centre(2,:) = (/0.0d0,-S_cyl/)
     blob_centre(3,:) = (/0.0d0,S_cyl/)
     blob_centre(4,:) = (/0.0d0,-2.0d0*S_cyl/)
     blob_centre(5,:) = (/0.0d0,2.0d0*S_cyl/)
        
     blob_centre(6,:) = (/-S_cyl*sqrt(3.0d0)/2.0d0,-0.5d0*S_cyl/) !! Column -1/2
     blob_centre(7,:) = (/-S_cyl*sqrt(3.0d0)/2.0d0,0.5d0*S_cyl/) 
     blob_centre(8,:) = (/-S_cyl*sqrt(3.0d0)/2.0d0,-1.5d0*S_cyl/) 
     blob_centre(9,:) = (/-S_cyl*sqrt(3.0d0)/2.0d0,1.5d0*S_cyl/)           
     
     blob_centre(10,:) = (/S_cyl*sqrt(3.0d0)/2.0d0,-0.5d0*S_cyl/) !! Column 1/2
     blob_centre(11,:) = (/S_cyl*sqrt(3.0d0)/2.0d0,0.5d0*S_cyl/) 
     blob_centre(12,:) = (/S_cyl*sqrt(3.0d0)/2.0d0,-1.5d0*S_cyl/) 
     blob_centre(13,:) = (/S_cyl*sqrt(3.0d0)/2.0d0,1.5d0*S_cyl/) 

     blob_centre(14,:) = (/S_cyl*sqrt(3.0d0),0.0d0/)   !! Column 1
     blob_centre(15,:) = (/S_cyl*sqrt(3.0d0),-S_cyl/)
     blob_centre(16,:) = (/S_cyl*sqrt(3.0d0),S_cyl/)
     blob_centre(17,:) = (/S_cyl*sqrt(3.0d0),-2.0d0*S_cyl/)
     blob_centre(18,:) = (/S_cyl*sqrt(3.0d0),2.0d0*S_cyl/)

     blob_centre(19,:) = (/-S_cyl*sqrt(3.0d0),0.0d0/)   !! Column -1
     blob_centre(20,:) = (/-S_cyl*sqrt(3.0d0),-S_cyl/)
     blob_centre(21,:) = (/-S_cyl*sqrt(3.0d0),S_cyl/)
     blob_centre(22,:) = (/-S_cyl*sqrt(3.0d0),-2.0d0*S_cyl/)
     blob_centre(23,:) = (/-S_cyl*sqrt(3.0d0),2.0d0*S_cyl/)


     blob_centre(24,:) = (/-3.0d0*S_cyl*sqrt(3.0d0)/2.0d0,-0.5d0*S_cyl/) !! Column -3/2
     blob_centre(25,:) = (/-3.0d0*S_cyl*sqrt(3.0d0)/2.0d0,0.5d0*S_cyl/) 
     blob_centre(26,:) = (/-3.0d0*S_cyl*sqrt(3.0d0)/2.0d0,-1.5d0*S_cyl/) 
     blob_centre(27,:) = (/-3.0d0*S_cyl*sqrt(3.0d0)/2.0d0,1.5d0*S_cyl/)           
     
     blob_centre(28,:) = (/3.0d0*S_cyl*sqrt(3.0d0)/2.0d0,-0.5d0*S_cyl/) !! Column 3/2
     blob_centre(29,:) = (/3.0d0*S_cyl*sqrt(3.0d0)/2.0d0,0.5d0*S_cyl/) 
     blob_centre(30,:) = (/3.0d0*S_cyl*sqrt(3.0d0)/2.0d0,-1.5d0*S_cyl/) 
     blob_centre(31,:) = (/3.0d0*S_cyl*sqrt(3.0d0)/2.0d0,1.5d0*S_cyl/)      
                          
     
     do i=1,nb_blobs
        blob_coeffs(i,:)=h0*(/1.0d0,1.0d0,0.0d0,0.0d0,0.0d0,0.0d0/);blob_rotation(i)=-pi/9.0d0;blob_ellipse(i)=1
     end do

     call make_boundary_edge_vectors

     xb_min = minval(b_node(:,1));xb_max = maxval(b_node(:,1));yb_min = minval(b_node(:,2));yb_max = maxval(b_node(:,2))

     varresratio = 3.0d0  !! Ratio for scaling near the solid objects
     dxmax = dx0  
     dxmin = dx0/varresratio
     dxb=dx0/varresratio;dx_in=2.0d0*dxmax;dx_out=dx_in  !! dx for solids and in/outs...!! Ratio for scaling far field...
     call make_boundary_particles
     call make_boundary_blobs               
     ipart = nb   
     
     b0 = 2.0d0;b1 = 40.0d0;b2 = 30.0d0
         
     !! Initialise a line of potential dot points   
     nsearch = ceiling(yb_max-yb_min)/dxmin/2.0d0
     allocate(pdp_x(10*nsearch),pdp_y(10*nsearch))
     npdps = nsearch
     y=yb_min-0.5d0*dxmin
     i=0
     do while (y.lt.yb_max)
        y = y + dxmin/2.0d0;i=i+1
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_x(i) = xb_min + temp
        call random_number(temp);temp = (temp -0.5d0)*dxmin;
        pdp_y(i) = y 
     end do          
     npdps=i         
     allocate(pdp_dist2(10*nsearch));pdp_dist2=0.0d0   
        
     minpdp = minval(pdp_x(1:npdps))
     do while (minpdp.le.xb_max)  !! Keep going until all PDPs have passed out the end of the domain
        
        !! Pick the left-most pdp
        j=minloc(pdp_x(1:npdps),DIM=1)
        keepgoing = .true.
        x=pdp_x(j);y=pdp_y(j)
        pdp_dist2(1:npdps) = (x-pdp_x(1:npdps))**2.0d0 + (y-pdp_y(1:npdps))**2.0d0
        
        !! How far are we from the boundaries?
        dist2bound = xb_max-xb_min + 100.0d0
        i=1
        sumdG=0.0d0;sumG=0.0d0
        do while(i.le.nb)
           tmpN(1) = x - xp(i);tmpN(2) = y - yp(i);temp = sqrt(dot_product(tmpN,tmpN))
           if(i.gt.nbio.and.temp.le.dist2bound) dist2bound = temp
           i=i+1
           if(temp.le.4.25*dxp(i-1)) keepgoing = .false. !! Too close to bound, leave it.           
        end do        
        
        !! And what is the spacing, based on dist2bound?
        if(dist2bound.le.b0*dx0) then  !! Close - set to dxmin
           dx = dxmin                     
        else if(dist2bound.le.b1*dx0)then  !! A bit further out, smoothly vary from dxmin to dx0
           dx = 0.5d0*(dx0+dxmin) - 0.5d0*(dx0-dxmin)*cos((dist2bound-b0*dx0)*pi/((b1-b0)*dx0))  
        else if(dist2bound.le.b1*dx0+b2*dx_in)then  !! Further still: linearly vary from dx0 to dx_in
           dx = dx0 + (dx_in-dx0)*((dist2bound-b1*dx0)/(b2*dx_in))
        else     !! Far out: set to dx_in
           dx = dx_in
        end if       
      
               
        !! Check whether to place a particle here, based on some criteria
        !! --------------------------------------------------------------
           !! Are we within the object!!?!?!
           do i=1,nb_blobs
              temp = sqrt((x-blob_centre(i,1))**2. + (y-blob_centre(i,2))**2.)
              if(x-blob_centre(i,1).ge.0.0d0)then
                 tmp2 = asin((y-blob_centre(i,2))/temp)
              else
                 tmp2 = pi-asin((y-blob_centre(i,2))/temp)
              endif              
              tmp2 = tmp2 - blob_rotation(i)
              if(blob_ellipse(i).eq.0)then !! blob is blob
                 r_mag = blob_coeffs(i,1) + blob_coeffs(i,2)*sin(tmp2) + blob_coeffs(i,3)*sin(2.0*tmp2) + &
                      blob_coeffs(i,4)*sin(3.0*tmp2) + blob_coeffs(i,5)*sin(4.0*tmp2) + blob_coeffs(i,6)*sin(5.0*tmp2)
              else              !! blob is ellipse
                 r_mag = blob_coeffs(i,1)*blob_coeffs(i,2) &
                 /sqrt(blob_coeffs(i,2)*blob_coeffs(i,2)*cos(tmp2)**2. + &
                   blob_coeffs(i,1)*blob_coeffs(i,1)*sin(tmp2)**2.)                       
              end if                      
              if(temp.le.r_mag)then
                 keepgoing = .false.
              end if
           end do
           
           !! Are we too close to a boundary?
           if(x-0.5d0*dx.le.xb_min.or.x+0.5d0*dx.ge.xb_max.or.y-0.5d0*dx.le.yb_min.or.y+0.5d0*dx.ge.yb_max)then
              keepgoing = .false.
           end if            
                                  

        !! END CRITERIA
        !! --------------------------------------------------------------

        !! Place a particle here
        if(keepgoing) then
           ipart = ipart + 1
           call random_number(temp);temp = temp -0.5d0;xp(ipart) = pdp_x(j) + temp*dxmin*0.5d0
           call random_number(temp);temp = temp -0.5d0;yp(ipart) = pdp_y(j) + temp*dxmin*0.5d0
           dxp(ipart) = dx     
        end if           
        
        !! Deactive all pdps within dx of this pdp
        !! Search down
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i-1;
           if(i.eq.0) then
              idown = 0
              thdown = -0.49999999d0*pi
              keepgoing = .false.
           else if(pdp_dist2(i).ge.dx*dx) then             
              idown = i
              thdown = atan2((pdp_y(idown)-y),(pdp_x(idown)-x))
              keepgoing = .false.    
 
           end if       
        end do
        i=j;keepgoing = .true.
        do while(keepgoing)
           i=i+1
           if(i.eq.npdps+1) then
              iup = npdps+1
              thup = 0.4999999999d0*pi
              keepgoing = .false.              
           else if(pdp_dist2(i).ge.dx*dx) then 
              iup = i
              thup = atan2((pdp_y(iup)-y),(pdp_x(iup)-x))              
              keepgoing = .false.
           
           end if
        end do     
        
        !! Temporary store for the new pdps
        th_increment = (thup-thdown)/5.0d0
        inew = 0
        do i=1,5
           temp = y + dx*sin(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)    
           if(temp.ge.yb_min.and.temp.le.yb_max) then
              inew = inew + 1
              pdp_x_new(inew) = x + dx*cos(thdown + 0.1*(thup-thdown) + dble(i-1)*th_increment)
              pdp_y_new(inew) = temp
           end if            
        end do
 

        block_new = idown + 1
        block_right = block_new + inew
        block_delete = iup - idown - 1
        npdps_new = npdps + inew - block_delete

        !! Shunt indices above pdp of interest
        if(block_delete.gt.inew) then !! Shunt to LEFT
           do i=block_right,npdps_new
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do
        end if
        if(block_delete.lt.inew) then !! Shunt to RIGHT
           do i=npdps_new,block_right,-1
              ii = i + block_delete - inew
              pdp_x(i) = pdp_x(ii);pdp_y(i) = pdp_y(ii)
           end do

        end if
                  
        !! Insert any new pdps
        if(inew.ne.0) then          
           pdp_x(block_new:block_new+inew-1) = pdp_x_new(1:inew)
           pdp_y(block_new:block_new+inew-1) = pdp_y_new(1:inew)        
        end if       

        npdps = npdps_new
                  
                
        !! How left is the left-most PDP?                 
        minpdp = minval(pdp_x(1:npdps))
                
        write(6,*) ipart,(minval(pdp_x(1:npdps))-xb_min)/(xb_max-xb_min)
     end do  
                         
                                                  
                        
     npfb = ipart
     dx0 = maxval(dxp(1:npfb))
     
     write(*,*) 'nb,npfb= ', nb,npfb,nbio

  
  

  end select
!! ------------------------------------------------------------------------------------------------
!! Re-order nodes (from left to right)
  
   write(6,*) "About to quicksort"
   call quicksort(xp,1,npfb)
   write(6,*) "Quicksorted nodes ordered increasing x"

!! ------------------------------------------------------------------------------------------------
  !      ** write data out **
  ! 
  open(13,file='./IPART')
  write(13,*) nb,npfb,dx0
  write(13,*) xb_min,xb_max,yb_min,yb_max
  write(13,*) xbcond,ybcond
  do i=1,npfb
     if(node_type(i).ge.0.and.node_type(i).le.2) then
        write(13,*) xp(i), yp(i),node_type(i),xnorm(i),ynorm(i),dxp(i)
     else
        write(13,*) xp(i), yp(i),999,0.0d0,0.0d0,dxp(i)
     end if
  end do
  close(13)
  !! we can safely deallocate here
  deallocate(xp, yp,dxp,node_type,xnorm,ynorm)

  deallocate(b_node,b_edge)
  deallocate(b_type)
     ! end boundary

     write(*,*) 'END of DATCLASS'
200  stop
   end program datgen

! ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
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
!! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
  subroutine swap_nodes(i,j)
     use kind_parameters
     use common_parameter
     use global_variables
     implicit none  
     integer :: i,j
     double precision :: tmp
     integer :: itmp

     !! xp is already swapped by sub-routine quicksort     
     tmp = yp(j);yp(j)=yp(i);yp(i)=tmp
     tmp = xnorm(j);xnorm(j)=xnorm(i);xnorm(i)=tmp
     tmp = ynorm(j);ynorm(j)=ynorm(i);ynorm(i)=tmp
     tmp = dxp(j);dxp(j)=dxp(i);dxp(i)=tmp
     itmp = node_type(j);node_type(j)=node_type(i);node_type(i)=itmp                    
     

     return
  end subroutine swap_nodes
!! +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
   subroutine make_boundary_edge_vectors
     use kind_parameters
     use common_parameter
     use global_variables
     implicit none
     integer(ikind) ib,ibp1

     do ib = 1,nb_patches ! loop over all boundary patches
        ibp1 = mod(ib,nb_patches) + 1   
        b_edge(ib,:) = b_node(ibp1,:) - b_node(ib,:)  ! calculate b_edge
     end do

     return 
   end subroutine make_boundary_edge_vectors
!! ------------------------------------------------------------------------------------------------
   subroutine make_boundary_particles
     use kind_parameters
     use common_parameter
     use global_variables 

     implicit none
     integer(ikind) :: ipart,ib,ibm1
     real(rkind) :: x,y,m_be,tmp,tmp2
     integer(ikind) :: nround,iround
     real(rkind),dimension(2) :: nrm
   
     ipart=0

     !! we only allocate memory when we need
     allocate(xp(npar), yp(npar))
     allocate(xnorm(npar),ynorm(npar))
     allocate(dxp(npar))
     allocate(node_type(npar));node_type=999     

     !! Wall particles
     do ib=1,nb_patches  ! loop over all boundary patches
        ibm1 = mod(ib+nb_patches-2,nb_patches)+1
        if(abs(b_type(ib)).ne.3)then  ! if it is a wall, inflow or outflow patch
           m_be = dsqrt(dot_product(b_edge(ib,:),b_edge(ib,:)))
           nrm(1)=-b_edge(ib,2)/m_be;nrm(2)=b_edge(ib,1)/m_be
           if(b_type(ib).eq.1) dxio = dx_in
           if(b_type(ib).eq.2) dxio = dx_out
           if(b_type(ib).eq.0) dxio = dx_wall
           tmp = 0.5d0*dxio/m_be
           do while(tmp.lt.1.0-1.0d-10)   ! move along the patch in increments of dx
              ipart = ipart + 1
              x = b_node(ib,1) + tmp*b_edge(ib,1)
              y = b_node(ib,2) + tmp*b_edge(ib,2)
              xp(ipart) = x;yp(ipart) = y
              tmp = tmp + dxio/m_be  ! note, in future we should allow for dx.ne.dy here
              xnorm(ipart) = nrm(1);ynorm(ipart) = nrm(2)
              dxp(ipart)=dxio
              node_type(ipart) = b_type(ib) !! Set node-type so we know if it's a wall, inflow or outflow...
           end do
        end if
     end do
     nbio = ipart
          
     nb=ipart    
!!
     write(6,*) 'no. of solid boundary particles: ',nb
     return
   end subroutine make_boundary_particles
!! ------------------------------------------------------------------------------------------------
   subroutine make_boundary_blobs
     use kind_parameters
     use common_parameter
     use global_variables 

     implicit none
     integer(ikind) :: ipart,ib
     real(rkind) :: x,y,x0,y0,tmp2,r_mag,th,th_sh,dxlocal,th0
     real(rkind) :: a0,a1,a2,a3,a4,a5,y_stretch
   
     ipart=nb


     if(nb_blobs.ne.0)then
        do ib=1,nb_blobs
          
           !! Blobby blobs
           if(blob_ellipse(ib).eq.0)then
              a0 = blob_coeffs(ib,1);a1 = blob_coeffs(ib,2);a2 = blob_coeffs(ib,3)
              a3 = blob_coeffs(ib,4);a4 = blob_coeffs(ib,5);a5 = blob_coeffs(ib,6)       
              x0 = blob_centre(ib,1);y0=blob_centre(ib,2)
              th = 0.0d0 !! theta
              do while(th.le.2.0d0*pi-0.25*dxb/a0)
                 !! Position node
                 th_sh = th - blob_rotation(ib)
                 r_mag = (a0 +a1*sin(th_sh) + a2*sin(2.0*th_sh) + a3*sin(3.0*th_sh) + a4*sin(4.0*th_sh) + a5*sin(5.0*th_sh))
                 x = r_mag*cos(th);y = r_mag*sin(th)

                 if(x+x0.ge.xb_min.and.x+x0.le.xb_max.and.y+y0.ge.yb_min.and.y+y0.le.yb_max)then              
                    ipart = ipart + 1 
                    xp(ipart)=x0 + x;yp(ipart)=y0 + y;dxp(ipart)=dxb
         
                    !! Normals
                    tmp2=a1*cos(th_sh)+2.0*a2*cos(2.0*th_sh)+3.0*a3*cos(3.0*th_sh) &
                         +4.0*a4*cos(4.0*th_sh)+5.0*a5*cos(5.0*th_sh)
                    tmp2 = -atan(tmp2/r_mag)
                    xnorm(ipart) = (x*cos(tmp2) - y*sin(tmp2))/r_mag
                    ynorm(ipart) = (x*sin(tmp2) + y*cos(tmp2))/r_mag
                    node_type(ipart) = 0 !! Identify it as a wall                    
                 end if                 
      
                 !! Increment angle...
                 th  = th + cos(tmp2)*dxb/r_mag     
              end do   
                                                       
           !! Elliptical blobs
           else
              a0 = blob_coeffs(ib,1);a1 = blob_coeffs(ib,2)        
              x0 = blob_centre(ib,1);y0=blob_centre(ib,2)
           
           
              th = 0.0d0 !! theta
              do while(th.le.2.0d0*pi-0.25*dxb/a0)

                 !! Position node
                 th_sh = th - blob_rotation(ib)
                 r_mag = a0*a1/sqrt(a1*a1*cos(th_sh)**2. + a0*a0*sin(th_sh)**2.)      
                 
                 x = r_mag*cos(th);y = r_mag*sin(th)

                 if(x+x0.ge.xb_min.and.x+x0.lt.xb_max.and.y+y0.ge.yb_min.and.y+y0.lt.yb_max)then              
                    ipart = ipart + 1 
                    xp(ipart)=x0 + x;yp(ipart)=y0 + y
                    dxp(ipart)=dxb
         
                    !! Normals
                    tmp2 = (a0*a1**3. -a1*a0**3.)*cos(th_sh)*sin(th_sh) &
                       /(a1*a1*cos(th_sh)**2. + a0*a0*sin(th_sh)**2.)**(3./2.)
                    tmp2 = -atan(tmp2/r_mag)
                    xnorm(ipart) = (x*cos(tmp2) - y*sin(tmp2))/r_mag
                    ynorm(ipart) = (x*sin(tmp2) + y*cos(tmp2))/r_mag
                    node_type(ipart) = 0  !! Identify it as a wall
                 end if                 
      
                 !! Increment angle...
                 
                 th  = th + cos(tmp2)*dxb/r_mag     
              end do                                            
           end if                                                  
        end do              
     end if

     nb=ipart    
!!
     write(6,*) 'no. of solid boundary particles: ',nb
     return
   end subroutine make_boundary_blobs  
