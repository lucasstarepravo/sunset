module rhs
  !! This module contains routines to construct the RHS of all evolution equations
  !! These RHS' are built using calls to routines from the derivatives module, and
  !! calls to specific thermodynamics routines.
  
  !! Although separate subroutines, they must be called in correct order, as they rely
  !! on each other (e.g. divvel calculated in calc_rhs_lnro, also used in calc_rhs_vel)
  
  !! We use the lists internal_list and boundary_list to loop through internal and boundary
  !! nodes respectively. The L arrays for boundaries run from 1 to nb, and so use index j
  !! when within a loop over boundary nodes.
  use kind_parameters
  use common_parameter
  use common_2d
  use derivatives
  use thermodynamics
  use boundaries
  implicit none
  
 
  private
  public calc_all_rhs,filter_variables

  !! Allocatable arrays for 1st and 2nd gradients
  real(rkind),dimension(:,:),allocatable :: gradu,gradv,gradw,gradlnro,gradp
  real(rkind),dimension(:,:),allocatable :: graddivvel
  real(rkind),dimension(:),allocatable :: lapu,lapv,lapw
  real(rkind),dimension(:),allocatable :: lapYspec
  real(rkind),dimension(:,:),allocatable :: gradYspec
  
  real(rkind),dimension(:),allocatable :: store_diff_E
  
  real(rkind) :: dlnrodn,dlnrodt,dundn,dundt,dutdn,dutdt,dpdn,dpdt
  real(rkind) :: xn,yn,un,ut
  

contains
!! ------------------------------------------------------------------------------------------------
  subroutine calc_all_rhs
  
     !! Some initial allocation of space
#ifndef ms
     if(nb.ne.0) allocate(L(nb,5)) 
#else
     if(nb.ne.0) allocate(L(nb,5+nspec)) 
#endif  

     !! Initialise to zero
     rhs_lnro=zero;rhs_u=zero;rhs_v=zero;rhs_w=zero;rhs_roE=zero;rhs_Yspec=zero
     
     !! Call individual rhs routines
     call calc_rhs_lnro
     call calc_rhs_vel
     call calc_rhs_roE
     call calc_rhs_Yspec
     if(nb.ne.0) call calc_rhs_nscbc     
     
  
     return
  end subroutine calc_all_rhs
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_lnro
     !! Construct the RHS for lnro-equation
     integer(ikind) :: i,j
     real(rkind),dimension(dims) :: tmp_vec
     real(rkind) :: tmp_scal


     segment_tstart = omp_get_wtime()

  

     !! Allocate some arrays for gradients
     allocate(gradlnro(npfb,dims));gradlnro=zero  
!     allocate(divvel(npfb));divvel=zero
     allocate(gradu(npfb,dims),gradv(npfb,dims),gradw(npfb,dims));gradw=zero
      
     !! Evaluate gradients
     call calc_gradient(lnro,gradlnro)     
     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)     
#ifdef dim3
     call calc_gradient(w,gradw)
#endif     
     
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,tmp_scal)
     do j=1,npfb-nb
        i=internal_list(j)
        tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3)= w(i)
        tmp_scal = dot_product(tmp_vec,gradlnro(i,:))
!        divvel(i) = gradu(i,1) + gradv(i,2) + gradw(i,3)

        rhs_lnro(i) = -divvel(i) - tmp_scal  
     end do
     !$omp end parallel do

     !! For any boundary nodes, begin population of L2 and make transverse part of rhs
     if(nb.ne.0)then
        !$omp parallel do private(i,tmp_scal,xn,yn,un,ut,dutdt)
        do j=1,nb
           i=boundary_list(j)
           tmp_scal = exp(lnro(i))
           if(node_type(i).eq.0)then  !! in bound norm coords for walls
              xn=rnorm(i,1);yn=rnorm(i,2)
              un = u(i)*xn + v(i)*yn; ut = -u(i)*yn + v(i)*xn  !! Normal and transverse components of velocity           
              dutdt = -yn*gradu(i,2)+xn*gradv(i,2) !! Transverse derivative of transverse velocity...
              

              rhs_lnro(i) = - dutdt - gradw(i,3)
           else !! In x-y coords for inflow, outflow
              rhs_lnro(i) = -v(i)*gradlnro(i,2) - gradv(i,2) - w(i)*gradlnro(i,3) - gradw(i,3)
           end if
#ifndef isoT           
           L(j,2) = tmp_scal*gradlnro(i,1)  !! First bit of L2 (where non-isothermal)          
#endif           
        end do
        !$omp end parallel do 
     end if       

     !! Deallocate any stores no longer required      
#ifndef ms     
     if(nb.eq.0) deallocate(gradlnro)
#endif

     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(7) = segment_time_local(7) + segment_tend - segment_tstart
     return
  end subroutine calc_rhs_lnro
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_vel
     !! Construct the RHS for u-equation
     integer(ikind) :: i,j
     real(rkind),dimension(dims) :: tmp_vec
     real(rkind) :: tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w
     real(rkind) :: tmpro,body_force_u,body_force_v,body_force_w
     real(rkind) :: c


     segment_tstart = omp_get_wtime()          
     
     !! Gradient of velocity divergence
     allocate(graddivvel(npfb,dims));graddivvel=zero
     call calc_gradient(divvel,graddivvel)
          
     !! Determine secondary thermodynamic quantities and transport properties
     call pressure_from_primary_vars     
     call temp_from_primary_vars
     call visc_from_temp

     !! Allocate memory for spatial derivatives and stores
     allocate(lapu(npfb),lapv(npfb),lapw(npfb))
     lapu=zero;lapv=zero;lapw=zero
     allocate(gradp(npfb,dims));gradp=zero
     allocate(store_diff_E(npfb))

     !! Calculate spatial derivatives
     call calc_laplacian(u,lapu)
     call calc_laplacian(v,lapv)
#ifdef dim3
     call calc_laplacian(w,lapw)          
#endif     

     call calc_gradient(p,gradp)
         
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w,tmpro &
     !$omp ,body_force_u,body_force_v,body_force_w)
     do j=1,npfb-nb
        i=internal_list(j)
        tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i) !! tmp_vec holds (u,v,w) for node i
        tmp_scal_u = dot_product(tmp_vec,gradu(i,:))      !! convective term for u
        tmp_scal_v = dot_product(tmp_vec,gradv(i,:))      !! Convective term for v   
        tmp_scal_w = dot_product(tmp_vec,gradw(i,:))      !! Convective term for w    
        
        !! Viscous term (Lap(U) + (1/3)grad(div.U) formulation - avoids explicit calculation of cross derivs)
        f_visc_u = visc(i)*(lapu(i) + onethird*graddivvel(i,1))
        f_visc_v = visc(i)*(lapv(i) + onethird*graddivvel(i,2))
        f_visc_w = visc(i)*(lapw(i) + onethird*graddivvel(i,3))
      
        !! Local density 
        tmpro = exp(lnro(i))  

        !! Body force
        body_force_u = grav(1) + driving_force(1)/tmpro
        body_force_v = grav(2) + driving_force(2)/tmpro
        body_force_w = grav(3) + driving_force(3)/tmpro 

        !! Store u.(F_visc + F_body/ro) for use in energy eqn later 
        store_diff_E(i) = u(i)*(f_visc_u + body_force_u*tmpro) &
                        + v(i)*(f_visc_v + body_force_v*tmpro) &
                        + w(i)*(f_visc_w + body_force_w*tmpro)

        !! RHS 
        rhs_u(i) = -tmp_scal_u - gradp(i,1)/tmpro + body_force_u + f_visc_u/tmpro 
        rhs_v(i) = -tmp_scal_v - gradp(i,2)/tmpro + body_force_v + f_visc_v/tmpro
#ifdef dim3
        rhs_w(i) = -tmp_scal_w - gradp(i,3)/tmpro + body_force_w + f_visc_w/tmpro
#else
        rhs_w(i) = zero
#endif                       
     end do
     !$omp end parallel do

        
     !! Finalise L2, make L1,L3,L4, and populate viscous + body force part of rhs' and save transverse
     !! parts of convective terms for later...
     if(nb.ne.0)then
        !$omp parallel do private(i,tmpro,c,xn,yn,un,ut,f_visc_u,f_visc_v,body_force_u,body_force_v &
        !$omp ,dpdn,dundn,dutdn)
        do j=1,nb
           i=boundary_list(j)
           tmpro = exp(lnro(i))
           c=calc_sound_speed(p(i),tmpro) 
           dpdn = gradp(i,1)
           if(node_type(i).eq.0)then !! walls are in bound norm coords
              xn=rnorm(i,1);yn=rnorm(i,2)  !! Bound normals
              un = u(i)*xn + v(i)*yn; ut = -u(i)*yn + v(i)*xn  !! Normal and transverse components of velocity

              dundn = xn*gradu(i,1)+yn*gradv(i,1)
              dutdn = -yn*gradu(i,1)+xn*gradv(i,1)

              L(j,2) = zero !! Because the boundary is stationary, and L2 has u in it.

              L(j,1) = half*(un-c)*(dpdn - tmpro*c*dundn) !! L1 
              L(j,5) = half*(un+c)*(dpdn + tmpro*c*dundn) !! L5 
              L(j,3) = un*dutdn !! L3 
              L(j,4) = un*(gradw(i,1)*xn + gradw(i,2)*yn) !! L4

              !! Should evaluate viscous terms here, but a) rhs_u,v =0, and b) they aren't needed for
              !! subsequent energy equation 
              rhs_u(i) = zero
              rhs_v(i) = zero    
              rhs_w(i) = zero       
              
           else    !! In/out is in x-y coord system
#ifdef isoT           
              L(j,2) = zero !! No entropy wave in isothermal flow
#else           
              L(j,2) = (c*c*L(j,2) - dpdn*u(i))/c/c !! L2 now completely built
#endif
              L(j,1) = half*(u(i)-c)*(dpdn - tmpro*c*gradu(i,1))
              L(j,5) = half*(u(i)+c)*(dpdn + tmpro*c*gradu(i,1))
              L(j,3) = u(i)*gradv(i,1)
              L(j,4) = u(i)*gradw(i,1)

              !! Viscous
              f_visc_u = visc(i)*(lapu(i) + onethird*graddivvel(i,1))
              f_visc_v = visc(i)*(lapv(i) + onethird*graddivvel(i,2))
              f_visc_w = visc(i)*(lapw(i) + onethird*graddivvel(i,3))
             
     
              !! Body force
              body_force_u = grav(1) + driving_force(1)/tmpro
              body_force_v = grav(2) + driving_force(2)/tmpro
              body_force_w = grav(3) + driving_force(3)/tmpro              
                
              !! Viscous dissipation term for energy equation   
              store_diff_E(i) = u(i)*(f_visc_u + body_force_u*tmpro) &
                              + v(i)*(f_visc_v + body_force_v*tmpro) &
                              + w(i)*(f_visc_w + body_force_w*tmpro)

              !! Transverse + visc + source terms only
              rhs_u(i) = -v(i)*gradu(i,2) - w(i)*gradu(i,3) + f_visc_u/tmpro + body_force_u  
              rhs_v(i) = -v(i)*gradv(i,2) - w(i)*gradv(i,3) - gradp(i,2)/tmpro + f_visc_v/tmpro + body_force_v
              rhs_w(i) = -v(i)*gradw(i,2) - w(i)*gradw(i,3) - gradp(i,3)/tmpro + f_visc_w/tmpro + body_force_w              
           end if
        end do
        !$omp end parallel do          
     end if
         
     !! Deallocate any stores no longer required
     deallocate(lapu,lapv,lapw)
     deallocate(graddivvel) 
      
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(6) = segment_time_local(6) + segment_tend - segment_tstart
     
     return
  end subroutine calc_rhs_vel
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_rhs_roE
     !! Construct the RHS for energy equation. Do (almost) nothing if isothermal
#ifndef isoT     
     integer(ikind) :: i,j
     real(rkind) :: tmp_conv,tmp_visc,tmpro,tmp_p,tmp_E
     real(rkind),dimension(:,:),allocatable :: gradE,grad2T
     real(rkind),dimension(:),allocatable :: lapT
     
     segment_tstart = omp_get_wtime()
     
     !! Allocation and calculation of gradients     
     allocate(gradE(npfb,dims))
     allocate(lapT(npfb))
     call calc_gradient(roE,gradE)
     call calc_laplacian(T,lapT)
     
     
     !! Build RHS
     !$omp parallel do private(i,tmpro,tmp_p,tmp_conv,tmp_visc,tmp_E)
     do j=1,npfb-nb
        i=internal_list(j)
        tmpro = exp(lnro(i))

        !! Convection term: -u.grad(roE)
        tmp_conv = -u(i)*gradE(i,1) - v(i)*gradE(i,2) - w(i)*gradE(i,3)

        !! Pressure term: div.(pu)
        tmp_p = -p(i)*divvel(i) + u(i)*gradp(i,1) + v(i)*gradp(i,2) + w(i)*gradp(i,3) 
        
        !! Energy dilation term: -E*div.u
        tmp_E = -roE(i)*divvel(i)
     
        !! Viscous energy term: div.(tau u)
        tmp_visc = (fourthirds*gradu(i,1) - twothirds*gradv(i,2) - twothirds*gradw(i,3))*gradu(i,1) &
                 + (fourthirds*gradv(i,2) - twothirds*gradw(i,3) - twothirds*gradu(i,1))*gradv(i,2) &
                 + (fourthirds*gradw(i,3) - twothirds*gradu(i,1) - twothirds*gradv(i,2))*gradw(i,3) &
                 + (gradu(i,2)+gradv(i,1))**two + (gradu(i,3)+gradw(i,1))**two + (gradv(i,3)+gradw(i,2))**two        
        store_diff_E(i) = store_diff_E(i) + visc(i)*tmp_visc
        
        !! Thermal diffusion term: div.(lambda grad T)
        store_diff_E(i) = store_diff_E(i) + lambda_cond*lapT(i)  ! assuming const lambda_cond for now...
        
        !! Heat sink term
!        store_diff_E(i) = store_diff_E(i) + tmpro*u(i)*heat_sink_mag 
        
        !! Build final RHS
        rhs_roE(i) = tmp_conv + tmp_p + tmp_E + store_diff_E(i)
     end do
     !$omp end parallel do
     
     !! Populate viscous and transverse parts of rhs
     if(nb.ne.0)then
        allocate(grad2T(nb,3))
        call calc_grad2bound(T,grad2T)          
        !$omp parallel do private(i,tmp_visc)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0)then !! Wall
              !! Viscous energy term: div.(tau u) (some bits omitted as u,v,w=0,0 and grad(u/v/w,2/3)=0 )
              tmp_visc = fourthirds*gradu(i,1)*gradu(i,1) + gradv(i,1)*gradv(i,1) 
                        
              !! RHS (visc + cond + transverse + source). N.B. normal element of lap(T)=0... (adiabatic!)
              rhs_roE(i) = -roE(i)*gradv(i,2) - roE(i)*gradw(i,3) &
                           -p(i)*gradv(i,2) - p(i)*gradw(i,3) &
                           + visc(i)*tmp_visc + lambda_cond*(grad2T(j,2)+grad2T(j,3))              
           else !! Inflow/outflow  (is in the x-y coordinate system)
              !! Viscous energy term: div.(tau u)
              tmp_visc = (fourthirds*gradu(i,1) - twothirds*gradv(i,2) - twothirds*gradw(i,3))*gradu(i,1) &
                       + (fourthirds*gradv(i,2) - twothirds*gradw(i,3) - twothirds*gradu(i,1))*gradv(i,2) &
                       + (fourthirds*gradw(i,3) - twothirds*gradu(i,1) - twothirds*gradv(i,2))*gradw(i,3) &
                       + (gradu(i,2)+gradv(i,1))**two + (gradu(i,3)+gradw(i,1))**two &
                       + (gradv(i,3)+gradw(i,2))**two         
              store_diff_E(i) = store_diff_E(i) + visc(i)*tmp_visc
                        
              !! RHS (visc + cond + transverse + source)
              rhs_roE(i) = - v(i)*gradE(i,2) - roE(i)*gradv(i,2) - w(i)*gradE(i,3) - roE(i)*gradw(i,3) &
                         - p(i)*gradv(i,2) - v(i)*gradp(i,2) - p(i)*gradw(i,3) - w(i)*gradp(i,3) &
                         + store_diff_E(i) + lambda_cond*lapT(i)
           end if
        end do
        !$omp end parallel do
        deallocate(grad2T)
     end if
     deallocate(gradE,lapT)
     
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(5) = segment_time_local(5) + segment_tend - segment_tstart             
#else
     rhs_roE(1:npfb) = zero
#endif

     !! Deallocation of spatial derivatives
!     deallocate(divvel)
     deallocate(visc)
     deallocate(store_diff_E)
     if(nb.eq.0) deallocate(T,p)
     if(nb.eq.0) deallocate(gradp,gradu,gradv,gradw)
           
     return
  end subroutine calc_rhs_roE
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_Yspec
     !! Construct the RHS for species Yspec equation
#ifdef ms
     real(rkind),dimension(:,:),allocatable :: grad2Yspec
     integer(ikind) :: i,j,ispec
     real(rkind),dimension(dims) :: tmp_vec
     real(rkind) :: tmp_scal,lapYspec_tmp,tmpY

     segment_tstart = omp_get_wtime()     
     
     allocate(gradYspec(npfb,dims))
     allocate(lapYspec(npfb))
     
     do ispec=1,nspec
     
        call calc_gradient(Yspec(:,ispec),gradYspec)
        call calc_laplacian(Yspec(:,ispec),lapYspec)
          
        !$omp parallel do private(i,tmp_vec,tmp_scal,tmpY)
        do j=1,npfb-nb
           i=internal_list(j)
           tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i)
           tmp_scal = dot_product(tmp_vec,gradYspec(i,:))
         

           rhs_Yspec(i,ispec) = -tmp_scal + MD*(lapYspec(i) + dot_product(gradYspec(i,:),gradlnro(i,:)))
        end do
        !$omp end parallel do

        !! Make L5+ispec and boundary RHS
        if(nb.ne.0)then
           allocate(grad2Yspec(nb,dims))
           call calc_grad2bound(Yspec(:,ispec),grad2Yspec)

           !$omp parallel do private(i,tmp_scal,xn,yn,un,ut,dutdt,lapYspec_tmp,tmpY)
           do j=1,nb
              i=boundary_list(j)
              tmp_scal = exp(lnro(i))
                        
              if(node_type(i).eq.0)then  !! walls in bound norm coords

                 lapYspec_tmp = grad2Yspec(j,2) + grad2Yspec(j,3) !! Transverse terms only (no diffusion of Yspec through walls!)

                 !! diffusive terms only (u=0 on wall)
                 rhs_Yspec(i,ispec) = MD*(lapYspec_tmp + dot_product(gradYspec(i,2:3),gradlnro(i,2:3)))

                 L(j,5+ispec) = zero !! No transport through walls 
              else !! inflow/outflow in x-y coords
                 lapYspec_tmp = grad2Yspec(j,1) + grad2Yspec(j,2) + grad2Yspec(j,3)

                 rhs_Yspec(i,ispec) = -v(i)*gradYspec(i,2) - w(i)*gradYspec(i,3) &
                           + MD*(lapYspec_tmp + dot_product(gradYspec(i,:),gradlnro(i,:))) !! Transverse and diffusive only


                 L(j,5+ispec) = u(i)*gradYspec(i,1)                

              end if       
           end do
           !$omp end parallel do 
           deallocate(grad2Yspec)
        end if       


     end do
  

     !! Deallocate any stores no longer required    
     deallocate(gradYspec,lapYspec)
     if(nb.eq.0) deallocate(gradlnro)    

     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(8) = segment_time_local(8) + segment_tend - segment_tstart              
#endif          
     return
  end subroutine calc_rhs_Yspec  
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_nscbc
    !! This routine asks boundaries module to prescribe L as required, then builds the final 
    !! rhs for each equation. It should only be called if nb.ne.0
    
    !! Boundary framework follows (with newest taking precedence):: 
    !! Sutherland & Kennedy (2003)
    !! Yoo & Im (2007)
    !! Coussement et al. (2012)    
    integer(ikind) :: i,j,ispec
    real(rkind) :: tmpro,c,tmp_scal,cv
    
    segment_tstart = omp_get_wtime()
       
    !! Loop over boundary nodes and specify L as required
    !$omp parallel do private(i,tmpro,c,ispec)
    do j=1,nb
       i=boundary_list(j)
       tmpro = exp(lnro(i))
       c=calc_sound_speed(p(i),tmpro)   
       
       !! WALL BOUNDARY +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(node_type(i).eq.0) then 
          L(j,5)=L(j,1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro) 
          L(j,3) = zero
          L(j,4) = zero
       
!! ISOTHERMAL FLOWS =================================================
#ifdef isoT
          L(j,2) = zero
          
!! THERMAL FLOWS ====================================================
#else          
#ifdef wall_isoT
   !! ISOTHERMAL WALL =====================================
          L(j,2) = gammagasm1*(L(j,5)+L(j,1))/c/c &
                 - gammagasm1*tmpro*gradv(i,2) - gammagasm1*tmpro*gradw(i,3) !! Compatible with fixing dT/dt=0?  
#else          
   !! IMPOSED HEAT FLUX WALL ==============================
          L(j,2) = zero
#endif                             
#endif          
#ifdef ms
          do ispec=1,nspec
             L(j,5+ispec) = zero
          end do
#endif          
!! ==================================================================
       !! INFLOW BOUNDARY +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       else if(node_type(i).eq.1) then 
!! ISOTHERMAL FLOWS =================================================
#ifdef isoT
   !! (Partially-) NON-REFLECTING INFLOW ==================
#ifndef hardinf
          L(j,5) = (u(i)-u_inflow)*0.278d0*(one-u(i)/c)*c*c*one/(xmax-xmin) &      !! Track u_inflow
                 - half*(v(i)*gradp(i,2)+p(i)*gradv(i,2)+tmpro*c*v(i)*gradu(i,2)) &    !! transverse 1 conv. terms
                 - half*(w(i)*gradp(i,3)+p(i)*gradw(i,3)+tmpro*c*w(i)*gradu(i,3))      !! transverse 2 conv. terms 
          L(j,3) = v(i)*0.278d0*c/(xmax-xmin) &             !! track v=zero
                   + rhs_v(i)                                !! rhs_v contains transverse and visc terms needed
          L(j,4) = w(i)*0.278d0*c/(xmax-xmin) &             !! track w=zero
                   + rhs_w(i)                                !! rhs_w contains transverse and visc terms needed
   !! HARD INFLOW =========================================          
#else
          L(j,5) = L(j,1)
          L(j,3) = zero        
          L(j,4) = zero
#endif
          L(j,2) = zero
          
!! THERMAL FLOWS ====================================================
#else

   !! (Partially-) NON-REFLECTING INFLOW ==================
#ifndef hardinf
          L(j,5) = (u(i)-u_inflow)*0.278d0*(one-u(i)/c)*c*c*one/(xmax-xmin) &      !! Track u_inflow
                 - half*(v(i)*gradp(i,2)+gammagas*p(i)*gradv(i,2)+tmpro*c*v(i)*gradu(i,2))  & !! transverse 1 conv. terms
                 - half*(w(i)*gradp(i,3)+gammagas*p(i)*gradw(i,3)+tmpro*c*w(i)*gradu(i,3))    !! transverse 2 conv. terms  
          L(j,3) = v(i)*0.278d0*c/(xmax-xmin) &             !! track v=zero
                   + rhs_v(i)                                !! rhs_v contains transverse and visc terms needed
          L(j,4) = w(i)*0.278d0*c/(xmax-xmin) &             !! track w=zero
                   + rhs_w(i)                                !! rhs_w contains transverse and visc terms needed

          L(j,2) = (T0-T(i))*c*0.278d0/(xmax-xmin)/gammagas &
                 - (v(i)*gradlnro(i,2)*tmpro + tmpro*gradv(i,2) + v(i)*gradp(i,2)/c/c + gammagas*p(i)*gradv(i,2)/c/c) &
                 - (w(i)*gradlnro(i,3)*tmpro + tmpro*gradw(i,3) + w(i)*gradp(i,3)/c/c + gammagas*p(i)*gradw(i,3)/c/c)
                 !! Need to add on visc+source terms
   !! HARD INFLOW =========================================          
#else
          L(j,5) = L(j,1)
          L(j,3) = zero        
          L(j,4) = zero
          L(j,2) = gammagasm1*(L(j,1)+L(j,5))/c/c &
                 - gammagasm1*tmpro*gradv(i,2) &   !! trans 1 term
                 - gammagasm1*tmpro*gradw(i,3)     !! Trans 2 term 
                 ! + dT/dy term?
#endif
#endif     
!! ==================================================================     
#ifdef ms
          do ispec=1,nspec
             L(j,5+ispec) = zero
          end do
#endif

       !! OUTFLOW BOUNDARY ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       else if(node_type(i).eq.2) then   !! L1 incoming, L2,L3,L4,L5,L5+ispec outgoing from definition
!! ISOTHERMAL FLOWS =================================================
#ifdef isoT
          L(j,2) = zero
          if(u(i).le.c) then !! Subsonic. If supersonic, just use L1 from definition...
             L(j,1) = (p(i)-p_infinity)*0.278d0*c*(one)/two/(xmax-xmin) &               !! track p_infinity
                    - (one-u(i)/c)*half*(v(i)*gradp(i,2)+p(i)*gradv(i,2)-tmpro*c*v(i)*gradu(i,2)) & !!transverse 1 conv. terms
                    - (one-u(i)/c)*half*(w(i)*gradp(i,3)+p(i)*gradw(i,3)-tmpro*c*w(i)*gradu(i,3)) !! transverse 2 conv. terms
          end if
          if(u(i).le.zero) then
             L(j,3)=zero !! no incoming shear if outflow velocity is zero...
             L(j,4)=zero
          end if
!! THERMAL FLOWS ====================================================
#else
          if(u(i).le.c) then !! Subsonic. If supersonic, just use L1 from definition...
             L(j,1) = (p(i)-p_infinity)*0.278d0*c*(one)/two/(xmax-xmin) &               !! track p_infinity
                    - (one-u(i)/c)*half*(v(i)*gradp(i,2)+gammagas*p(i)*gradv(i,2)-tmpro*c*v(i)*gradu(i,2)) &!!trans1 conv.
                    - (one-u(i)/c)*half*(w(i)*gradp(i,3)+gammagas*p(i)*gradw(i,3)-tmpro*c*w(i)*gradu(i,3)) !! trasn2 conv.
          end if
          if(u(i).le.zero) then
             L(j,3)=zero !! no incoming shear if outflow velocity is zero...
             L(j,4)=zero
          end if
#endif            
!! ==================================================================
          !! Want reflecting, set L1=-L5
!          L(j,1) = -L(j,5) !! Reflecting!!
       end if          
    end do
    !$omp end parallel do
    
    !! Deallocate remaining gradient arrays
    deallocate(gradp,gradu,gradv,gradw,gradlnro)    


    !! ==================================================================================
    !! Use L to update the rhs on boundary nodes
    !$omp parallel do private(i,tmpro,c,xn,yn,un,ut,tmp_scal,cv,ispec)
    do j=1,nb
       i=boundary_list(j)
       tmpro = exp(lnro(i))
       c=calc_sound_speed(p(i),tmpro) 
       xn=rnorm(i,1);yn=rnorm(i,2)  !! Bound normals
       un = u(i)*xn + v(i)*yn; ut = -u(i)*yn + v(i)*xn  !! Normal and transverse components of velocity

       !! This quantity appears in rhs_lnro and rhs_roE, so save in tmp_scal
       tmp_scal = (c*c*L(j,2) + L(j,5) + L(j,1))/c/c
       
       !! The RHS of the density equation. Independent of BC type
       rhs_lnro(i) = rhs_lnro(i) - tmp_scal/tmpro
       
       !! WALL BOUNDARY +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       if(node_type(i).eq.0)then      
          rhs_u(i) = zero
          rhs_v(i) = zero
          rhs_w(i) = zero
#ifdef isoT          
          rhs_roE(i) = zero
#else          
#ifdef wall_isoT
          rhs_roE(i) = rhs_lnro(i)*roE(i)  !! This if isothermal wall...
#else
          rhs_roE(i) = rhs_roE(i) - (L(j,5)+L(j,1))/gammagasm1 !This if adiabatic wall
#endif
#endif
#ifdef ms
          do ispec=1,nspec
             rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) - L(j,5+ispec)
          end do
#endif          
       !! INFLOW BOUNDARY +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       else if(node_type(i).eq.1)then 
          rhs_u(i) = rhs_u(i) - (L(j,5)-L(j,1))/(tmpro*c)       
          rhs_v(i) = rhs_v(i) - L(j,3)
          rhs_w(i) = rhs_w(i) - L(j,4)

          cv = Rs0/gammagasm1
          rhs_roE(i) = rhs_roE(i) - u(i)*(L(j,5)-L(j,1))/c - tmpro*v(i)*L(j,3) - tmpro*w(i)*L(j,4) &
                                - (roE(i)/tmpro - cv*T(i))*tmp_scal - (L(j,5)+L(j,1))/gammagasm1         


#ifdef hardinf
          rhs_u(i) = zero
          rhs_v(i) = zero
          rhs_w(i) = zero
#endif
#ifdef ms
          do ispec=1,nspec
             rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) - L(j,5+ispec) 
          end do
#endif          
       !! OUTFLOW BOUNDARY ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
       else if(node_type(i).eq.2)then    
          rhs_u(i) = rhs_u(i) - (L(j,5)-L(j,1))/(tmpro*c)       
          rhs_v(i) = rhs_v(i) - L(j,3)
          rhs_w(i) = rhs_w(i) - L(j,4)
#ifdef isoT
          rhs_roE(i) = zero
#else                    
          cv = Rs0/gammagasm1
          rhs_roE(i) = rhs_roE(i) - u(i)*(L(j,5)-L(j,1))/c - tmpro*v(i)*L(j,3) -tmpro*w(i)*L(j,4) &
                                - (roE(i)/tmpro - cv*T(i))*tmp_scal - (L(j,5)+L(j,1))/gammagasm1
#endif           
#ifdef ms
          do ispec=1,nspec
             rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) - L(j,5+ispec) 
          end do
#endif          
       end if
    end do
    !$omp end parallel do
  
    !! De-allocate Secondary properties and L
    deallocate(T,p)
    deallocate(L)
    
    !! Profiling
    segment_tend = omp_get_wtime()
    segment_time_local(4) = segment_time_local(4) + segment_tend - segment_tstart
    return  
  end subroutine calc_rhs_nscbc
!! ------------------------------------------------------------------------------------------------  
  subroutine filter_variables
     !! This routine calls the specific filtering routine (within derivatives module) for each
     !! variable - lnro,u,v,w,roE,Yspec - and forces certain values on boundaries as required.
     integer(ikind) :: i,j,ispec
     real(rkind) :: tmp
      
     segment_tstart = omp_get_wtime()
        
     !! Filter variables
     call calc_filtered_var(lnro)
     call calc_filtered_var(u)
     call calc_filtered_var(v)
#ifdef dim3
     call calc_filtered_var(w)
#endif     
#ifndef isoT     
     call calc_filtered_var(roE)       
#endif     
#ifdef ms
     do ispec=1,nspec
       call calc_filtered_var(Yspec(:,ispec))
     end do
#endif   

     !! Reset values (velocity and sometimes Energy) on boundaries as required...  
     if(nb.ne.0)then   
        !$omp parallel do private(i)
        do j=1,nb
           i=boundary_list(j)
           
           if(node_type(i).eq.0)then       !! Walls
              u(i)=zero;v(i)=zero;w(i)=zero
#ifdef wall_isoT
              roE(i)=exp(lnro(i))*T_bound(j)*Rs0/gammagasm1   !! Keep roE such that T=T_bound
#endif
           else if(node_type(i).eq.1)then  !! Inflow   
#ifdef hardinf
              u(i)=u_inflow;v(i)=zero;w(i)=zero
              roE(i)=exp(lnro(i))*T_bound(j)*Rs0/gammagasm1   !! Keep roE such that T=T_bound
#else
              !! DO NOTHING              
#endif                            
           else if(node_type(i).eq.2)then  !! Outflow
              !! DO NOTHING
           end if
        end do
        !$omp end parallel do
     end if
     
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(3) = segment_time_local(3) + segment_tend - segment_tstart
     
     return
  end subroutine filter_variables  
!! ------------------------------------------------------------------------------------------------    
end module rhs
