module rhs
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
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
  use common_vars
  use derivatives
  use thermodynamics
  use characteristic_boundaries
  use chemistry
  implicit none
  
 
  private
  public calc_all_rhs,filter_variables

  !! Allocatable arrays for 1st and 2nd gradients
  real(rkind),dimension(:,:),allocatable :: gradu,gradv,gradw,gradlnro,gradp
  real(rkind),dimension(:,:),allocatable :: gradroE
  real(rkind),dimension(:,:),allocatable :: graddivvel
  real(rkind),dimension(:),allocatable :: lapu,lapv,lapw
  real(rkind),dimension(:,:),allocatable :: gradT,gradcp
  
  real(rkind),dimension(:),allocatable :: store_diff_E
  
  real(rkind),dimension(:),allocatable :: enth
  real(rkind),dimension(:,:),allocatable :: grad_enth
  
  real(rkind) :: dlnrodn,dlnrodt,dundn,dundt,dutdn,dutdt,dpdn,dpdt
  real(rkind) :: xn,yn,un,ut
  
  !! Characteristic boundary condition formulation
  real(rkind),dimension(:,:),allocatable :: L  !! The "L" in NSCBC formulation    
  

contains
!! ------------------------------------------------------------------------------------------------
  subroutine calc_all_rhs
     !! Control routine for calculating right hand sides. Does thermodynamic evaluations, finds
     !! gradients, and then calls property-specific RHS routines
     integer(ikind) :: k,i
     real(rkind) :: segment_tstart_rhs,segment_tend_rhs,segment_tstart_subtract,segment_tend_subtract
         
     !! Profiling
     segment_tstart_rhs = omp_get_wtime()     
     segment_tstart_subtract = segment_time_local(4) + segment_time_local(5) + segment_time_local(7)
  
     !! Some initial allocation of space for boundaries
#ifndef ms
     if(nb.ne.0) allocate(L(nb,5)) 
#else
     if(nb.ne.0) allocate(L(nb,5+nspec)) 
#endif  

     !! Determine secondary thermodynamic quantities and transport properties
     !! N.B. These are also calculated when setting time step, so not necessary for the first call
     !! to calc_all_rhs in the RK scheme.
     if(iRKstep.ne.1) then   
        call evaluate_temperature_and_pressure
        call evaluate_transport_properties      
     end if

     !! Initialise right hand sides to zero
     rhs_lnro=zero;rhs_u=zero;rhs_v=zero;rhs_w=zero;rhs_roE=zero;rhs_Yspec=zero
     
     !! Calculate derivatives of primary variables
     allocate(gradlnro(npfb,dims));gradlnro=zero  
     allocate(gradu(npfb,dims),gradv(npfb,dims),gradw(npfb,dims));gradw=zero
     call calc_gradient(lnro,gradlnro)     
     call calc_gradient(u,gradu)
     call calc_gradient(v,gradv)     
#ifdef dim3
     call calc_gradient(w,gradw)
#endif    
#ifndef isoT
     allocate(gradroE(npfb,dims))
     call calc_gradient(roE,gradroE)

     !! Temperature gradient... a lot depends on this
     allocate(gradT(npfb,dims))
     call calc_gradient(T,gradT)   
     
#endif          

     !! Call individual routines to build the RHSs
     !! N.B. second derivatives and derivatives of secondary variables are calculated within
     !! these subroutines
     call calc_rhs_lnro
     call calc_rhs_Yspec
     call calc_rhs_vel
     call calc_rhs_roE

     !! Calculate chemical production rates and add these to rhs of species equation
#ifdef react     
     call calculate_chemical_production_rate 
#endif     

     !! Evaluate RHS for boundaries
     if(nb.ne.0) then 
        call calc_rhs_nscbc   
#ifdef react
        deallocate(sumoverspecies_homega,reaction_rate_bound)
#endif                
     end if
     
     !! Clear space no longer required
     deallocate(gradlnro,gradu,gradv,gradw,gradp)
#ifndef isoT     
     deallocate(gradroE)
     deallocate(gradT)
#endif     
  
     !! Profiling - time spent doign RHS minus time spent doing gradients for RHS
     segment_tend_rhs = omp_get_wtime()
     segment_tend_subtract = segment_time_local(4) + segment_time_local(5) + segment_time_local(7)     
     segment_time_local(8) = segment_time_local(8) + segment_tend_rhs - segment_tstart_rhs
     segment_time_local(8) = segment_time_local(8) + segment_tstart_subtract - segment_tend_subtract

     return
  end subroutine calc_all_rhs
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_lnro
     !! Construct the RHS for lnro-equation
     integer(ikind) :: i,j
     real(rkind),dimension(dims) :: tmp_vec
     real(rkind) :: tmp_scal
     
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

     !! For any boundary nodes make transverse part of rhs
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
        end do
        !$omp end parallel do 
     end if       

     return
  end subroutine calc_rhs_lnro
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_rhs_Yspec
     !! Construct the RHS for species Yspec equation
     real(rkind),dimension(:),allocatable :: lapYspec
     real(rkind),dimension(:,:,:),allocatable :: gradYspec
     real(rkind),dimension(:,:),allocatable :: grad2Yspec
     integer(ikind) :: i,j,ispec
     real(rkind),dimension(dims) :: tmp_vec,gradroMdiff,grad_enthalpy
     real(rkind) :: tmp_scal,lapYspec_tmp,tmpY,divroDgradY,enthalpy,dcpdT,cpispec,tmpro
     real(rkind),dimension(:,:),allocatable :: speciessum_roDgradY,speciessum_hgradY
     real(rkind),dimension(:),allocatable :: speciessum_divroDgradY,speciessum_hY

     !! Allocate space for gradients and stores
     allocate(store_diff_E(npfb));store_diff_E = zero
#ifndef isoT     
     allocate(gradcp(npfb,dims));gradcp = zero
#endif     
#ifdef ms     

     allocate(gradYspec(npfb,dims,nspec))
     allocate(lapYspec(npfb))
     allocate(speciessum_divroDgradY(npfb));speciessum_divroDgradY = zero
     allocate(speciessum_roDgradY(npfb,dims));speciessum_roDgradY = zero
     allocate(speciessum_hY(npfb));speciessum_hY = zero
     allocate(speciessum_hgradY(npfb,dims));speciessum_hgradY = zero
     
     !! Loop over all species
     do ispec=1,nspec
    
        !! Calculate gradient and Laplacian for Yspec for this species     
        call calc_gradient(Yspec(:,ispec),gradYspec(:,:,ispec))
        call calc_laplacian(Yspec(:,ispec),lapYspec)
             
      
        !$omp parallel do private(i,tmp_vec,tmp_scal,tmpY,gradroMdiff,divroDgradY,enthalpy,grad_enthalpy, &
        !$omp dcpdT,cpispec,tmpro)
        do j=1,npfb-nb
           i=internal_list(j)
           tmpro = exp(-lnro(i)) !! tmpro contains 1/ro
           tmp_vec(1) = u(i);tmp_vec(2) = v(i);tmp_vec(3) = w(i)

           !! Advection term                    
           tmp_scal = -dot_product(tmp_vec,gradYspec(i,:,ispec))
         
           !! Molecular diffusion
#ifndef isoT           
           gradroMdiff(:) = r_temp_dependence*roMdiff(i,ispec)*gradT(i,:)/T(i)
#else
           gradroMdiff(:) = zero           
#endif           
           !!divroDgradY = ro*D*lapY + grad(ro*D)*gradY
           divroDgradY = roMdiff(i,ispec)*lapYspec(i) &
                      + dot_product(gradYspec(i,:,ispec),gradroMdiff(:))
           
           
           !! Sum roDgradY and div(roDgradY) over species
           speciessum_divroDgradY(i) = speciessum_divroDgradY(i) + divroDgradY
           speciessum_roDgradY(i,:) = speciessum_roDgradY(i,:) + roMdiff(i,ispec)*gradYspec(i,:,ispec)

#ifndef isoT                                 
           !! Evaluate enthalpy, cp and dcp/dT for species ispec
           call evaluate_enthalpy_at_node(T(i),ispec,enthalpy,cpispec,dcpdT)

           !! Add h*div.(ro*D*gradY) for species ispec to the energy diffusion store
           store_diff_E(i) = store_diff_E(i) + divroDgradY*enthalpy   
               
           !! Add gradh.ro*D*gradY for species ispec 
           grad_enthalpy = cpispec*gradT(i,:)
           store_diff_E(i) = store_diff_E(i) + roMdiff(i,ispec)* &
                                               dot_product(gradYspec(i,:,ispec),grad_enthalpy)    
                                               
           !! Add this species contribution to the mixture enthalpy
           speciessum_hY(i) = speciessum_hY(i) + enthalpy*Yspec(i,ispec) 
           
           !! Add this species contribution to hgradY
           speciessum_hgradY(i,:) = speciessum_hgradY(i,:) + enthalpy*gradYspec(i,:,ispec)
           
           !! Add this species contrib to gradcp(mix) = YdcpdT*gradT + cp(ispec)gradY
           gradcp(i,:) = gradcp(i,:) + Yspec(i,ispec)*dcpdT*gradT(i,:) &
                                     + cpispec*gradYspec(i,:,ispec)
#endif                    
                                          
           !! Add convective and diffusive terms to the RHS
           rhs_Yspec(i,ispec) = tmp_scal + divroDgradY*tmpro 
        end do
        !$omp end parallel do

        !! Make L5+ispec and boundary RHS
        if(nb.ne.0)then
           allocate(grad2Yspec(nb,dims))
           call calc_grad2bound(Yspec(:,ispec),grad2Yspec)
           !$omp parallel do private(i,xn,yn,un,ut,dutdt,lapYspec_tmp,tmpY &
           !$omp ,gradroMdiff,divroDgradY,enthalpy,grad_enthalpy,tmpro,cpispec,dcpdT)
           do j=1,nb
              i=boundary_list(j)
              tmpro = exp(-lnro(i))  !! tmpro contains 1/ro
              
              !! Evaluate gradient of molecular diffusivity              
#ifndef isoT                 
              gradroMdiff(:) = r_temp_dependence*roMdiff(i,ispec)*gradT(i,:)/T(i)
#else
              gradroMdiff(:) = zero
#endif                              
                           
              !! Slight difference here depending on type of boundary           
              if(node_type(i).eq.0)then  !! WALLS: no diffusive mass flux through walls

                 lapYspec_tmp = grad2Yspec(j,2) + grad2Yspec(j,3) !! Transverse terms only 
                 !!divroDgradY = ro*D*lapY + grad(ro*D)*gradY                 
                 divroDgradY = roMdiff(i,ispec)*lapYspec_tmp &
                            + dot_product(gradYspec(i,2:3,ispec),gradroMdiff(2:3))
           
              else            !! Inflow or outflow
                 lapYspec_tmp = grad2Yspec(j,1) + grad2Yspec(j,2) + grad2Yspec(j,3)
      
                 !!divroDgradY = ro*D*lapY + grad(ro*D)*gradY
                 divroDgradY = roMdiff(i,ispec)*lapYspec_tmp &                            
                            + dot_product(gradYspec(i,:,ispec),gradroMdiff(:))   
              end if              

              !! Sum roDgradY and div(roDgradY) over species
              speciessum_divroDgradY(i) = speciessum_divroDgradY(i) + divroDgradY
              speciessum_roDgradY(i,:) = speciessum_roDgradY(i,:) + roMdiff(i,ispec)*gradYspec(i,:,ispec)
              if(node_type(i).eq.0) speciessum_roDgradY(i,1) = zero !! zero normal contribution on walls

#ifndef isoT
              !! Evaluate enthalpy, cp and dcp/dT for species ispec
              call evaluate_enthalpy_at_node(T(i),ispec,enthalpy,cpispec,dcpdT)

              !! Add h*div.(ro*D*gradY) for species ispec to the energy diffusion store
              store_diff_E(i) = store_diff_E(i) + divroDgradY*enthalpy
          
              !! Add gradh(ispec).ro*D*gradY for species ispec 
              grad_enthalpy = cpispec*gradT(i,:)
              store_diff_E(i) = store_diff_E(i) + roMdiff(i,ispec)* &
                                                  dot_product(gradYspec(i,:,ispec),grad_enthalpy)

              !! Add this species contribution to the mixture enthalpy
              speciessum_hY(i) = speciessum_hY(i) + enthalpy*Yspec(i,ispec)

              !! Add this species contribution to hgradY
              speciessum_hgradY(i,:) = speciessum_hgradY(i,:) + enthalpy*gradYspec(i,:,ispec)
          
              !! Add this species contrib to gradcp(mix) = YdcpdT*gradT + cp(ispec)gradY
              gradcp(i,:) = gradcp(i,:) + Yspec(i,ispec)*dcpdT*gradT(i,:) &
                                           + cpispec*gradYspec(i,:,ispec)
#endif                   
                 
              !! Construct RHS (transverse convective and diffusive) (v(i)=w(i)=zero if WALL)
              rhs_Yspec(i,ispec) = -v(i)*gradYspec(i,2,ispec) - w(i)*gradYspec(i,3,ispec) &
                                    + divroDgradY*tmpro 

              !! Build the characteristic
              L(j,5+ispec) = u(i)*gradYspec(i,1,ispec)    !! (u(i)=zero if WALL)
              
           end do
           !$omp end parallel do 
           deallocate(grad2Yspec)
                       
        end if       

     end do

     
     !! Deallocate any stores no longer required    
     deallocate(lapYspec)
 
     !! Run through species again and finalise rhs and diffusion store for energy
     !$omp parallel do private(enthalpy,grad_enthalpy,tmpro,ispec)
     do i=1,npfb
        tmpro = exp(-lnro(i)) !! tmpro contains 1/ro
        do ispec=1,nspec                                

           !! Add the diffusion correction term to the rhs
           rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) - tmpro*Yspec(i,ispec)*speciessum_divroDgradY(i) &
                              - tmpro*dot_product(gradYspec(i,:,ispec),speciessum_roDgradY(i,:))          
           
        end do
        !! Additional terms for energy equation
#ifndef isoT                                     
        !! Add h*diffusion_correction_term to the energy diffusion store
        store_diff_E(i) = store_diff_E(i) - speciessum_divroDgradY(i)*speciessum_hY(i)
        
        !! Add sum_over_species(h*gradY).sum_over_species(roDgradY) 
        store_diff_E(i) = store_diff_E(i) - dot_product(speciessum_hgradY(i,:),speciessum_roDgradY(i,:))

        !! Add ro*cp(mix)*gradT*sum_over_species(DgradY) 
        grad_enthalpy = cp(i)*gradT(i,:)
        store_diff_E(i) = store_diff_E(i) - dot_product(grad_enthalpy,speciessum_roDgradY(i,:))
                                               
#endif                                               
     end do
     !$omp end parallel do          

     deallocate(speciessum_divroDgradY,speciessum_roDgradY,gradYspec)
     deallocate(speciessum_hY,speciessum_hgradY)
          
#endif          

     return
  end subroutine calc_rhs_Yspec    
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_vel
     !! Construct the RHS for u-equation
     integer(ikind) :: i,j
     real(rkind),dimension(dims) :: tmp_vec,gradvisc
     real(rkind) :: tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w
     real(rkind) :: tmpro,body_force_u,body_force_v,body_force_w,one_over_ro
     real(rkind) :: c
     real(rkind),dimension(:,:),allocatable :: grad2uvec,grad2ucross

     
     !! Gradient of velocity divergence
     allocate(graddivvel(npfb,dims));graddivvel=zero
     call calc_gradient(divvel,graddivvel)
          
     !! Allocate memory for spatial derivatives and stores
     allocate(lapu(npfb),lapv(npfb),lapw(npfb))
     lapu=zero;lapv=zero;lapw=zero
     allocate(gradp(npfb,dims));gradp=zero

     !! Calculate spatial derivatives
     call calc_laplacian(u,lapu)
     call calc_laplacian(v,lapv)
#ifdef dim3
     call calc_laplacian(w,lapw)          
#endif     

     !! Pressure gradient (method depends on whether isoT or not)
#ifndef isoT
     !! Evaluate pressure gradient directly (LABFM...)
     call calc_gradient(p-p_outflow,gradp) 
#else
     !! Evaluate the pressure gradient (from density gradient
     !$omp parallel do
     do i=1,npfb  
        gradp(i,:) = p(i)*gradlnro(i,:)  !! N.B. not precisely correct for isothermal multispec
     end do
     !$omp end parallel do
#endif                           
         
     !! Build RHS for internal nodes
     !$omp parallel do private(i,tmp_vec,tmp_scal_u,tmp_scal_v,tmp_scal_w,f_visc_u,f_visc_v,f_visc_w,tmpro &
     !$omp ,body_force_u,body_force_v,body_force_w,gradvisc,one_over_ro)
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
        
        !! Viscous forces due to non-uniform viscosity
#ifdef tdtp
        gradvisc(:) = r_temp_dependence*visc(i)*gradT(i,:)/T(i)
        f_visc_u = f_visc_u + gradvisc(1)*(fourthirds*gradu(i,1) - twothirds*(gradv(i,2)+gradw(i,3))) &
                            + gradvisc(2)*(gradu(i,2)+gradv(i,1)) &
                            + gradvisc(3)*(gradu(i,3)+gradw(i,1))
        f_visc_v = f_visc_v + gradvisc(1)*(gradu(i,2)+gradv(i,1)) &
                            + gradvisc(2)*(fourthirds*gradv(i,2) - twothirds*(gradu(i,1)+gradw(i,3))) &
                            + gradvisc(3)*(gradv(i,3)+gradw(i,2))
        f_visc_w = f_visc_w + gradvisc(1)*(gradu(i,3)+gradw(i,1)) &
                            + gradvisc(2)*(gradv(i,3)+gradw(i,2)) &
                            + gradvisc(3)*(fourthirds*gradw(i,3) - twothirds*(gradu(i,1)+gradv(i,2))) 
#endif        
      
        !! Local density 
        tmpro = exp(lnro(i));one_over_ro = one/tmpro

        !! Body force
        body_force_u = grav(1) + driving_force(1)*one_over_ro
        body_force_v = grav(2) + driving_force(2)*one_over_ro
        body_force_w = grav(3) + driving_force(3)*one_over_ro

        !! Store u.(F_visc + F_body/ro) for use in energy eqn later 
        store_diff_E(i) = store_diff_E(i) + u(i)*(f_visc_u + body_force_u*tmpro) &
                                          + v(i)*(f_visc_v + body_force_v*tmpro) &
                                          + w(i)*(f_visc_w + body_force_w*tmpro)
                        
                        
        !! RHS 
        rhs_u(i) = -tmp_scal_u - gradp(i,1)*one_over_ro + body_force_u + f_visc_u*one_over_ro 
        rhs_v(i) = -tmp_scal_v - gradp(i,2)*one_over_ro + body_force_v + f_visc_v*one_over_ro
#ifdef dim3
        rhs_w(i) = -tmp_scal_w - gradp(i,3)*one_over_ro + body_force_w + f_visc_w*one_over_ro
#else
        rhs_w(i) = zero
#endif                       
     end do
     !$omp end parallel do
         
     !! Make L1,L2,L3,L4,L5 and populate viscous + body force part of rhs' and save transverse
     !! parts of convective terms for later...
     if(nb.ne.0)then
     
        !! Evaluate d2u/dx2,d2v/dx2,d2w/dx2, and d2udxy, d2udxz
        allocate(grad2uvec(nb,3),grad2ucross(nb,2))
        call calc_grad2vecbound(u,v,w,grad2uvec)  
        call calc_grad2crossbound(gradu,grad2ucross)        
   
     
        !$omp parallel do private(i,tmpro,c,xn,yn,un,ut,f_visc_u,f_visc_v,body_force_u,body_force_v &
        !$omp ,dpdn,dundn,dutdn,gradvisc,one_over_ro,f_visc_w)
        do j=1,nb
           i=boundary_list(j)
           tmpro = exp(lnro(i));one_over_ro = one/tmpro
#ifndef isoT           
           c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
#else
           c=sqrt(csq)
#endif            
           dpdn = gradp(i,1)
           if(node_type(i).eq.0)then !! walls are in bound norm coords
              xn=rnorm(i,1);yn=rnorm(i,2)  !! Bound normals
              un = u(i)*xn + v(i)*yn; ut = -u(i)*yn + v(i)*xn  !! Normal and transverse components of velocity

              dundn = xn*gradu(i,1)+yn*gradv(i,1)
              dutdn = -yn*gradu(i,1)+xn*gradv(i,1)

              L(j,1) = half*(un-c)*(dpdn - tmpro*c*dundn) !! L1 
              L(j,2) = un*(tmpro*gradlnro(i,1) - gradp(i,1)/c/c)
              L(j,3) = un*dutdn !! L3 
              L(j,4) = un*(gradw(i,1)*xn + gradw(i,2)*yn) !! L4
              L(j,5) = half*(un+c)*(dpdn + tmpro*c*dundn) !! L5 
              
              !! Should evaluate viscous terms here, but a) rhs_u,v =0, and b) they aren't needed for
              !! subsequent energy equation 
              rhs_u(i) = zero
              rhs_v(i) = zero    
              rhs_w(i) = zero     
              
              !! Don't augment store_diff_E as velocity is zero
!              store_diff_E(i) = store_diff_E(i) + zero 
              
           else    !! In/out is in x-y coord system

              L(j,1) = half*(u(i)-c)*(dpdn - tmpro*c*gradu(i,1))
              L(j,2) = un*(tmpro*gradlnro(i,1) - gradp(i,1)/c/c)
              L(j,3) = u(i)*gradv(i,1)
              L(j,4) = u(i)*gradw(i,1)
              L(j,5) = half*(u(i)+c)*(dpdn + tmpro*c*gradu(i,1))

              !! Build initial stress divergence 
              f_visc_u = visc(i)*(lapu(i) + onethird*graddivvel(i,1)) 
              f_visc_v = visc(i)*(lapv(i) + onethird*graddivvel(i,2)) 
              f_visc_w = visc(i)*(lapw(i) + onethird*graddivvel(i,3))
#ifdef tdtp
               !! non-uniform viscosity terms. N.B. dtau_nn/dn=zero for both inflow and outflow, so applied here.
              gradvisc(:) = r_temp_dependence*visc(i)*gradT(i,:)/T(i)
              f_visc_u = f_visc_u + zero                           !! dtau_nn/dn=0
              f_visc_v = f_visc_v + gradvisc(1)*(gradu(i,2)+gradv(i,1)) &
                                  + gradvisc(2)*(fourthirds*gradv(i,2) - twothirds*(gradu(i,1)+gradw(i,3))) &
                                  + gradvisc(3)*(gradv(i,3)+gradw(i,2))
              f_visc_w = f_visc_w + gradvisc(1)*(gradu(i,3)+gradw(i,1)) &
                                  + gradvisc(2)*(gradv(i,3)+gradw(i,2)) &
                                  + gradvisc(3)*(fourthirds*gradw(i,3) - twothirds*(gradu(i,1)+gradv(i,2)))     
#endif 
              
              !! Zero various terms depending on inflow or outflow. This is done by subtracting terms as
              !! appropriate.
              if(node_type(i).eq.1) then !! inflow
                 !! dtau_nn/dn = 0
                 f_visc_u = f_visc_u + visc(i)*(twothirds*graddivvel(i,1) - two*grad2uvec(j,1))
              else        !! Outflow
                 !! dtau_nn/dn = 0
                 f_visc_u = f_visc_u + visc(i)*(twothirds*graddivvel(i,1) - two*grad2uvec(j,1))
                             
                 !! dtau_jn/dn = 0 for j=2,3 
                 f_visc_v = f_visc_v -grad2uvec(j,2) - grad2ucross(j,1)
                 f_visc_w = f_visc_w -grad2uvec(j,3) - grad2ucross(j,2)

                 !! zero non-uniform viscosity contribution to dtau_jn/dn for j=2,3 
#ifdef tdtp
                 f_visc_v = f_visc_v - gradvisc(1)*(gradu(i,2)+gradv(i,1)) 
                 f_visc_w = f_visc_w - gradvisc(1)*(gradu(i,3)+gradw(i,1)) 
#endif 
              
             
              end if
     
              !! Body force
              body_force_u = grav(1) + driving_force(1)*one_over_ro
              body_force_v = grav(2) + driving_force(2)*one_over_ro
              body_force_w = grav(3) + driving_force(3)*one_over_ro            
                
              !! Viscous dissipation term for energy equation   
              store_diff_E(i) = store_diff_E(i) + u(i)*(f_visc_u + body_force_u*tmpro) &
                                                + v(i)*(f_visc_v + body_force_v*tmpro) &
                                                + w(i)*(f_visc_w + body_force_w*tmpro)

              !! Transverse + visc + source terms only
              rhs_u(i) = -v(i)*gradu(i,2) - w(i)*gradu(i,3) + f_visc_u*one_over_ro + body_force_u  
              rhs_v(i) = -v(i)*gradv(i,2) - w(i)*gradv(i,3) - gradp(i,2)*one_over_ro + f_visc_v*one_over_ro + body_force_v
              rhs_w(i) = -v(i)*gradw(i,2) - w(i)*gradw(i,3) - gradp(i,3)*one_over_ro + f_visc_w*one_over_ro + body_force_w 
           end if
        end do
        !$omp end parallel do 
        deallocate(grad2uvec,grad2ucross)         
     end if
         
     !! Deallocate any stores no longer required
     deallocate(lapu,lapv,lapw)
     deallocate(graddivvel) 

     return
  end subroutine calc_rhs_vel
!! ------------------------------------------------------------------------------------------------  
  subroutine calc_rhs_roE
     !! Construct the RHS for energy equation. Do (almost) nothing if isothermal
#ifndef isoT     
     integer(ikind) :: i,j
     real(rkind) :: tmp_conv,tmp_visc,tmp_p,tmp_E
     real(rkind),dimension(dims) :: gradlambda
     real(rkind),dimension(:,:),allocatable :: grad2T
     real(rkind),dimension(:),allocatable :: lapT
     real(rkind) :: lapT_tmp
         
     !! Allocation and calculation of temperature laplacian     
     allocate(lapT(npfb))
     call calc_laplacian(T,lapT)
     
     !! Build RHS
     !$omp parallel do private(i,tmp_p,tmp_conv,tmp_visc,tmp_E,gradlambda)
     do j=1,npfb-nb
        i=internal_list(j)

        !! Convection term: -u.grad(roE)
        tmp_conv = -u(i)*gradroE(i,1) - v(i)*gradroE(i,2) - w(i)*gradroE(i,3)

        !! Pressure term: -div.(pu)
        tmp_p = -p(i)*divvel(i) - u(i)*gradp(i,1) - v(i)*gradp(i,2) - w(i)*gradp(i,3) 
        
        !! Energy dilation term: -roE*div.u
        tmp_E = -roE(i)*divvel(i)
     
        !! Viscous heating term: div.(tau u)
        tmp_visc = (fourthirds*gradu(i,1) - twothirds*gradv(i,2) - twothirds*gradw(i,3))*gradu(i,1) &
                 + (fourthirds*gradv(i,2) - twothirds*gradw(i,3) - twothirds*gradu(i,1))*gradv(i,2) &
                 + (fourthirds*gradw(i,3) - twothirds*gradu(i,1) - twothirds*gradv(i,2))*gradw(i,3) &
                 + (gradu(i,2)+gradv(i,1))**two + (gradu(i,3)+gradw(i,1))**two + (gradv(i,3)+gradw(i,2))**two        
        store_diff_E(i) = store_diff_E(i) + visc(i)*tmp_visc
        
        !! Thermal diffusion term: div.(lambda*gradT + heating due to molecular diffusion)
        store_diff_E(i) = store_diff_E(i) + lambda_th(i)*lapT(i)  
        gradlambda(:) = lambda_th(i)*r_temp_dependence*gradT(i,:)/T(i)  + lambda_th(i)*gradcp(i,:)/cp(i)
        store_diff_E(i) = store_diff_E(i) + dot_product(gradlambda(:),gradT(i,:))

                
        !! Build final RHS
        rhs_roE(i) = tmp_conv + tmp_p + tmp_E + store_diff_E(i)
     end do
     !$omp end parallel do
     
     
     !! Populate viscous and transverse parts of rhs
     if(nb.ne.0)then
        allocate(grad2T(nb,3))
        call calc_grad2bound(T,grad2T)          
        !$omp parallel do private(i,tmp_visc,gradlambda,lapT_tmp)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0)then !! Wall
              !! Viscous energy term: div.(tau u) (some bits omitted as u,v,w=0,0 and grad(u/v/w,2/3)=0 )
              tmp_visc = fourthirds*gradu(i,1)*gradu(i,1) &
                       + gradv(i,1)*gradv(i,1) &
                       + gradw(i,1)*gradw(i,1)
              store_diff_E(i) = store_diff_E(i) + visc(i)*tmp_visc

              lapT_tmp = grad2T(j,2) + grad2T(j,3)
#ifdef wall_isoT
              lapT_tmp = lapT_tmp + grad2T(j,1)   !! Isothermal walls can have heat flux in normal
#endif              

              !! Add thermal diffusion term: div.(lambda*gradT)
              store_diff_E(i) = store_diff_E(i) + lambda_th(i)*lapT_tmp
              gradlambda(:) = lambda_th(i)*r_temp_dependence*gradT(i,:)/T(i)  + lambda_th(i)*gradcp(i,:)/cp(i)
#ifdef wall_isoT
              store_diff_E(i) = store_diff_E(i) + dot_product(gradlambda(:),gradT(i,:))!! isoT-walls can have normal heat flux
#else
              store_diff_E(i) = store_diff_E(i) + dot_product(gradlambda(2:3),gradT(i,2:3))!! adiabatic
#endif              
                        
              !! RHS (visc + cond + transverse + source)
              rhs_roE(i) = -roE(i)*gradv(i,2) - roE(i)*gradw(i,3) &
                           -p(i)*gradv(i,2) - p(i)*gradw(i,3) &
                           + store_diff_E(i)             
           else !! Inflow/outflow  (is in the x-y coordinate system)
              !! Viscous energy term: div.(tau u)
              tmp_visc = (fourthirds*gradu(i,1) - twothirds*gradv(i,2) - twothirds*gradw(i,3))*gradu(i,1) &
                       + (fourthirds*gradv(i,2) - twothirds*gradw(i,3) - twothirds*gradu(i,1))*gradv(i,2) &
                       + (fourthirds*gradw(i,3) - twothirds*gradu(i,1) - twothirds*gradv(i,2))*gradw(i,3) &
                       + (gradu(i,2)+gradv(i,1))**two + (gradu(i,3)+gradw(i,1))**two &
                       + (gradv(i,3)+gradw(i,2))**two         
              store_diff_E(i) = store_diff_E(i) + visc(i)*tmp_visc
              
              !! Add thermal diffusion term: div.(lambda*gradT)
              store_diff_E(i) = store_diff_E(i) + lambda_th(i)*lapT(i)  
              gradlambda(:) = lambda_th(i)*r_temp_dependence*gradT(i,:)/T(i)  + lambda_th(i)*gradcp(i,:)/cp(i)
              store_diff_E(i) = store_diff_E(i) + dot_product(gradlambda(:),gradT(i,:))
                        
              !! RHS (visc + cond + transverse + source)
              rhs_roE(i) = - v(i)*gradroE(i,2) - roE(i)*gradv(i,2) - w(i)*gradroE(i,3) - roE(i)*gradw(i,3) &
                         - p(i)*gradv(i,2) - v(i)*gradp(i,2) - p(i)*gradw(i,3) - w(i)*gradp(i,3) &
                         + store_diff_E(i) 
           end if
        end do
        !$omp end parallel do
        deallocate(grad2T)
     end if
     deallocate(lapT)
     deallocate(gradcp)
               
#else
     rhs_roE(1:npfb) = zero
#endif

     !! Deallocation of spatial derivatives
     deallocate(store_diff_E)
           
     return
  end subroutine calc_rhs_roE
!! ------------------------------------------------------------------------------------------------
  subroutine calc_rhs_nscbc
    !! This routine asks boundaries module to prescribe L as required, then builds the final 
    !! rhs for each equation. It should only be called if nb.ne.0
    integer(ikind) :: i,j,ispec
    real(rkind) :: tmpro,c,tmp_scal,cv,gammagasm1
           
    !! Loop over boundary nodes and specify L as required
    !$omp parallel do private(i)
    do j=1,nb
       i=boundary_list(j)  
       
       !! WALL BOUNDARY 
       if(node_type(i).eq.0) then      
          call specify_characteristics_wall(j,L(j,:),gradv(i,:),gradw(i,:))

       !! INFLOW BOUNDARY
       else if(node_type(i).eq.1) then 
          call specify_characteristics_inflow(j,L(j,:),gradlnro(i,:),gradp(i,:),gradu(i,:),gradv(i,:),gradw(i,:))       
 
       !! OUTFLOW BOUNDARY 
       else if(node_type(i).eq.2) then   
          call specify_characteristics_outflow(j,L(j,:),gradp(i,:),gradu(i,:),gradv(i,:),gradw(i,:))       
       end if          
    end do
    !$omp end parallel do
    
    !! ==================================================================================
    !! Use L to update the rhs on boundary nodes
    !$omp parallel do private(i,tmpro,c,tmp_scal,cv,ispec,gammagasm1)
    do j=1,nb
       i=boundary_list(j)
       tmpro = exp(lnro(i))
#ifndef isoT       
       c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
       gammagasm1 = Rgas_mix(i)/(cp(i)-Rgas_mix(i))
#else
       c=sqrt(csq)
#endif       

       !! This quantity appears in rhs_lnro and rhs_roE, so save in tmp_scal
       tmp_scal = (c*c*L(j,2) + L(j,5) + L(j,1))/c/c
       
       !! RHS for density logarithm
       rhs_lnro(i) = rhs_lnro(i) - tmp_scal/tmpro

       !! Velocity components       
       rhs_u(i) = rhs_u(i) - (L(j,5)-L(j,1))/(tmpro*c)       
       rhs_v(i) = rhs_v(i) - L(j,3)
       rhs_w(i) = rhs_w(i) - L(j,4)

       !! Energy
#ifndef isoT
       cv = cp(i) - Rgas_mix(i)
       rhs_roE(i) = rhs_roE(i) - u(i)*(L(j,5)-L(j,1))/c - tmpro*v(i)*L(j,3) -tmpro*w(i)*L(j,4) &
                               - (roE(i)/tmpro - cv*T(i))*tmp_scal - (L(j,5)+L(j,1))/gammagasm1
#endif           
       
       !! Species mass fractions
#ifdef ms
       do ispec=1,nspec
          rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) - L(j,5+ispec) 
       end do
#endif          
    end do
    !$omp end parallel do
    !! De-allocate L
    deallocate(L)

    return  
  end subroutine calc_rhs_nscbc
!! ------------------------------------------------------------------------------------------------  
  subroutine filter_variables
     !! This routine calls the specific filtering routine (within derivatives module) for each
     !! variable - lnro,u,v,w,roE,Yspec - and forces certain values on boundaries as required.
     integer(ikind) :: ispec
      
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
   
     
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(3) = segment_time_local(3) + segment_tend - segment_tstart
     
     return
  end subroutine filter_variables  
!! ------------------------------------------------------------------------------------------------    
end module rhs
