module characteristic_boundaries
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2021 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines which take in characteristic waves and some
  !! gradient terms, and modify the characteristics to specify the desired boundary condition.

    
  !! Boundary framework follows a mixture of:: 
  !! Sutherland & Kennedy (2003)  <--------- probably most useful ref
  !! Yoo & Im (2007)              <--------- probably the formulation most closely followed...
  !! Coussement et al. (2012)   
  
  !! 
  use kind_parameters
  use common_parameter
  use common_vars
  use thermodynamics
  implicit none
contains
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_wall(j,Lchar,gradb_v,gradb_w)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     real(rkind),dimension(:),intent(in) :: gradb_v,gradb_w
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,psource,Ysource
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

     !! Store the sound speed
#ifndef isoT       
     c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
#else
     c=sqrt(csq)
#endif 

     !! ISOTHERMAL FLOWS
#ifdef isoT
     !Lchar(1) is outgoing, and so is unchanged
     Lchar(2) = zero !! as there is no entropy wave
     Lchar(3) = zero
     Lchar(4) = zero
     Lchar(5)= Lchar(1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro) 
#ifdef ms     
     Lchar(5+1:5+nspec) = zero
#endif     
#else            

#ifdef wall_isoT
     !! THERMAL FLOWS, ISOTHERMAL WALLS
     gammagasm1 = Rgas_mix(i)/(cp(i)-Rgas_mix(i))
     !Lchar(1) is outgoing, and so is unchanged    
     !Lchar(3) is zero, and unchanged
     !Lchar(4) is zero, and unchanged
     Lchar(5)= Lchar(1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro)
     
     !! Pressure source term:
     psource = - half*gammagasm1*sumoverspecies_homega(j)
     
     !! Loop over species
#ifdef react
     Ysource = zero
     do ispec = 1,nspec
        !Lchar(5+ispec) is zero (hopefully!) and unchanged
        Ysource = Ysource + one_over_molar_mass(ispec)* &
                          ( reaction_rate_bound(j,ispec) - &
                            Lchar(5+ispec))
     end do
     Ysource = Ysource*Rgas_universal/Rgas_mix(i)                                       
#else
     psource = zero;ysource=zero
#endif     
     Lchar(2) = gammagasm1*(Lchar(1)+Lchar(5))/c/c &
              + tmpro*psource/p(i) &            !! Pressure source terms
              + tmpro*Ysource         !! species source terms
              !+ additional dT/dt source terms TBC
                           
#else          
     !! THERMAL FLOWS, ADIABATIC WALLS (imposed heat flux)
     !Lchar(1) is outgoing, and so is unchanged    
     !Lchar(2) is zero, and unchanged
     !Lchar(3) is zero, and unchanged
     !Lchar(4) is zero, and unchanged
     Lchar(5)= Lchar(1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro)
#ifdef ms     
     !Lchar(5+1:5+nspec) is zero, and unchanged
#endif     
#endif                             
#endif          

  end subroutine specify_characteristics_wall
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_inflow(j,Lchar,gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     real(rkind),dimension(:),intent(in) :: gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

     !! Store the sound speed
#ifndef isoT       
     c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
#else
     c=sqrt(csq)
#endif 

#ifdef isoT
   
#ifndef hardinf
     !! ISOTHERMAL FLOWS, PARTIALLY NON-REFLECTING
     !Lchar(1) is outgoing, so doesn't require modification
     Lchar(2) = zero
     Lchar(3) = v(i)*0.278d0*c/L_domain_x &        !! track v=zero
              + rhs_rov(i)                           !! rhs_v contains transverse and visc terms needed
     Lchar(4) = w(i)*0.278d0*c/L_domain_x &        !! track w=zero
              + rhs_row(i)                           !! rhs_w contains transverse and visc terms needed
     Lchar(5) = (u(i)-u_inflow)*0.278d0*(one-u(i)/c)*c*c*one/L_domain_x &     !! Track u_inflow
              - half*(v(i)*gradb_p(2)+p(i)*gradb_v(2)+tmpro*c*v(i)*gradb_u(2)) &    !! transverse 1 conv. terms
              - half*(w(i)*gradb_p(3)+p(i)*gradb_w(3)+tmpro*c*w(i)*gradb_u(3))      !! transverse 2 conv. terms 
#ifdef ms
     Lchar(5+1:5+nspec) = zero
#endif     
#else
     !! ISOTHERMAL FLOWS, HARD INFLOW
     !Lchar(1) is outgoing, so doesn't require modification
     Lchar(2) = zero
     Lchar(3) = zero        
     Lchar(4) = zero
     Lchar(5) = Lchar(1)   !! Acoustically reflecting
#ifdef ms
     Lchar(5+1:5+nspec) = zero          
#endif     
#endif
         
#else


#ifndef hardinf
    !! THERMAL FLOWS, PARTIALLY NON-REFLECTING
    gammagas = cp(i)/(cp(i)-Rgas_mix(i))
    gammagasm1 = gammagas - one
    Lchar(2) = (T_bound(j)-T(i))*c*0.278d0/L_domain_x/gammagas &  !! Track T_bound(j)
             - (v(i)*gradb_ro(2) + tmpro*gradb_v(2) + &
                v(i)*gradb_p(2)/c/c + gammagas*p(i)*gradb_v(2)/c/c) &
             - (w(i)*gradb_ro(3) + tmpro*gradb_w(3) + &
                w(i)*gradb_p(3)/c/c + gammagas*p(i)*gradb_w(3)/c/c)
             !! + visc + source terms TBC
    Lchar(3) = v(i)*0.278d0*c/L_domain_x &        !! track v=zero
             + rhs_rov(i)                           !! rhs_v contains transverse and visc terms needed
    Lchar(4) = w(i)*0.278d0*c/L_domain_x &        !! track w=zero
             + rhs_row(i)                           !! rhs_w contains transverse and visc terms needed

    Lchar(5) = (u(i)-u_inflow)*0.278d0*(one-u(i)/c)*c*c*one/L_domain_x &      !! Track u_inflow
             - half*(v(i)*gradb_p(2)+gammagas*p(i)*gradb_v(2)+tmpro*c*v(i)*gradb_u(2))  & !! transverse 1 conv. terms
             - half*(w(i)*gradb_p(3)+gammagas*p(i)*gradb_w(3)+tmpro*c*w(i)*gradb_u(3))    !! transverse 2 conv. terms  

    !! Add source terms for chemical reactions
#ifdef react
    Lchar(2) = Lchar(2) + (gammagas-one)*sumoverspecies_homega(j)/c/c
    Lchar(5) = Lchar(5) - half*(gammagas-one)*sumoverspecies_homega(j)
#endif
             
#ifdef ms
    Lchar(5+1:5+nspec) = zero                 
#endif    
#else
    !! THERMAL FLOWS, HARD INFLOW 
    !Lchar(1) is outgoing, so doesn't require modification
    Lchar(5) = Lchar(1)
    Lchar(3) = zero        
    Lchar(4) = zero
    
    !! Fixed temperature option       
!    Lchar(2) = gammagasm1*(Lchar(1)+Lchar(5))/c/c &
!             - gammagasm1*tmpro*gradb_v(2) &   !! trans 1 term
!             - gammagasm1*tmpro*gradb_w(3)     !! Trans 2 term 
             !+dT/dy term?
             
    !! Fixed density (and hence mass flux) option             
    Lchar(2) = -Lchar(1)/c/c &
               -v(i)*gradb_ro(2) - tmpro*gradb_v(2) &  !! trans 1 term
               -w(i)*gradb_ro(3) - tmpro*gradb_w(3)    !! trans 2 term

#ifdef ms
    Lchar(5+1:5+nspec) = zero                                  
#endif    
#endif
#endif       
         

  end subroutine specify_characteristics_inflow
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_outflow(j,Lchar,gradb_p,gradb_u,gradb_v,gradb_w)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     real(rkind),dimension(:),intent(in) :: gradb_p,gradb_u,gradb_v,gradb_w
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)
     
     !! Store the sound speed
#ifndef isoT       
     c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
#else
     c=sqrt(csq)
#endif 

#ifdef isoT
     !! ISOTHERMAL FLOWS, PARTIALLY NON-REFLECTING
     if(u(i).lt.c) then
        Lchar(1) = (p(i)-p_outflow)*0.278d0*c*(one)/two/L_domain_x &                      !! track p_outflow
                 - (one-u(i)/c)*half*(v(i)*gradb_p(2) + &
                                      p(i)*gradb_v(2)-tmpro*c*v(i)*gradb_u(2)) & !!transverse 1 conv. terms
                 - (one-u(i)/c)*half*(w(i)*gradb_p(3) + &
                                      p(i)*gradb_w(3)-tmpro*c*w(i)*gradb_u(3))   !! transverse 2 conv. terms
     end if
     Lchar(2) = zero   !! No entropy in isothermal flows
     !Lchar(3) is outgoing
     !Lchar(4) is outgoing
     !Lchar(5) is outgoing
     !Lchar(5+1:5+nspec) is outgoing          
        
#else
     !! THERMAL FLOWS, PARTIALLY NON-REFLECTING
     if(u(i).le.c) then !! Subsonic. If supersonic, just use L1 from definition...
     
        gammagas = cp(i)/(cp(i)-Rgas_mix(i))     
        Lchar(1) = (p(i)-p_outflow)*outflow_coeff*c*(one)/two/L_domain_x &                               !! track p_outflow
                 - (one-u(i)/c)*half*(v(i)*gradb_p(2)+gammagas*p(i)*gradb_v(2) - &
                                      tmpro*c*v(i)*gradb_u(2)) & !! trans1 conv.
                 - (one-u(i)/c)*half*(w(i)*gradb_p(3)+gammagas*p(i)*gradb_w(3) - &
                                      tmpro*c*w(i)*gradb_u(3))   !! trasn2 conv.
     
       !! Add source terms for reacting flows
#ifdef react
        Lchar(1) = Lchar(1) - half*(gammagas-one)*sumoverspecies_homega(j)
#endif
      
if(iRKstep.eq.1.and.j.eq.1) write(711,*) time,p(i)              
     end if
     !Lchar(2) is outgoing
     !Lchar(3) is outgoing
     !Lchar(4) is outgoing
     !Lchar(5) is outgoing
     !Lchar(5+1:5+nspec) is outgoing          
#endif            

     !! Sometimes we have a little inflow at an outflow boundary. In this case, set Lchar(3)=Lchar(4)=zero
     !! to suppress shear 
     if(u(i).le.zero) then
        Lchar(2)=zero !! no incoming entropy
        Lchar(3)=zero !! no incoming shear if outflow velocity is zero...
        Lchar(4)=zero
#ifdef ms
        Lchar(6:5+nspec)=zero !! No incoming composition waves??
#endif        
     end if   
         

  end subroutine specify_characteristics_outflow
!! ------------------------------------------------------------------------------------------------
  subroutine apply_time_dependent_bounds
     integer(ikind) :: i,j
  
     !! Loop over all boundary nodes
     !$omp parallel do private(i)
     do j=1,nb
        i=boundary_list(j)
        
        !! Wall boundaries
        if(node_type(i).eq.0) then
           !! In all cases, velocity on wall is zero
           u(i) = zero
           v(i) = zero
           w(i) = zero

           !! For isothermal walls, evaluate the energy given the prescribed temperature
#ifdef wall_isoT        
           call set_energy_on_bound(i,T_bound(j))        
#endif           
        
        !! Inflow boundaries
        else if(node_type(i).eq.1) then 
#ifdef hardinf
              !! Prescribed velocity for hardinflow
              u(i)=u_inflow
              v(i)=zero
              w(i)=zero        
              
              !! Optional fixed T or fixed ro.
#endif
        
        !! Outflow boundaries
        else if(node_type(i).eq.2) then
           !! Do nothing, generally
        end if
        
     end do
     !$omp end parallel do
  
  
     return
  end subroutine apply_time_dependent_bounds
!! ------------------------------------------------------------------------------------------------  
end module characteristic_boundaries
