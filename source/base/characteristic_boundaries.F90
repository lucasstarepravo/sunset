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

    
  !! Boundary framework follows (with newest taking precedence):: 
  !! Sutherland & Kennedy (2003)
  !! Yoo & Im (2007)
  !! Coussement et al. (2012)    
  use kind_parameters
  use common_parameter
  use common_vars
  use thermodynamics
  implicit none
contains
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_wall(j,Lchar,gradv,gradw)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     real(rkind),dimension(:),intent(in) :: gradv,gradw
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = exp(lnro(i))

     !! Store the sound speed
#ifndef isoT       
     c=calc_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
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
     Lchar(3) = zero
     Lchar(4) = zero
     Lchar(5)= Lchar(1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro) 
#ifdef ms
     Lchar(5+1:5+nspec) = zero     
#endif     
     Lchar(2) = gammagasm1*(Lchar(5)+Lchar(1))/c/c &
                 - gammagasm1*tmpro*gradv(2) - gammagasm1*tmpro*gradw(3) !! Compatible with fixing dT/dt=0?  
#else          
     !! THERMAL FLOWS, ADIABATIC WALLS (imposed heat flux)
     !Lchar(1) is outgoing, and so is unchanged    
     Lchar(2) = zero
     Lchar(3) = zero
     Lchar(4) = zero
     Lchar(5)= Lchar(1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro) 
#ifdef ms     
     Lchar(5+1:5+nspec) = zero     
#endif     
#endif                             
#endif          

  end subroutine specify_characteristics_wall
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_inflow(j,Lchar,gradlnro,gradp,gradu,gradv,gradw)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     real(rkind),dimension(:),intent(in) :: gradlnro,gradp,gradu,gradv,gradw
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = exp(lnro(i))

     !! Store the sound speed
#ifndef isoT       
     c=calc_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
#else
     c=sqrt(csq)
#endif 

#ifdef isoT
   
#ifndef hardinf
     !! ISOTHERMAL FLOWS, PARTIALLY NON-REFLECTING
     !Lchar(1) is outgoing, so doesn't require modification
     Lchar(2) = zero
     Lchar(3) = v(i)*0.278d0*c/L_domain_x &        !! track v=zero
              + rhs_v(i)                           !! rhs_v contains transverse and visc terms needed
     Lchar(4) = w(i)*0.278d0*c/L_domain_x &        !! track w=zero
              + rhs_w(i)                           !! rhs_w contains transverse and visc terms needed
     Lchar(5) = (u(i)-u_inflow)*0.278d0*(one-u(i)/c)*c*c*one/L_domain_x &     !! Track u_inflow
              - half*(v(i)*gradp(2)+p(i)*gradv(2)+tmpro*c*v(i)*gradu(2)) &    !! transverse 1 conv. terms
              - half*(w(i)*gradp(3)+p(i)*gradw(3)+tmpro*c*w(i)*gradu(3))      !! transverse 2 conv. terms 
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
             - (v(i)*gradlnro(2)*tmpro + tmpro*gradv(2) + v(i)*gradp(2)/c/c + gammagas*p(i)*gradv(2)/c/c) &
             - (w(i)*gradlnro(3)*tmpro + tmpro*gradw(3) + w(i)*gradp(3)/c/c + gammagas*p(i)*gradw(3)/c/c)
             !! + visc + source terms TBC
    Lchar(3) = v(i)*0.278d0*c/L_domain_x &        !! track v=zero
             + rhs_v(i)                           !! rhs_v contains transverse and visc terms needed
    Lchar(4) = w(i)*0.278d0*c/L_domain_x &        !! track w=zero
             + rhs_w(i)                           !! rhs_w contains transverse and visc terms needed

    Lchar(5) = (u(i)-u_inflow)*0.278d0*(one-u(i)/c)*c*c*one/L_domain_x &      !! Track u_inflow
             - half*(v(i)*gradp(2)+gammagas*p(i)*gradv(2)+tmpro*c*v(i)*gradu(2))  & !! transverse 1 conv. terms
             - half*(w(i)*gradp(3)+gammagas*p(i)*gradw(3)+tmpro*c*w(i)*gradu(3))    !! transverse 2 conv. terms  
#ifdef ms
    Lchar(5+1:5+nspec) = zero                 
#endif    
#else
    !! THERMAL FLOWS, HARD INFLOW 
    !Lchar(1) is outgoing, so doesn't require modification
    Lchar(5) = Lchar(1)
    Lchar(3) = zero        
    Lchar(4) = zero
    Lchar(2) = gammagasm1*(Lchar(1)+Lchar(5))/c/c &
             - gammagasm1*tmpro*gradv(2) &   !! trans 1 term
             - gammagasm1*tmpro*gradw(3)     !! Trans 2 term 
                 ! + dT/dy term?
#ifdef ms
    Lchar(5+1:5+nspec) = zero                                  
#endif    
#endif
#endif       
         

  end subroutine specify_characteristics_inflow
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_outflow(j,Lchar,gradp,gradu,gradv,gradw)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     real(rkind),dimension(:),intent(in) :: gradp,gradu,gradv,gradw
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = exp(lnro(i))

     !! Store the sound speed
#ifndef isoT       
     c=calc_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
#else
     c=sqrt(csq)
#endif 

#ifdef isoT
     !! ISOTHERMAL FLOWS, PARTIALLY NON-REFLECTING
     if(u(i).lt.c) then
        Lchar(1) = (p(i)-p_outflow)*0.278d0*c*(one)/two/L_domain_x &                      !! track p_outflow
                 - (one-u(i)/c)*half*(v(i)*gradp(2)+p(i)*gradv(2)-tmpro*c*v(i)*gradu(2)) & !!transverse 1 conv. terms
                 - (one-u(i)/c)*half*(w(i)*gradp(3)+p(i)*gradw(3)-tmpro*c*w(i)*gradu(3))   !! transverse 2 conv. terms
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
        Lchar(1) = (p(i)-p_outflow)*0.278d0*c*(one)/two/L_domain_x &                               !! track p_outflow
                 - (one-u(i)/c)*half*(v(i)*gradp(2)+gammagas*p(i)*gradv(2)-tmpro*c*v(i)*gradu(2)) & !! trans1 conv.
                 - (one-u(i)/c)*half*(w(i)*gradp(3)+gammagas*p(i)*gradw(3)-tmpro*c*w(i)*gradu(3))   !! trasn2 conv.
                 
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
        Lchar(3)=zero !! no incoming shear if outflow velocity is zero...
        Lchar(4)=zero
     end if   
         

  end subroutine specify_characteristics_outflow
!! ------------------------------------------------------------------------------------------------
end module characteristic_boundaries
