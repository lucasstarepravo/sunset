module characteristic_boundaries
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2021 onwards     |Main developer    
  !! JRCK               |Nov 2022         |Including combustion terms                 
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines which take in characteristic waves and some
  !! gradient terms, and modify the characteristics to specify the desired boundary condition.

    
  !! Boundary framework follows a mixture of:: 
  !! SK03: Sutherland & Kennedy (2003) 
  !! YI07: Yoo & Im (2007)             
  
  !! 
  use kind_parameters
  use common_parameter
  use common_vars
  use thermodynamics
  implicit none
contains
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_isothermal_wall(j,Lchar)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,Ysource
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

#ifdef isoT
     !! ISOTHERMAL FLOWS 
     c=sqrt(csq)

     !Lchar(1) is outgoing, and so is unchanged
     !Lchar(2) = zero !! as there is no entropy wave
     !Lchar(3) = zero and unchanged
     !Lchar(4) = zero and unchanged
     
     !! Acoustic reflection
     Lchar(5)= Lchar(1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro) 
#ifdef ms     
     !Lchar(5+1:5+nspec) = zero and unchanged
#endif     

#else            
     !! THERMAL FLOWS

     c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i))      
     gammagasm1 = Rgas_mix(i)/(cp(i)-Rgas_mix(i))

     !Lchar(1) is outgoing, and so is unchanged    
     !Lchar(3) is zero, and unchanged
     !Lchar(4) is zero, and unchanged
     Lchar(5)= Lchar(1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro)
     
     !! Y source terms if any
     Ysource = zero     
#ifdef react
     !! Loop over species and build source terms
     do ispec = 1,nspec
        !Lchar(5+ispec) is zero (hopefully!) and unchanged
        Ysource = Ysource + one_over_molar_mass(ispec)* &
                          ( reaction_rate_bound(j,ispec) - &
                            Lchar(5+ispec))
     end do
     Ysource = Ysource*Rgas_universal/Rgas_mix(i)                                       
#endif     
   
     Lchar(2) = gammagasm1*(Lchar(1)+Lchar(5))/c/c &
#ifdef react
              + tmpro*gammagasm1*sumoverspecies_homega(j)/p(i) &              
#endif     
              + tmpro*Ysource         !! species source terms
              ! + (ro/T)*dT/dt
                           
#endif          

  end subroutine specify_characteristics_isothermal_wall
!! ------------------------------------------------------------------------------------------------  
  subroutine specify_characteristics_adiabatic_wall(j,Lchar)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,psource,Ysource

#ifdef isoT
     !! If this subroutine is called when running isothermal flows, it throws up an error
     write(6,*) "Compiled for isothermal flows, but requested adiabatic"
     write(6,*) "wall boundaries in control.in. Please change to "
     write(6,*) "isothermal walls. Stopping."     
#endif          

     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

     !! Store the sound speed
     c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
 
      
     !! THERMAL FLOWS, ADIABATIC WALLS (imposed heat flux)
     !Lchar(1) is outgoing, and so is unchanged    
     !Lchar(2) is zero, and unchanged - no entropy input
     !Lchar(3) is zero, and unchanged 
     !Lchar(4) is zero, and unchanged
    
     !! Acoustic reflections    
     Lchar(5)= Lchar(1) + tmpro*c*dot_product(rnorm(i,:),grav+driving_force/tmpro)
#ifdef ms     
     !Lchar(5+1:5+nspec) is zero, and unchanged (no surface reactions)
#endif     

  end subroutine specify_characteristics_adiabatic_wall  
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_soft_inflow(j,Lchar,gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w)
     !! Specify characteristic waves for a (partially) non-reflecting inflow boundary.
     !! Follows something between YI07 and SK03 in that we start from SK03, and add the transverse
     !! terms as in YI07, but not the viscous terms.
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     real(rkind),dimension(:),intent(in) :: gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

#ifdef isoT   
     !! ISOTHERMAL FLOWS
     c=sqrt(csq)     
     
     !Lchar(1) is outgoing, so doesn't require modification
     Lchar(2) = zero
     Lchar(3) = v(i)*nscbc_coeff*c/L_domain_x &        !! track v=zero
              - v(i)*gradb_v(2) - w(i)*gradb_v(3) - gradb_p(2)/ro(i) !! transverse terms
     Lchar(4) = w(i)*nscbc_coeff*c/L_domain_x &        !! track w=zero
              + rhs_row(i)                           !! rhs_w contains transverse and visc terms needed
     Lchar(5) = (u(i)-u_inflow_local(j))*nscbc_coeff*(one-u(i)/c)*c*c*one/L_domain_x &     !! Track u_inflow
              - half*(v(i)*gradb_p(2)+p(i)*gradb_v(2)+tmpro*c*v(i)*gradb_u(2)) &    !! transverse 1 conv. terms
              - half*(w(i)*gradb_p(3)+p(i)*gradb_w(3)+tmpro*c*w(i)*gradb_u(3))      !! transverse 2 conv. terms 
#ifdef ms
     Lchar(5+1:5+nspec) = zero
#endif     

         
#else
     !! THERMAL FLOWS    
     c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
     gammagas = cp(i)/(cp(i)-Rgas_mix(i))

     !Lchar(1) is outgoing, don't modify
     Lchar(2) = (T_bound(j)-T(i))*c*nscbc_coeff/L_domain_x/gammagas &  !! Track T_bound(j)
#ifdef react
              + (gammagas-one)*sumoverspecies_homega(j)/c/c &     !! Source terms
#endif
              - (v(i)*gradb_ro(2) + tmpro*gradb_v(2) + &                !! Transverse terms
                 v(i)*gradb_p(2)/c/c + gammagas*p(i)*gradb_v(2)/c/c) &
              - (w(i)*gradb_ro(3) + tmpro*gradb_w(3) + &
                 w(i)*gradb_p(3)/c/c + gammagas*p(i)*gradb_w(3)/c/c)
               
     Lchar(3) = v(i)*nscbc_coeff*c/L_domain_x &        !! track v=zero
              - v(i)*gradb_v(2) - w(i)*gradb_v(3) - gradb_p(2)/ro(i) !! transverse terms
     Lchar(4) = w(i)*nscbc_coeff*c/L_domain_x &        !! track w=zero
              - v(i)*gradb_w(2) - w(i)*gradb_w(3) - gradb_p(3)/ro(i) !! transverse terms

     Lchar(5) = (u(i)-u_inflow_local(j))*nscbc_coeff*(one-u(i)/c)*c*c*one/L_domain_x &      !! Track u_inflow
#ifdef react
              - half*(gammagas-one)*sumoverspecies_homega(j) &                 !! Source terms
#endif
              - half*(v(i)*gradb_p(2)+gammagas*p(i)*gradb_v(2)+tmpro*c*v(i)*gradb_u(2))  & !! transverse 1 conv. terms
              - half*(w(i)*gradb_p(3)+gammagas*p(i)*gradb_w(3)+tmpro*c*w(i)*gradb_u(3))    !! transverse 2 conv. terms  
             
#ifdef ms
     !! TBC, include tracking, transverse and source terms for flame-inflow interactions
     Lchar(5+1:5+nspec) = zero !! At present just assume inflow is far enough from flame...
#endif    

#endif                   
         

  end subroutine specify_characteristics_soft_inflow
!! ------------------------------------------------------------------------------------------------  
  subroutine specify_characteristics_hard_inflow(j,Lchar)
     !! Specifies characteristics for hard-inflows with prescribed temperature.
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagasm1,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

#ifdef isoT
     !! ISOTHERMAL FLOWS
     c=sqrt(csq)
     !Lchar(1) is outgoing, so doesn't require modification
     Lchar(2) = zero  !! No entropy for isothermal flows
     Lchar(3) = zero  !! v,w are prescribed
     Lchar(4) = zero
     Lchar(5) = Lchar(1)  !-tmpro*c*du_inflowdt !! Acoustically reflecting
#ifdef ms
     Lchar(5+1:5+nspec) = zero ! presume no inflow composition variation         
#endif     
         
#else
     !! THERMAL FLOWS, prescribed inflow temperature
     c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
     gammagas = cp(i)/(cp(i)-Rgas_mix(i))     

     !Lchar(1) is outgoing, so doesn't require modification

     !! Acoustically reflecting
     Lchar(5) = Lchar(1)  !- tmpro*c*du_inflowdt

     !! v and w are prescribed
     Lchar(3) = zero
     Lchar(4) = zero
    
     !! Fixed temperature option       
     Lchar(2) = (gammagas-one)*(Lchar(1)+Lchar(5))/c/c &
#ifdef react
              - tmpro*(gammagas-one)*sumoverspecies_homega(j)/p(i) !&
#endif     
!              - gammagasm1*tmpro*gradb_v(2) &   !! trans 1 term
!              - gammagasm1*tmpro*gradb_w(3)     !! Trans 2 term 
              ! + (ro/T)*dT/dt
              ! + sumoveralpha(dY/dt) terms
             
     !! Fixed density (and hence mass flux) option             
!     Lchar(2) = - (Lchar(1)+Lchar(5))/c/c !&
!                - v(i)*gradb_ro(2) - tmpro*gradb_v(2) &  !! trans 1 term
!                - w(i)*gradb_ro(3) - tmpro*gradb_w(3) &   !! trans 2 term
!                + (gammagas-one)*sumoverspecies_homega(j)/c/c
                !-dro/dt
#ifdef ms
     do ispec=1,nspec
        Lchar(5+1:5+nspec) = zero! reaction_rate_bound(j,ispec)
                           ! - dY(ispec)/dt
     end do
#endif    


#endif                   
         

  end subroutine specify_characteristics_hard_inflow  
!! ------------------------------------------------------------------------------------------------  
  subroutine specify_characteristics_inflow_outflow(j,Lchar,gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w)
     !! Specify characteristic waves for a (partially) non-reflecting inflow boundary.
     !! Follows something between YI07 and SK03 in that we start from SK03, and add the transverse
     !! terms as in YI07, but not the viscous terms.
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(in) :: gradb_ro,gradb_p,gradb_u,gradb_v,gradb_w     
     real(rkind),dimension(:),intent(inout) :: Lchar
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,c,gammagas
     
     !! Index of this boundary node
     i = boundary_list(j)

     !! Store the density
     tmpro = ro(i)

#ifdef isoT   
     !! ISOTHERMAL FLOWS
     c=sqrt(csq)     

     if(u(i).gt.zero) then !! Inflow options
    
        !Lchar(1) is outgoing, so doesn't require modification
        Lchar(2) = zero
        Lchar(3) = v(i)*nscbc_coeff*c/L_domain_x &      !! track v=zero
        Lchar(4) = w(i)*nscbc_coeff*c/L_domain_x       !! track w=zero
        Lchar(5) = (u(i)-u_inflow_local(j))*nscbc_coeff*(one-u(i)/c)*c*c*one/L_domain_x     !! Track u_inflow
#ifdef ms
        Lchar(5+1:5+nspec) = zero
#endif     
     else               !! Outflow options
        !Lchar(1) is outgoing, unchanged
        Lchar(2) = zero   !! No entropy in isothermal flows
        !Lchar(3) is outgoing
        !Lchar(4) is outgoing
        Lchar(5) (p(i)-p_inflow)*nscbc_coeff*c*(one)/two/L_domain_x     !! track p_outflow
        !Lchar(5+1:5+nspec) is outgoing     
     end if        
#else
     !! THERMAL FLOWS    
     c=evaluate_sound_speed_at_node(cp(i),Rgas_mix(i),T(i)) 
     gammagas = cp(i)/(cp(i)-Rgas_mix(i))

     if(u(i).gt.zero) then !! inflow options
        !Lchar(1) is outgoing, don't modify
        Lchar(2) = (T_bound(j)-T(i))*c*nscbc_coeff/L_domain_x/gammagas &  !! Track T_bound(j)
#ifdef react
                 + (gammagas-one)*sumoverspecies_homega(j)/c/c &     !! Source terms
#endif
                 - (v(i)*gradb_ro(2) + tmpro*gradb_v(2) + &                !! Transverse terms
                    v(i)*gradb_p(2)/c/c + gammagas*p(i)*gradb_v(2)/c/c) &
                 - (w(i)*gradb_ro(3) + tmpro*gradb_w(3) + &
                    w(i)*gradb_p(3)/c/c + gammagas*p(i)*gradb_w(3)/c/c)
               
        Lchar(3) = v(i)*nscbc_coeff*c/L_domain_x &        !! track v=zero
                 - v(i)*gradb_v(2) - w(i)*gradb_v(3) - gradb_p(2)/ro(i) !! transverse terms
        Lchar(4) = w(i)*nscbc_coeff*c/L_domain_x &        !! track w=zero
                 - v(i)*gradb_w(2) - w(i)*gradb_w(3) - gradb_p(3)/ro(i) !! transverse terms

        Lchar(5) = (p(i)-P_inflow)*nscbc_coeff*c*(one)/two/L_domain_x &      !! Track p_inflow
#ifdef react
                 - half*(gammagas-one)*sumoverspecies_homega(j) !&                 !! Source terms
#endif
!                 - half*(v(i)*gradb_p(2)+gammagas*p(i)*gradb_v(2)+tmpro*c*v(i)*gradb_u(2))  & !! transverse 1 conv. terms
!                 - half*(w(i)*gradb_p(3)+gammagas*p(i)*gradb_w(3)+tmpro*c*w(i)*gradb_u(3))    !! transverse 2 conv. terms  
             
#ifdef ms
        !! TBC, include tracking, transverse and source terms for flame-inflow interactions
        Lchar(5+1:5+nspec) = zero !! At present just assume inflow is far enough from flame...
#endif    

     else        !! Outflow options
        !Lchar(1) is outgoing
        !Lchar(2) is outgoing
        !Lchar(3) is outgoing
        !Lchar(4) is outgoing
        !Lchar(5+1:5+nspec) is outgoing 
        Lchar(5) = (p(i)-p_inflow)*nscbc_coeff*c*(one)/two/L_domain_x &        !! track p_outflow
#ifdef react
                 - half*(gammagas-one)*sumoverspecies_homega(j) &
#endif
                 + zero !! Neglecting transverse terms
     end if
#endif      
             
         

  end subroutine specify_characteristics_inflow_outflow
!! ------------------------------------------------------------------------------------------------
  subroutine specify_characteristics_outflow(j,Lchar)
     integer(ikind),intent(in) :: j
     real(rkind),dimension(:),intent(inout) :: Lchar
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
        Lchar(1) = (p(i)-p_outflow)*nscbc_coeff*c*(one)/two/L_domain_x! &                      !! track p_outflow
                 !- (one-u(i)/c)*half*(v(i)*gradb_p(2) + &
                 !                     p(i)*gradb_v(2)-tmpro*c*v(i)*gradb_u(2)) & !!transverse 1 conv. terms
                 !- (one-u(i)/c)*half*(w(i)*gradb_p(3) + &
                 !                     p(i)*gradb_w(3)-tmpro*c*w(i)*gradb_u(3))   !! transverse 2 conv. terms
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
        Lchar(1) = (p(i)-p_outflow)*nscbc_coeff*c*(one)/two/L_domain_x &                               !! track p_outflow
#ifdef react
                 - half*(gammagas-one)*sumoverspecies_homega(j) &
#endif
                 !! N.B. It's more stable to just follow Sutherland 2003 and neglect transverse terms
                 + zero!- (one-u(i)/c)*half*(v(i)*gradb_p(2)+gammagas*p(i)*gradb_v(2) - &
                 !                     tmpro*c*v(i)*gradb_u(2)) & !! trans1 conv.
                 !- (one-u(i)/c)*half*(w(i)*gradb_p(3)+gammagas*p(i)*gradb_w(3) - &
                 !                     tmpro*c*w(i)*gradb_u(3))   !! trasn2 conv.
    
      

     !Lchar(2) is outgoing
     !Lchar(3) is outgoing
     !Lchar(4) is outgoing
     !Lchar(5) is outgoing
     !Lchar(5+1:5+nspec) is outgoing          
     end if
#endif            

     !! Sometimes we have a little inflow at an outflow boundary. In this case, set Lchar(3)=Lchar(4)=zero
     !! to suppress shear 
     if(u(i).le.zero) then
        Lchar(2)=zero !! no incoming entropy
        Lchar(3)=zero !! no incoming shear if outflow velocity is zero...
        Lchar(4)=zero
#ifdef ms
        Lchar(5+1:5+nspec)=zero !! No incoming composition waves??
#endif        
     end if   
         

  end subroutine specify_characteristics_outflow
!! ------------------------------------------------------------------------------------------------
  subroutine apply_time_dependent_bounds
     integer(ikind) :: i,j,ispec
     
     segment_tstart = omp_get_wtime()
     
     !! Update time-dependant inflow velocity
!     call update_u_inflow
                     
     !! Loop over all boundary nodes
     !$omp parallel do private(i)
     do j=1,nb
        i=boundary_list(j)
        
        !! Wall boundaries
        if(node_type(i).eq.0) then
           !! In all cases, velocity on wall is zero
           rou(i) = zero
           rov(i) = zero
           row(i) = zero

           !! For isothermal walls, evaluate the energy given the prescribed temperature
           if(wall_type.eq.1) then
              call set_energy_on_bound(i,T_bound(j))        
           end if
        
        !! Inflow boundaries
        else if(node_type(i).eq.1) then 

           !! Hard or non-reflecting
           if(inflow_type.eq.1.or.inflow_type.eq.0)then
              do ispec = 1,nspec
                 Yspec(i,ispec) = Yspec_inflow(ispec)*ro(i)
              end do
           end if
        
        
           if(inflow_type.eq.1) then !! Hard inflow
              !! Prescribed velocity               
              rou(i)=u_inflow_local(j)*ro(i)
              rov(i)=zero
              row(i)=zero        

              ! Prescribed temperature              
              call set_energy_on_bound(i,T_bound(j))                      
           endif
           
           !! Don't need to do anything for soft inflow or inout
        
        !! Outflow boundaries
        else if(node_type(i).eq.2) then
           !! Do nothing, generally
        end if
        
     end do
     !$omp end parallel do
     
     
     !! For inflow-outflow types, adjust the zero-normal-flux flags dynamically depending on velocity sign
     if(inflow_type.eq.2) then
        !$omp parallel do private(i)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then  !! It's an inflow-outflow type
              if(u(i).gt.zero) then     !! Incoming flow, set diffusive flux flags for inflow
                 znf_mdiff(j) = .false.        
                 znf_tdiff(j) = .false.
                 znf_vdiff(j) = .true.     !! No normal viscous diffusion through inflows                
                 znf_vtdiff(j) = .false.                           
              else                      !! Outgoing flow, set diffusive flux flags for outflow
                 znf_mdiff(j) = .true.      !! No mass diffusion through outflow (N.B. not imposed...)
                 znf_tdiff(j) = .true.      !! No thermal diffusion through outflow
                 znf_vdiff(j) = .false.      
                 znf_vtdiff(j) = .true.      !! No tangential viscous diffusion through outflow       
              end if
           end if
        end do
        !$omp end parallel do     
     endif           

     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(2) = segment_time_local(2) + segment_tend - segment_tstart  
  
     return
  end subroutine apply_time_dependent_bounds 
!! ------------------------------------------------------------------------------------------------
  subroutine update_u_inflow
     !! Hard-coded routine to apply time-dependent inflow velocity. Currently hard-coded to ramp 
     !! the inflow over a prescribed time
     integer(ikind) :: i,j
     real(rkind) :: y,u_inflow_start,u_inflow_end
     real(rkind) :: ramp_time,u_inflow_mean
     
     !! Start and end inflow speeds, and ramp time
     u_inflow_start = u_char
     u_inflow_end = u_char
     ramp_time = 1.0d-4   
     
     !! Set the desired mean inflow velocity
     if(time.le.ramp_time) then
        u_inflow_mean = u_inflow_start + (u_inflow_end-u_inflow_start)*time/ramp_time
     else
        u_inflow_mean = u_inflow_end
     end if
     
     !! Only update u_inflow_local if it has changed (i.e. time<ramp_time)
     if(time.le.ramp_time) then
        !! Loop over all boundary nodes
        !$omp parallel do private(i,y)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then !! Inflows only
              y = rp(i,2)
              u_inflow_local(j) = u_inflow_mean*six*(half-y)*(half+y)
           end if
        end do
        !$omp end parallel do
     end if             
  
     return
  end subroutine update_u_inflow
!! ------------------------------------------------------------------------------------------------    
end module characteristic_boundaries
