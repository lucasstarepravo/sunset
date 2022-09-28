module chemistry
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2022 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to evaluate chemical production rates for reacting flows.
  !! If the flow is inert, it does nothing. 
  !! Reaction rates are added to the rhs_Yspec arrays.
  use kind_parameters
  use common_parameter
  use common_vars
  use thermodynamics
  implicit none


contains
#ifdef react
!! ------------------------------------------------------------------------------------------------
  subroutine calculate_chemical_production_rates
     integer(ikind) :: i,ispec
     real(rkind) :: tmpro,production_rate,arrhenius_rate
     
     !$omp parallel do private(ispec,tmpro,production_rate,arrhenius_rate)
     do i=1,npfb
        tmpro = exp(lnro(i))
!        do ispec=1,nspec

           arrhenius_rate = arrhenius_coefs(1,1) + arrhenius_coefs(1,2)*log(T(i)) - arrhenius_coefs(1,3)/T(i)
           arrhenius_rate = exp(arrhenius_rate)

           production_rate = arrhenius_rate*(tmpro*max(Yspec(i,1),zero)/molar_mass(1))**one !! molar production rate
           production_rate = production_rate*molar_mass(1) !! Mass production rate
!write(6,*) T(i),tmpro,Yspec(i,1),arrhenius_rate,production_rate           
!if(arrhenius_rate.ne.zero)then
!write(6,*) iproc,i,itime,arrhenius_rate,production_rate        
!end if
!if(itime.eq.1.and.production_rate.gt.1d-5) write(6,*) T(i),production_rate

           rhs_Yspec(i,1) = rhs_Yspec(i,1) - production_rate/tmpro
           rhs_Yspec(i,2) = rhs_Yspec(i,2) + production_rate/tmpro
!        end do
alpha_out(i) = production_rate
     end do
     !$omp end parallel do
     
  
     return
  end subroutine calculate_chemical_production_rates
!! ------------------------------------------------------------------------------------------------  
  subroutine calculate_chemical_production_rate
     integer(ikind) :: i,ispec,istep,jspec,jstep,ithirdbody
     real(rkind) :: tmpro,production_rate,arrhenius_rate,tmpT,logT
     real(rkind) :: arrhenius_rate_back,gibbs,logP0ovR0
     real(rkind) :: molar_production_rate,forward_rate,arrhenius_rate0
     real(rkind) :: backward_rate,third_body_conc,p_reduced,lnkf
     real(rkind),dimension(:),allocatable :: logYspec_local,Yspec_local
     real(rkind),dimension(:),allocatable :: rateYspec
     
     !! Part of gibbs term (saves time doing outside loop)
     logP0ovR0 = log(P_ref/Rgas_universal)
     
     !! space of Y and ln(Y) local
     allocate(logYspec_local(nspec));logYspec_local = zero
     allocate(Yspec_local(nspec));Yspec_local=one
     allocate(rateYspec(nspec));rateYspec = zero
     
     !! Loop over all nodes
     !$omp parallel do private(istep,jstep,ispec,jspec,tmpro,tmpT,logT,production_rate,arrhenius_rate, &
     !$omp molar_production_rate,arrhenius_rate_back,gibbs,forward_rate,ithirdbody,arrhenius_rate0, &
     !$omp backward_rate,logYspec_local,Yspec_local,third_body_conc,rateYspec,p_reduced,lnkf)
     do i=1,npfb
     
        !! Density, temperature
        tmpro = exp(lnro(i))        
        tmpT = T(i)
        logT = log(tmpT)
        
        !! Store logarithm of species mass fraction for each species
        do ispec = 1,nspec
           Yspec_local(ispec) = max(Yspec(i,ispec),verysmall)
           logYspec_local(ispec) = log(Yspec_local(ispec))
        end do
        
        !! Initialise rate to zero
        rateYspec(:) = zero
        
        !! Loop over all steps ================================================
        do istep = 1,nsteps

           !! Forward rate ====================================================

           !! Evaluate Arrhenius rate for this step
           arrhenius_rate = arrhenius_coefs(istep,1) &
                          + arrhenius_coefs(istep,2)*logT &
                          - arrhenius_coefs(istep,3)/tmpT

           !! Loop over all reactants and build log of forward rate
           forward_rate = arrhenius_rate
           !! forward_rate contains ln(k_{f,m})
           
           
           !! Third bodies ====================================================
           if(third_body_flag(istep).ne.0) then
              ithirdbody = third_body_flag(istep) !! This is the type of third body for step istep
           
              !! Build third-body concentration
              third_body_conc = zero
              do ispec = 1,nspec
                 third_body_conc = third_body_conc + third_body_efficiencies(ithirdbody,ispec)* &
                                                     Yspec_local(ispec)/ &
                                                     molar_mass(ispec)
              end do
              third_body_conc = third_body_conc*tmpro !! Multiply out by ro
           else
              third_body_conc = one
           end if     
           
           !! Lindemann steps =================================================
           if(lindemann_form_flag(istep).ne.0) then
              jstep = lindemann_form_flag(istep)  !! This is the jstep-th Lindemann step
              
              !! Evaluate k0
              arrhenius_rate0 = lindemann_coefs(jstep,1) + &
                                lindemann_coefs(jstep,2)*logT - &
                                lindemann_coefs(jstep,3)/tmpT

              !! Reduced pressure
              p_reduced = exp(arrhenius_rate0 - arrhenius_rate)*third_body_conc
              
              !! rate constant for this step
              forward_rate = arrhenius_rate + &
                             log(p_reduced/(one+p_reduced)) + &
                             lindemann_coefs(jstep,4)              
           else
              forward_rate = arrhenius_rate
           end if   
           
           !! Store log(k_forward) for use in any backward steps
           lnkf = forward_rate 
           
           !! At this stage, forward rate contains ln(k) for both regular and Lindemann steps
           
           !! Finalise forward rate ===========================================
           !! Multiply up each reactant contrib (in log space)
           do jspec = 1,num_reactants(istep)
              ispec = reactant_list(istep,jspec)  !! ispec is the jspec-th reactant of step istep
              
              forward_rate = forward_rate + nu_dash(istep,ispec)* &
                            (logYspec_local(ispec) + lnro(i) - log(molar_mass(ispec)))
           end do
           !! exponential
           forward_rate = exp(forward_rate)                       

           !! Backward rate ===================================================
           if(gibbs_rate_flag(istep).eq.1) then

              !! Evaluate backwards rate                           
              !! Initialise ln(backwards rate) = ln(forwards rate)
              backward_rate = lnkf
           
              !! Loop over all species in step
              do jspec = 1,num_reactants(istep) + num_products(istep)
                 ispec = stepspecies_list(istep,jspec)
                 
                 call evaluate_gibbs_at_node(tmpT,logT,ispec,gibbs)
                                 
                 backward_rate = backward_rate + delta_nu(istep,ispec)*( &
                                  gibbs)! See Senga?? + logP0ovR0 - logT)                                
              end do          
              !! backward_rate contains k_{b,m}
              
              !! Loop over all products and build log of backward_rate
              do jspec = 1,num_products(istep)
                 ispec = product_list(istep,jspec)
                 
                 backward_rate = backward_rate + nu_ddash(istep,ispec)* &
                                (logYspec_local(ispec) + lnro(i) - log(molar_mass(ispec)))
              
              end do
              !! exponential
              backward_rate = exp(backward_rate)
           else                    
              !! No backward production otherwise       
              backward_rate = zero                               
           end if
                      
                      
           !! Net production rate =============================================
           do jspec = 1,num_reactants(istep) + num_products(istep)
              ispec = stepspecies_list(istep,jspec)
           
              !! Net molar production rate of ispec by step istep
              molar_production_rate = delta_nu(istep,ispec)* &
                                      (forward_rate - backward_rate)* &
                                      third_body_conc
          
              !! Mass production rate
              production_rate = molar_production_rate*molar_mass(ispec)

!if(tmpT.gt.2300.0d0.and.ispec.eq.3) write(6,*) tmpT,istep,ispec,production_rate
              
              !! Add to rate for this species
              rateYspec(ispec) = rateYspec(ispec) + molar_production_rate
!              rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) + production_rate/tmpro
           end do
           
        end do
        !! End steps loop =====================================================

        !! Loop over all species and add rate to rhs
        do ispec = 1,nspec
           rhs_yspec(i,ispec) = rhs_Yspec(i,ispec) + molar_mass(ispec)*rateYspec(ispec)/tmpro           
        end do
alpha_out(i) = rateYspec(3)*molar_mass(3)/tmpro           
     end do
     !$omp end parallel do
     
     deallocate(logYspec_local,Yspec_local,rateYspec)
     
  
     return
  end subroutine calculate_chemical_production_rate  
!! ------------------------------------------------------------------------------------------------  
#endif
end module chemistry

