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
     integer(ikind) :: i,ispec,istep,jspec
     real(rkind) :: tmpro,production_rate,arrhenius_rate,tmpT,logT
     real(rkind) :: arrhenius_rate_back,gibbs,logP0ovR0
     real(rkind) :: molar_production_rate,product_of_reactants,kforward_prod
     real(rkind) :: product_of_products,kbackward_prod
     real(rkind),dimension(:),allocatable :: logYspec_local
     
     !! Part of gibbs term (saves time doing outside loop)
     logP0ovR0 = log(P_ref/Rgas_universal)
     
     !! space of log(Y)
     allocate(logYspec_local(nspec));logYspec_local = zero
     
     !! Loop over all nodes
     !$omp parallel do private(istep,ispec,jspec,tmpro,tmpT,logT,production_rate,arrhenius_rate, &
     !$omp molar_production_rate,product_of_reactants,arrhenius_rate_back,gibbs,kforward_prod, &
     !$omp product_of_products,kbackward_prod,logYspec_local)
     do i=1,npfb
     
        !! Density, temperature
        tmpro = exp(lnro(i))        
        tmpT = T(i)
        logT = log(tmpT)
        
        !! Store logarithm of species mass fraction for each species
        do ispec = 1,nspec
           logYspec_local(ispec) = log(max(Yspec(i,ispec),verysmall))
        end do
        
        !! Loop over all steps
        do istep = 1,nsteps

           !! Evaluate Arrhenius rate for this step
           arrhenius_rate = arrhenius_coefs(istep,1) &
                          + arrhenius_coefs(istep,2)*logT &
                          - arrhenius_coefs(istep,3)/tmpT

           !! Loop over all reactants and build log of product_of_reactants
           product_of_reactants = zero
           do jspec = 1,num_reactants(istep)
              ispec = reactant_list(istep,jspec)  !! ispec is the jspec-th reactant of step istep
              
              product_of_reactants = product_of_reactants + nu_dash(istep,ispec)* &
                                     (logYspec_local(ispec) + lnro(i) - log(molar_mass(ispec)))
           end do
           !! forward rate * product_of_reactants
           kforward_prod = exp(arrhenius_rate + product_of_reactants)


           !! Backwards reactions
           if(gibbs_rate_flag(istep).eq.1) then

              !! Evaluate backwards rate                           
              !! Initialise ln(backwards rate) = ln(forwards rate)
              arrhenius_rate_back = arrhenius_rate
           
              !! Loop over all species in step
              do jspec = 1,num_reactants(istep) + num_products(istep)
                 ispec = stepspecies_list(istep,jspec)
                 
                 call evaluate_gibbs_at_node(tmpT,logT,ispec,gibbs)
                 
                 arrhenius_rate_back = arrhenius_rate_back + delta_nu(istep,ispec)*( &
                                       gibbs + logP0ovR0 - logT)                                
              end do          
              
              !! Loop over all products and build log of product_of_products
              product_of_products = zero
              do jspec = 1,num_products(istep)
                 ispec = product_list(istep,jspec)
                 
                 product_of_products = product_of_products + nu_ddash(istep,ispec)* &
                                       (logYspec_local(ispec) + lnro(i) - log(molar_mass(ispec)))
              
              end do
              !! backward rate * product_of_products
              kbackward_prod = exp(arrhenius_rate_back + product_of_products)     
           else                    
              !! No backward production otherwise       
              kbackward_prod = zero                               
           end if
           
           !! Loop over all species in step and calculate production rate of ispec
           do jspec = 1,num_reactants(istep) + num_products(istep)
              ispec = stepspecies_list(istep,jspec)
           
              !! Molar production rate   
              molar_production_rate = delta_nu(istep,ispec)*(kforward_prod - kbackward_prod)
           
              !! Mass production rate
              production_rate = molar_production_rate*molar_mass(ispec)

alpha_out(i) = -molar_mass(ispec)*delta_nu(istep,ispec)*kbackward_prod
              
              !! Add to to right hand side
              rhs_Yspec(i,ispec) = rhs_Yspec(i,ispec) + production_rate/tmpro           
           end do
           
        end do

     end do
     !$omp end parallel do
     
     deallocate(logYspec_local)
     
  
     return
  end subroutine calculate_chemical_production_rate  
!! ------------------------------------------------------------------------------------------------  
#endif
end module chemistry

