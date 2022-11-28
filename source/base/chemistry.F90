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
  use omp_lib
  implicit none


contains
#ifdef react
!! ------------------------------------------------------------------------------------------------  
  subroutine calculate_chemical_production_rate
     integer(ikind) :: j,i,ispec,istep,jspec,jstep,ithirdbody,kspec
     real(rkind) :: arrhenius_rate_back,gibbs_tmp,tmpro
     real(rkind) :: mass_production_rate,arrhenius_rate0,arrhenius_rate
     real(rkind) :: p_reduced,net_rate,logroYovW,enthalpy,heat_release
     real(rkind),dimension(:,:),allocatable :: rateYspec,gibbs
     real(rkind),dimension(:),allocatable :: rate,third_body_conc,backward_rate
         
     segment_tstart = omp_get_wtime()                  
     
     allocate(rateYspec(npfb,nspec));rateYspec = zero
     allocate(rate(npfb));rate=zero
     allocate(third_body_conc(npfb));third_body_conc=one
     allocate(backward_rate(npfb));backward_rate = zero
     
     !! Pre-evaluate gibbs functions as required ==============================
     if(num_gibbs_species.ne.0) then
        !! Allocate space
        allocate(gibbs(npfb,num_gibbs_species));gibbs=zero
        
        !! Loop over all species 
        do ispec = 1,nspec
           if(gibbs_flag_species(ispec).ne.0) then !! If this species needs gibbs evaluation
              jspec = gibbs_flag_species(ispec)  !! jspec is the species index within the gibbs-list
              
              !! Loop over all nodes
              !$omp parallel do private(gibbs_tmp)
              do i=1,npfb
                 !! Evaluate gibbs
                 call evaluate_gibbs_at_node(T(i),log(T(i)),ispec,gibbs_tmp)                 
                 
                 !! Store in gibbs array
                 gibbs(i,jspec) = gibbs_tmp
              end do
              !$omp end parallel do
           
           end if
        end do            
     end if
     !! =======================================================================

     !! Loop over all steps ===================================================
     do istep = 1,nsteps       

        !! Forward rate =======================================================

        !! Evaluate Arrhenius rate for this step
        !$omp parallel do private(arrhenius_rate)
        do i=1,npfb        
           arrhenius_rate = arrhenius_coefs(istep,1) &
                          + arrhenius_coefs(istep,2)*log(T(i)) &
                          - arrhenius_coefs(istep,3)/T(i)
 
           !! Loop over all reactants and build log of forward rate
           rate(i) = arrhenius_rate
        end do
        !$omp end parallel do
        !! rate contains ln(k_{f,m})
                    
        !! Third bodies =======================================================
        if(third_body_flag(istep).ne.0) then
           ithirdbody = third_body_flag(istep) !! This is the type of third body for step istep
           
           !! Build third-body concentration
           !$omp parallel do private(ispec)
           do i=1,npfb
              third_body_conc(i) = zero
              do ispec = 1,nspec
                 third_body_conc(i) = third_body_conc(i) &
                                    + third_body_efficiencies(ithirdbody,ispec)* &
                                      Yspec(i,ispec)* &
                                      one_over_molar_mass(ispec)
              end do
              third_body_conc(i) = third_body_conc(i)*ro(i)
           end do
           !$omp end parallel do           
        else
           third_body_conc(:) = one
        end if   

           
        !! Lindemann steps ====================================================
        if(lindemann_form_flag(istep).ne.0) then
           jstep = lindemann_form_flag(istep)  !! This is the jstep-th Lindemann step


           !$omp parallel do private(arrhenius_rate0,p_reduced)
           do i=1,npfb              
              !! Evaluate k0
              arrhenius_rate0 = lindemann_coefs(jstep,1) + &
                                lindemann_coefs(jstep,2)*log(T(i)) - &
                                lindemann_coefs(jstep,3)/T(i)

              !! Reduced pressure
              p_reduced = exp(arrhenius_rate0 - rate(i))*third_body_conc(i)
              
              !! rate constant for this step
              rate(i) = rate(i) + log(p_reduced/(one+p_reduced)) &
                                + lindemann_coefs(jstep,4)       
              
              !! Reset third body here
              third_body_conc(i) = one               
           end do
           !$omp end parallel do
        end if                     
           
        !! Store log(k_forward) for use in any backward steps
        backward_rate(:) = rate(:)

        !! At this stage, forward rate contains ln(k) for both regular and Lindemann steps
                   
        !! Finalise forward rate ==============================================
        !$omp parallel do private(jspec,ispec,logroYovW)
        do i=1,npfb
           !! Multiply up each reactant contrib (in log space)
           do jspec = 1,num_reactants(istep)
              ispec = reactant_list(istep,jspec)  !! ispec is the jspec-th reactant of step istep
              
              logroYovW = log(ro(i)*max(Yspec(i,ispec),verysmall)*one_over_molar_mass(ispec))
              
              rate(i) = rate(i) + nu_dash(istep,ispec)*logroYovW
           end do
           !! exponential 
           rate(i) = exp(rate(i))  
        end do
        !$omp end parallel do

        !! Backward rate ======================================================
        if(gibbs_rate_flag(istep).eq.1) then

           !! Evaluate backwards rate                           
           !$omp parallel do private(jspec,ispec,logroYovW)
           do i=1,npfb
           
              !! Loop over all species in step
              do jspec = 1,num_reactants(istep) + num_products(istep)
                 ispec = stepspecies_list(istep,jspec)
                                 
                 !! index of this species in gibbs list
                 kspec = gibbs_flag_species(ispec)
                                 
                 backward_rate(i) = backward_rate(i) + delta_nu(istep,ispec)*gibbs(i,kspec)
              end do          
              !! backward_rate contains k_{b,m}
              
              !! Loop over all products and build log of backward_rate
              do jspec = 1,num_products(istep)
                 ispec = product_list(istep,jspec)
                 
                 logroYovW = log(ro(i)*max(Yspec(i,ispec),verysmall)*one_over_molar_mass(ispec))
                 
                 backward_rate(i) = backward_rate(i) + nu_ddash(istep,ispec)*logroYovW           
              end do
              !! exponential 
              backward_rate(i) = exp(backward_rate(i))
           end do
           !$omp end parallel do
        else                    
           !! No backward production otherwise       
           backward_rate(:) = zero                               
        end if
                 
        !! Net production rate ================================================
        !$omp parallel do private(net_rate)
        do i=1,npfb
           
           !! net
           net_rate = rate(i) - backward_rate(i)
           
           !! Add third body concs
           rate(i) = net_rate*third_body_conc(i)
        end do
        !$omp end parallel do
        
        
        !! Loop over all species in this step, and add to the total rate for that species
        do jspec = 1,num_reactants(istep) + num_products(istep)
           ispec = stepspecies_list(istep,jspec)        
 
           !$omp parallel do 
           do i=1,npfb                                                        
              !! Add net molar production rate for this species
              rateYspec(i,ispec) = rateYspec(i,ispec) + delta_nu(istep,ispec)*rate(i)
           end do
           !$omp end parallel do
        end do
           
     end do !! End of steps loop ==============================================
         
     !! Add rate onto RHS =====================================================
     !$omp parallel do private(ispec,tmpro,heat_release,mass_production_rate)
     do i=1,npfb       
        tmpro = one/ro(i)    !! N.B. tmpro holds 1/ro here.
     
        !! Zero heat release
        heat_release = zero
     
        !! Loop over all species and add rate to rhs        
        do ispec = 1,nspec
        
           !! Convert from molar to mass production rate
           mass_production_rate = rateYspec(i,ispec)*molar_mass(ispec)
        
           !! Augment the RHS of Yspec for species ispec
           rhs_yspec(i,ispec) = rhs_Yspec(i,ispec) + mass_production_rate*tmpro           
                              
           !! Augment the heat release (production rate x enthalpy of formation)
           heat_release = heat_release - mass_production_rate*coef_h(ispec,polyorder_cp+2)

        end do
alpha_out(i) = heat_release


     end do
     !$omp end parallel do     

     !! Build contribution to source terms for boundary conditions ============
     if(nb.ne.0) then
        allocate(sumoverspecies_homega(nb))
        allocate(reaction_rate_bound(nb,nspec))
        !$omp parallel do private(j,ispec,enthalpy)
        do j=1,nb
           i=boundary_list(j)
              
           sumoverspecies_homega(j) = zero
           !! Loop over species
           do ispec=1,nspec                                    
              !! Store the reaction rate on the boundary
              reaction_rate_bound(j,ispec) = rateYspec(i,ispec)*molar_mass(ispec)                                     
              
              !! Evaluate enthalpy of species ispec
              call evaluate_enthalpy_only_at_node(T(i),ispec,enthalpy)           
           
              !! Evaluate "reduced enthalpy"
              enthalpy = enthalpy - cp(i)*T(i)*Rgas_universal/(Rgas_mix(i)*molar_mass(ispec))
           
              !! Augment sum of reduced_h*omega
              sumoverspecies_homega(j) = sumoverspecies_homega(j) + &
                                      enthalpy*reaction_rate_bound(j,ispec)        
              
           end do
        end do
        !$omp end parallel do
     end if           
     
     !! De-allocation of arrays
     if(num_gibbs_species.ne.0) deallocate(gibbs)
     deallocate(rateYspec,rate,backward_rate,third_body_conc)
  
     !! Profiling
     segment_tend = omp_get_wtime()
     segment_time_local(6) = segment_time_local(6) + segment_tend - segment_tstart  
  
     return
  end subroutine calculate_chemical_production_rate  
!! ------------------------------------------------------------------------------------------------ 
#endif
end module chemistry

