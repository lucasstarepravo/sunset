module thermodynamics
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to calculate thermodynamic properties (e.g. p,T, from ro,u,roE,Y)
  !! and temperature dependent transport properties (e.g. visc,lambda,D)
  
  !! Various thermodynamic options ::
  !! 1) isoT      - ISOTHERMAL FLOW. Don't solve an energy equation, p=ro*c*c with c a constant. Many
  !!                arrays are not used, and so not allocated.
  !! 2) not(isoT) - THERMAL FLOW: cp is a polynomial function of T (it can be a polynomial of 
  !!                order 0), and T is obtained from lnro,u,roE,Y via solution of a non-linear 
  !!                equation with a Newton-Raphson method. Within thermal framework, we have two 
  !!                options:
  !!     a) tdtp       - The temperature dependence of the viscosity, thermal conductivity, and
  !!                     molecular diffusivity is by a power scaling of the base values of 
  !!                     (T/T_ref)**r, with r a constant.
  !!     b) not(tdtp)  - Transport properties independent of temperature, although molecular diffusivity
  !!                     and thermal conductivity may be non-uniform, as they are functions of 
  !!                     composition and density.
  
  !! Evaluation of and direct reference to CHEMKIN polynomials (i.e. use of coef_cp,coef_h) should 
  !! only happen within the routines in this module.
  use kind_parameters
  use common_parameter
  use common_vars
  implicit none


contains
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_temperature_and_pressure
     !! Evaluate temperature numerically with cp = polynomial(T), and the pressure from P=ro*R*T
     !! Also evaluate the mixture gas constant, and the mixture specific heat capacity
     integer(ikind) :: i,ii,j,ispec,iorder,NRiters,maxiters,sumiters
     real(rkind) :: tmp_kinetic,deltaT,cp_tmp,tmpro
     real(rkind),dimension(:),allocatable :: fT_coef_C,dfT_coef_C
     real(rkind) :: fT_coef_C0
     real(rkind) :: T_tmp,fT,dfT
     real(rkind),parameter :: T_tolerance=1.0d-10
     integer(ikind),parameter :: NRiters_max=100
     logical :: keepgoing

   

#ifndef isoT
     allocate(fT_coef_C(polyorder_cp+1),dfT_coef_C(polyorder_cp+1))
     
     !! Loop over all nodes
     maxiters = 0;sumiters=0
     !$omp parallel do private(fT_coef_C0,fT_coef_C,ispec,iorder,T_tmp,NRiters,fT,dfT,deltaT,dfT_coef_C, &
     !$omp cp_tmp,tmpro) &
     !$omp reduction(max:maxiters) reduction(+:sumiters)
     do i=1,np
        
        !! For fluid/boundary nodes, and halo nodes, solve non-linear equation
!        if(i.le.npfb.or.i.gt.np_nohalo)then      
        
           !! Evaluate the gas constant for the mixture
           Rgas_mix(i) = zero
           do ispec = 1,nspec
              Rgas_mix(i) = Rgas_mix(i) + Yspec(i,ispec)/molar_mass(ispec)
           end do
           Rgas_mix(i) = Rgas_mix(i)*Rgas_universal     
          
           !! Density from its logarithm  
           tmpro = exp(lnro(i)) 
         
           !!Initialise coefficients:
           fT_coef_C0 = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i)) - roE(i)/tmpro
           fT_coef_C(1) = -Rgas_mix(i)
           fT_coef_C(2:polyorder_cp+1) = zero
        
           !! Build coefficients
           do iorder = 1,polyorder_cp+1
              do ispec=1,nspec
                 fT_coef_C(iorder) = fT_coef_C(iorder) + Yspec(i,ispec)*coef_h(ispec,iorder)
              end do       
           end do
           do ispec = 1,nspec
              fT_coef_C0 = fT_coef_C0 + Yspec(i,ispec)*coef_cp(ispec,polyorder_cp+2)
           end do
        
           !! Make coefficients for dfT
           do iorder = 1,polyorder_cp+1
              dfT_coef_C(iorder) = dble(iorder)*fT_coef_C(iorder)
           end do
        
           !! Initial guess for T is current temperature..
           T_tmp = T(i)
        
           !! Newton-Raphson iterations
           keepgoing = .true.
           NRiters = 0
           do while(keepgoing)
              NRiters = NRiters + 1
        
              !! Evaluate f(T) and f'(T)
              fT = fT_coef_C(polyorder_cp+1)*T_tmp
              dfT = dfT_coef_C(polyorder_cp+1)
              do iorder = polyorder_cp,1,-1
                 fT = (fT + fT_coef_C(iorder))*T_tmp
                 dfT = dfT*T_tmp + dfT_coef_C(iorder)
              end do
              fT = fT + fT_coef_C0
        
!if(iproc.eq.0.and.i.eq.100) then
!     write(6,*) roE(i),NRiters,T_tmp,fT,dfT
!end if        
              !! Calculate new T
              deltaT = - fT/dfT
              T_tmp = T_tmp + deltaT

              !! Check for convergence
              if(abs(deltaT).le.T_tolerance) then
                 keepgoing = .false.
              end if
              if(NRiters.ge.NRiters_max) then
                 keepgoing = .false.
              end if
           
           end do
        
           !! Pass new T back to temperature array
           T(i) = T_tmp
                
           !! Find the maximum number of iterations over this processor and sum of iterations
           maxiters = max(NRiters,maxiters)
           sumiters = sumiters + NRiters
        
           !! Evaluate the specific heat capacity of the mixture
           cp(i) = zero
           do ispec=1,nspec
              cp_tmp = coef_cp(ispec,polyorder_cp+1)
              do iorder=polyorder_cp,1,-1
                 cp_tmp = cp_tmp*T(i) + coef_cp(ispec,iorder)
              end do  
              cp(i) = cp(i) + Yspec(i,ispec)*cp_tmp          
           end do      

           !! Evaluate the pressure        
           p(i) = tmpro*Rgas_mix(i)*T(i)
           
!        endif        
     end do
     !$omp end parallel do
     
     !! For mirrors, copy properties (as opposed to solving non-linear equation)
!     !$omp parallel do private(i)
!#ifdef mp
!     do j=npfb+1,np_nohalo
!#else
!     do j=npfb+1,np
!#endif     
!        i=irelation(j)
!        T(j) = T(i)
!        p(j) = p(i)
!        cp(j) = cp(i)
!        Rgas_mix(j) = Rgas_mix(i)
!     
!     end do
!     !$omp end parallel do


     deallocate(fT_coef_C,dfT_coef_C)
#else
     !! Isothermal, set constant T, and p proportional to density
     T(:) = zero
     !$omp parallel do private(tmpro)
     do i=1,np    !! N.B. this is over ALL nodes.
        tmpro = exp(lnro(i))
        p(i) = csq*tmpro
     end do
     !$omp end parallel do      
     
#endif     
          

     return
  end subroutine evaluate_temperature_and_pressure
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_enthalpy_at_node(Temp,ispec,enthalpy,cpispec,dcpdT)  
     !! Evaluate the enthalpy and cp of species ispec based on temperature at one node. Also
     !! evaluate the rate of change of cp with T.
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(in) :: Temp
     real(rkind),intent(out) :: enthalpy,cpispec,dcpdT
     integer(ikind) :: iorder
     
     !! Enthalpy and cp
     enthalpy = Temp*coef_h(ispec,polyorder_cp+1)    
     cpispec = coef_cp(ispec,polyorder_cp+1)
    
     do iorder=polyorder_cp,1,-1

        enthalpy = Temp*(enthalpy + coef_h(ispec,iorder))   
        cpispec = cpispec*Temp + coef_cp(ispec,iorder)
    
     end do    
     enthalpy = enthalpy + coef_h(ispec,polyorder_cp+2)
    
     !! Rate of change of cp with T
     dcpdT = coef_dcpdT(ispec,polyorder_cp+1)
     
     do iorder=polyorder_cp,2,-1
     
        dcpdT = dcpdT*Temp + coef_dcpdT(ispec,iorder)
        
     end do
                 
     return
  end subroutine evaluate_enthalpy_at_node
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_gibbs_at_node(Temp,logT,ispec,gibbs)
     !! Evaluate the gibbs function of species ispec at a node, given T,logT and ispec
     !! N.B. actually returns molar_gibbs/(R0*T)
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(in) :: Temp,logT
     real(rkind),intent(out) :: gibbs
     integer(ikind) :: iorder
     
     !! Polynomial (in T) terms            
     gibbs = coef_gibbs(ispec,polyorder_cp+1)
     do iorder=polyorder_cp,1,-1
        gibbs = coef_gibbs(ispec,iorder) + gibbs*Temp
     end do
     gibbs = coef_gibbs(ispec,polyorder_cp+2)/Temp &
           - coef_gibbs(ispec,polyorder_cp+3)*logT &
           - gibbs
     
     return
  end subroutine evaluate_gibbs_at_node
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_transport_properties
     !! Uses temperature, cp and density to evaluate thermal conductivity, viscosity and 
     !! molecular diffusivity. For isothermal flows, or if not(tdtp), use reference values.
     integer(ikind) :: ispec,i
     real(rkind) :: tmp   
       
#ifndef isoT     
     !$omp parallel do private(ispec,tmp)
     do i=1,npfb
     
        !! Viscosity
#ifdef tdtp
        visc(i) = visc_ref*(T(i)/T_ref)**r_temp_dependence
#else
        visc(i) = visc_ref
#endif        
     
        !! Thermal conductivity
        lambda_th(i) = cp(i)*visc(i)/Pr

#ifdef ms       
        !! Molecular diffusivity
        tmp = visc(i)/(exp(lnro(i))*Pr)
        do ispec=1,nspec
           Mdiff(i,ispec) = tmp*one_over_Lewis_number(ispec)
        end do        
#endif        
   
     end do
     !$omp end parallel do
#else
     visc(:) = visc_ref
#ifdef ms
     Mdiff(:,:) = Mdiff_ref
#endif     
#endif     
  
     return
  end subroutine evaluate_transport_properties
!! ------------------------------------------------------------------------------------------------
  function evaluate_sound_speed_at_node(cp_local,Rgm_local,T_local) result(c)
     !! Sound speed from cp,Rmix and T (or prescribed by csq if isoT)
     real(rkind), intent(in) :: cp_local,Rgm_local,T_local
     real(rkind) ::  c
#ifdef isoT
     c = sqrt(csq)
#else  
     c = dsqrt(cp_local*Rgm_local*T_local/(cp_local-Rgm_local))
#endif      
  end function evaluate_sound_speed_at_node
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_energy
     !! Evaluate roE based on u,lnro,T, over the whole domain. 
     !! This routine is only called at start-up, and it is because loading temperature is a more
     !! intuitive variable to use for input than roE...
     
     !! Additionally calculate the pressure on outflow boundary nodes, and set P_outflow to the
     !! average of this (it should be uniform along bound)
     integer(ikind) :: i,ispec,j,nsum
     real(rkind) :: enthalpy,Rgas_mix_local,p_local,psum,tmpro,cpispec,dummy_real
     
#ifndef isoT     
     !! Evaluate the energy (roE) and pressure.
     !$omp parallel do private(ispec,enthalpy,Rgas_mix_local,tmpro,cpispec)
     do i=1,npfb
          
        !! Initialise roE with K.E. term
        roE(i) = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        
        !! Evaluate density from its logarithm
        tmpro = exp(lnro(i))

        !! Loop over species
        Rgas_mix_local = zero
        do ispec=1,nspec       
           !! Evaluate local species enthalpy
           call evaluate_enthalpy_at_node(T(i),ispec,enthalpy,cpispec,dummy_real)
           
           !! Add species enthalpy contribution
           roE(i) = roE(i) + Yspec(i,ispec)*enthalpy
           
           !! Build the local mixture gas constant
           Rgas_mix_local = Rgas_mix_local + Yspec(i,ispec)*Rgas_universal/molar_mass(ispec)
        end do           

        !! Evaluate the pressure           
        p(i) = tmpro*Rgas_mix_local*T(i)

        !! Subtract RgasT
        roE(i) = roE(i) - Rgas_mix_local*T(i)     
        
!if(node_type(i).eq.1) write(6,*) "reactants",roE(i),tmpro,p(i)
!if(node_type(i).eq.2) write(6,*) "products",roE(i),tmpro,p(i)
           
        !! Multiply to get roE
        roE(i) = roE(i)*tmpro
        
     end do
     !$omp end parallel do
          

     !! Pressure on outflow nodes 
     psum = zero;nsum = 0
     if(nb.ne.0)then
        !$omp parallel do private(i) reduction(+:psum,nsum)
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.2) then
                          
              !! Augment the accumulators for sum of pressure and # outflow nodes
              psum = psum + p(i)
              nsum = nsum + 1
           end if
        end do
        !$omp end parallel do
        
        !! Find the average for processors with outflows. Set to zero otherwise
        if(nsum.ne.0) then
           P_outflow = psum/dble(nsum)
        else
           P_outflow = zero
        end if
        
     end if
#else
     P_outflow = csq*rho_char
#endif

       
     return
  end subroutine initialise_energy  
!! ------------------------------------------------------------------------------------------------  
end module thermodynamics
