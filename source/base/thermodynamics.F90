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
  
  !! Evaluation of and direct reference to CHEMKIN polynomials should only happen within the routines
  !! in this module.
  use kind_parameters
  use common_parameter
  use common_vars
  implicit none


contains
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_temperature
     !! Evaluate temperature, either as constant (isothermal), or analytically (assuming cp=const)
     !! or numerically (assuming cp=cp(T)).
     integer(ikind) :: i,ispec,iorder,NRiters,maxiters
     real(rkind) :: tmp_kinetic,deltaT
     real(rkind) :: T_coef_A,T_tmp,fT,dfT
     real(rkind),dimension(:),allocatable :: T_coef_B
     real(rkind),parameter :: T_tolerance=1.0d-10
     integer(ikind),parameter :: NRiters_max=100
     logical :: keepgoing

#ifndef isoT
     T(:)=zero
     allocate(T_coef_B(polyorder_cp+1))
     
     !! Loop over all nodes
     maxiters = 0
     !$omp parallel do private(T_coef_A,T_coef_B,ispec,iorder,T_tmp,NRiters,fT,dfT,deltaT) &
     !$omp reduction(max:maxiters)
     do i=1,np
        !!Initialise coefficients:
        T_coef_A = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i)) - roE(i)/exp(lnro(i))
        T_coef_B(1) = -Rgas_mix(i)
        T_coef_B(2:polyorder_cp+1) = zero
        
        !! Build coefficients
        do iorder = 1,polyorder_cp+1
           do ispec=1,nspec
              T_coef_B(iorder) = T_coef_B(iorder) + Yspec(i,ispec)*coef_cp(ispec,iorder)/dble(iorder)    
           end do       
        end do
        do ispec = 1,nspec
           T_coef_A = T_coef_A + Yspec(i,ispec)*coef_cp(ispec,polyorder_cp+2)
        end do
        
        !! Initial guess for T
        T_tmp = T(i)
        
        !! Newton-Raphson iterations
        keepgoing = .true.
        NRiters = 0
        do while(keepgoing)
           NRiters = NRiters + 1
        
           !! Evaluate f(T) and f'(T)
           fT = T_coef_B(polyorder_cp+1)*T_tmp
           dfT = dble(polyorder_cp+1)*T_coef_B(polyorder_cp+1)
           do iorder = polyorder_cp,1,-1
              fT = (fT + T_coef_B(iorder))*T_tmp
              dfT = dfT*T_tmp + T_coef_B(iorder)
           end do
           fT = fT + T_coef_A
        
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
        
        !! Find the maximum number of iterations over this processor
        maxiters = max(NRiters,maxiters)
     end do
     !$omp end parallel do
#endif        

!write(6,*) maxiters

     return
  end subroutine evaluate_temperature
!! ------------------------------------------------------------------------------------------------  
  subroutine evaluate_cp
     !! Evaluate the specific heat capacity for the mixture.
     integer(ikind) :: i,ispec,iorder
     real(rkind) :: cp_tmp

#ifndef isoT
     allocate(cp(np))
            
     !! Calculate  mixture cp
     !$omp parallel do private(ispec,cp_tmp,iorder)
     do i=1,np
        cp(i) = zero
        do ispec=1,nspec
           cp_tmp = coef_cp(ispec,polyorder_cp+1)
           do iorder=polyorder_cp,1,-1
              cp_tmp = cp_tmp*T(i) + coef_cp(ispec,iorder)
           end do  
        
           cp(i) = cp(i) + Yspec(i,ispec)*cp_tmp
           
        end do
     end do
     !$omp end parallel do
#endif        

     return
  end subroutine evaluate_cp
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_enthalpy_at_node(Temp,ispec,enthalpy)  
     !! Evaluate the enthalpy of species ispec based on temperature at one node: take in temperature
     !! and species flag, and return single value of enthalpy.
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(in) :: Temp
     real(rkind),intent(out) :: enthalpy
     integer(ikind) :: iorder
     
     !! Enthalpy = polynomial(T).    
     enthalpy = Temp*coef_cp(ispec,polyorder_cp+1)/dble(polyorder_cp+1)
     do iorder=polyorder_cp,1,-1
        enthalpy = Temp*(enthalpy + coef_cp(ispec,iorder)/dble(iorder) )
     end do
     enthalpy = enthalpy + coef_cp(ispec,polyorder_cp+2)
     
     !! Expensive form:
!     enthalpy = coef_cp(ispec,polyorder_cp+2)
!     do iorder=1,polyorder_cp+1
!        enthalpy = enthalpy + (coef_cp(ispec,iorder)/dble(iorder))*Temp**dble(iorder)
!     end do
                 
     return
  end subroutine evaluate_enthalpy_at_node
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_pressure
     !! Evaluate the pressure over the whole domain. If isoT, then based on sound speed and 
     !! density. If thermal, then based on R=ro*R*T.
     integer(ikind) :: i
     real(rkind) :: tmp_kinetic,tmpro
     
#ifdef isoT    
     !! Isothermal, calculate pressure from density
     !$omp parallel do private(tmpro)
     do i=1,np    !! N.B. this is over ALL particles incl. ghosts
        tmpro = exp(lnro(i))
        p(i) = csq*tmpro
     end do
     !$omp end parallel do      
#else
     !! semi-perfect gas - evaluate from ro,T (and indirectly Y, as Rgas_mix is function of Y)
     !$omp parallel do private(tmp_kinetic,tmpro)
     do i=1,np         !! N.B. this is over ALL particles incl. ghosts
        tmpro = exp(lnro(i))
        p(i) = tmpro*Rgas_mix(i)*T(i)
     end do
     !$omp end parallel do  
#endif     

  end subroutine evaluate_pressure
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_mixture_gas_constant
     !! Evaluate the mixture gas constant based on the composition and molar mass of each species.
     !! For isothermal flows, do nothing.
     integer(ikind) :: i,ispec
#ifndef isoT  
     allocate(Rgas_mix(np))
     
     !$omp parallel do private(ispec)
     do i=1,np
        Rgas_mix(i) = zero
        do ispec = 1,nspec
           Rgas_mix(i) = Rgas_mix(i) + Yspec(i,ispec)/molar_mass(ispec)
        end do
        Rgas_mix(i) = Rgas_mix(i)*Rgas_universal
     end do
     !$omp end parallel do    
#endif    
     return
  end subroutine evaluate_mixture_gas_constant
!! ------------------------------------------------------------------------------------------------  
  subroutine evaluate_transport_properties
     !! Uses temperature, cp and density to evaluate thermal conductivity, viscosity and 
     !! molecular diffusivity. For isothermal flows, or if not(tdtp), use reference values.
     integer(ikind) :: ispec,i
     real(rkind) :: tmp
     
     allocate(visc(npfb),Mdiff(npfb,nspec))
#ifndef isoT     
     allocate(lambda_th(npfb))
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
       
        !! Molecular diffusivity
        tmp = lambda_th(i)/(exp(lnro(i))*cp(i))
        do ispec=1,nspec
           Mdiff(i,ispec) = tmp/Lewis_number(ispec)
        end do        
   
     end do
     !$omp end parallel do
#else
     visc(:) = visc_ref
     Mdiff(:,:) = Mdiff_ref
#endif     
  
     return
  end subroutine evaluate_transport_properties
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_dcpdT_at_node(Temp,ispec,cpispec,dcpdT)
     !! Evaluate the rate of change of cp of species ispec given a temperature, and also the cp
     !! for that species, at a single node.
     real(rkind),intent(in) :: Temp
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(out) :: dcpdT,cpispec
     integer(ikind) :: iorder
              
     dcpdT = dble(polyorder_cp)*coef_cp(ispec,polyorder_cp+1)
     cpispec = coef_cp(ispec,polyorder_cp+1)
     do iorder=polyorder_cp,2,-1
        dcpdT = dcpdT*Temp + dble(iorder-1)*coef_cp(ispec,iorder)
        cpispec = cpispec*Temp + coef_cp(ispec,iorder)
     end do
     cpispec = cpispec*Temp + coef_cp(ispec,1)
  
     return
  end subroutine evaluate_dcpdT_at_node
!! ------------------------------------------------------------------------------------------------
  function calc_sound_speed_at_node(cp_local,Rgm_local,T_local) result(c)
     !! Sound speed from cp,Rmix and T (or prescribed by csq if isoT)
     real(rkind), intent(in) :: cp_local,Rgm_local,T_local
     real(rkind) ::  c
#ifdef isoT
     c = sqrt(csq)
#else  
     c = dsqrt(cp_local*Rgm_local*T_local/(cp_local-Rgm_local))
#endif      
  end function calc_sound_speed_at_node
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_energy
     !! Evaluate roE based on u,lnro,T, over the whole domain. 
     !! This routine is only called at start-up, and it is because loading temperature is a more
     !! intuitive variable to use for input than roE...
     
     !! Additionally calculate the pressure on outflow boundary nodes, and set P_outflow to the
     !! average of this (it should be uniform along bound)
     integer(ikind) :: i,ispec,j,nsum
     real(rkind) :: enthalpy,Rgas_mix_local,p_local,psum,tmpro
     
#ifndef isoT     
     !! Evaluate the energy (roE) and pressure.
     !$omp parallel do private(ispec,enthalpy,Rgas_mix_local,tmpro)
     do i=1,npfb
        !! Initialise roE with K.E. term
        roE(i) = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        
        !! Evaluate density from its logarithm
        tmpro = exp(lnro(i))

        !! Loop over species
        Rgas_mix_local = zero
        do ispec=1,nspec       
           !! Evaluate local species enthalpy
           call evaluate_enthalpy_at_node(T(i),ispec,enthalpy)
           
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
        
write(499,*) rp(i,1),roE(i),tmpro

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
