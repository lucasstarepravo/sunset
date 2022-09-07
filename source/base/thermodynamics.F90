module thermodynamics
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to calculate various thermodynamic
  !! and temperature dependent transport properties
  
  !! Various thermodynamic options ::
  !! 1) isoT      - ISOTHERMAL FLOW. Don't solve an energy equation, p=ro*c*c with c a constant. Many
  !!                arrays are not used, and so not allocated.
  !! 2) not(tdtp) - Transport properties independent of temperature. Still build mixture-values of 
  !!                cp and Rgas from their base molar values, but let visc,lambda_th and Mdiff take
  !!                base values.
  !! 3) tdtp      - Temperature dependent transport properties. cp = cp(T) (a pre-defined polynomial),
  !!                and so calculation of the temperature from the energy involves solution of a 
  !!                non-linear system. The temperature dependence of the viscosity etc is by a power 
  !!                scaling of the base values of (T/T_ref)**r, with r a constant.
  !!
  
  !! Evaluation of and direct reference to CHEMKIN polynomials should only happen within the routines
  !! in this module.
  use kind_parameters
  use common_parameter
  use common_vars
  implicit none


contains
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_temperature_and_cp
     !! Evaluate temperature, either as constant (isothermal), or analytically (assuming cp=const)
     !! or numerically (assuming cp=cp(T)).
     integer(ikind) :: i,ispec,iorder,NRiters,maxiters
     real(rkind) :: tmp_kinetic,cp_tmp,deltaT
     real(rkind) :: T_coef_A,T_tmp,fT,dfT
     real(rkind),dimension(:),allocatable :: T_coef_B
     real(rkind),parameter :: T_tolerance=1.0d-10
     integer(ikind),parameter :: NRiters_max=100
     logical :: keepgoing

#ifndef isoT
     allocate(cp(np))
     T(:)=zero
     allocate(T_coef_B(polyorder_cp))
     
     !! Loop over all nodes
     maxiters = 0
     !$omp parallel do private(T_coef_A,T_coef_B,ispec,iorder,T_tmp,NRiters,fT,dfT,deltaT) &
     !$omp reduction(max:maxiters)
     do i=1,np
        !!Initialise coefficients:
        T_coef_A = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i)) - roE(i)/exp(lnro(i))
        T_coef_B(1) = -Rgas_mix(i)
        T_coef_B(2:polyorder_cp) = zero
        
        !! Build coefficients
        do iorder = 1,polyorder_cp       
           do ispec=1,nspec
              T_coef_B(iorder) = T_coef_B(iorder) + Yspec(i,ispec)*coef_cp(ispec,iorder,1)/dble(iorder)    
           end do       
        end do
        do ispec = 1,nspec
           T_coef_A = T_coef_A + Yspec(i,ispec)*coef_cp(ispec,polyorder_cp+1,1)
        end do
        
        !! Initial guess for T
        T_tmp = T_ref!(i)
        
        !! Newton-Raphson iterations
        keepgoing = .true.
        NRiters = 0
        do while(keepgoing)
           NRiters = NRiters + 1
        
           !! Evaluate f(T) and f'(T)
           fT = T_coef_B(polyorder_cp)*T_tmp
           dfT = dble(polyorder_cp)*T_coef_B(polyorder_cp)
           do iorder = polyorder_cp-1,1,-1
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
         
     !! Calculate  mixture cp
     !$omp parallel do private(ispec,cp_tmp,iorder)
     do i=1,np
        cp(i) = zero
        do ispec=1,nspec
           cp_tmp = coef_cp(ispec,polyorder_cp,1)
           do iorder=polyorder_cp-1,1,-1
              cp_tmp = cp_tmp*T(i) + coef_cp(ispec,iorder,1)
           end do  
        
           cp(i) = cp(i) + Yspec(i,ispec)*cp_tmp
           
        end do
     end do
     !$omp end parallel do

#endif        

     return
  end subroutine evaluate_temperature_and_cp
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_enthalpy(Temp,ispec,enthalpy)  
     !! Evaluate the enthalpy of species ispec based on temperature.
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(in) :: Temp
     real(rkind),intent(out) :: enthalpy
     real(rkind) :: enth_tmp
     integer(ikind) :: iorder
     
     !! Enthalpy = polynomial(T).    
     enth_tmp = Temp*coef_cp(ispec,polyorder_cp,1)/dble(polyorder_cp)
     do iorder=polyorder_cp-1,1,-1
        enth_tmp = Temp*(enth_tmp + coef_cp(ispec,iorder,1)/dble(iorder) )
     end do
     enth_tmp = enth_tmp + coef_cp(ispec,polyorder_cp+1,1)
     enthalpy = enth_tmp   
                  
  
  
     return
  end subroutine evaluate_enthalpy
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_pressure
     integer(ikind) :: i
     real(rkind) :: tmp_kinetic,tmpro
     
     allocate(p(np))

#ifdef isoT    
     !! Isothermal, calculate pressure from density
     !$omp parallel do private(tmpro)
     do i=1,np    !! N.B. this is over ALL particles incl. ghosts
        tmpro = exp(lnro(i))
        p(i) = csq*tmpro
     end do
     !$omp end parallel do      
#else
     !! semi-perfect gas - evaluate from ro,T,and Y
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
     !! molecular diffusivity
     integer(ikind) :: ispec,i
     real(rkind) :: tmp
     
     allocate(visc(npfb),Mdiff(npfb,nspec))
#ifndef isoT     
     allocate(lambda_th(npfb))
     !$omp parallel do private(ispec,tmp)
     do i=1,npfb
     
        !! Viscosity
        visc(i) = visc_ref*(T(i)/T_ref)**r_temp_dependence
     
        !! Thermal conductivity
        lambda_th(i) = cp(i)*visc(i)/Pr
       
        !! Molecular diffusivity
        tmp = lambda_th(i)/(exp(lnro(i))*cp(i))
        do ispec=1,nspec
           Mdiff(i,ispec) = tmp/Lewis_number(ispec)
        end do        

!visc(i) = visc_ref
!lambda_th(i) = lambda_th_ref
!Mdiff(i,:) = Mdiff_ref

!write(6,*) visc(i),lambda_th(i),Mdiff(i,1)     
     end do
     !$omp end parallel do
#else
     visc(:) = visc_ref
     Mdiff(:,:) = Mdiff_ref
#endif     
  
     return
  end subroutine evaluate_transport_properties
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_dcpdT(Temp,ispec,cpispec,dcpdT)
     !! Evaluate the rate of change of cp of species ispec given a temperature, and also the cp
     !! for that species
     real(rkind),intent(in) :: Temp
     integer(ikind),intent(in) :: ispec
     real(rkind),intent(out) :: dcpdT,cpispec
     integer(ikind) :: iorder
     
         
     dcpdT = dble(polyorder_cp-1)*coef_cp(ispec,polyorder_cp,1)
     cpispec = coef_cp(ispec,polyorder_cp,1)
     do iorder=polyorder_cp-1,2,-1
        dcpdT = dcpdT*Temp + dble(iorder-1)*coef_cp(ispec,iorder,1)
        cpispec = cpispec*Temp + coef_cp(ispec,iorder,1)
     end do
     cpispec = cpispec*Temp + coef_cp(ispec,1,1)
  
     return
  end subroutine evaluate_dcpdT
!! ------------------------------------------------------------------------------------------------
  function calc_sound_speed(cp_local,Rgm_local,T_local) result(c)
     !! Sound speed from cp,Rmix and T (or prescribed by csq if isoT)
     real(rkind), intent(in) :: cp_local,Rgm_local,T_local
     real(rkind) ::  c
#ifdef isoT
     c = sqrt(csq)
#else  
     c = dsqrt(cp_local*Rgm_local*T_local/(cp_local-Rgm_local))
#endif      
  end function calc_sound_speed
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_energy
     !! Evaluate roE based on u,lnro,T. This routine is only called at start-up.
     integer(ikind) :: i,ispec
     real(rkind) :: enthalpy
     
     !$omp parallel do private(ispec,enthalpy)
     do i=1,npfb
        !! Initialise roE with K.E. term
        roE(i) = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        
        !! Loop over species
        do ispec=1,nspec       
           !! Evaluate local species enthalpy
           call evaluate_enthalpy(T(i),ispec,enthalpy)
           
           !! Add species contribution to energy
           roE(i) = roE(i) + Yspec(i,ispec)*(enthalpy - Rgas_universal*T(i)/molar_mass(ispec))
        end do           
           
        !! Multiply to get roE
        roE(i) = roE(i)*exp(lnro(i))

     end do
     !$omp end parallel do
       
     return
  end subroutine initialise_energy  
!! ------------------------------------------------------------------------------------------------  
end module thermodynamics
