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
  !!                scaling of the base values of (T/T0)**r, with r a constant.
  !!
  use kind_parameters
  use common_parameter
  use common_vars
  implicit none


contains
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_temperature_and_cp
     !! Evaluate temperature, either as constant (isothermal), or analytically (assuming cp=const)
     !! or numerically (assuming cp=cp(T)).
     integer(ikind) :: i,ispec
     real(rkind) :: tmp_kinetic

#ifndef isoT
     allocate(T(np),cp(np))
     T(:)=zero

     !! Calculate cp from molar cp for each species for now. This is how we do it in long run when
     !! not tdtp.
     !$omp parallel do private(ispec)
     do i=1,np
        cp(i) = zero
        do ispec=1,nspec
           cp(i) = cp(i) + Yspec(i,ispec)*cp0_molar(ispec)*Rgas_universal/molar_mass(ispec)
        end do
     end do
     !$omp end parallel do

#ifdef tdtp       
     !! Temperature dependent transport properties. cp(T) = poly(T). TBC.
     !$omp parallel do private(tmp_kinetic)
     do i=1,np
        tmp_kinetic = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        T(i) = (roE(i)/exp(lnro(i)) - tmp_kinetic)/(cp(i)-Rgas_mix(i))       
     end do
     !$omp end parallel do
#else
     !! cp = constant.
     !$omp parallel do private(tmp_kinetic)
     do i=1,np
        tmp_kinetic = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        T(i) = (roE(i)/exp(lnro(i)) - tmp_kinetic)/(cp(i)-Rgas_mix(i))       
     end do
     !$omp end parallel do
#endif    

#endif        

     return
  end subroutine evaluate_temperature_and_cp
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
     
     allocate(lambda_th(npfb),visc(npfb),Mdiff(npfb,nspec))
#ifndef isoT     
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
     lambda_th(:) = lambda_th_ref
     Mdiff(:,:) = Mdiff_ref
#endif     
  
     return
  end subroutine evaluate_transport_properties
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_viscosity_gradient(i,gradT,grad_visc)
     !! Evaluate the viscosity gradient based on the temperature gradient
     !! local to node i
     integer(ikind),intent(in) :: i
     real(rkind),dimension(:),intent(in) :: gradT
     real(rkind),dimension(:),intent(out) :: grad_visc
     
     grad_visc(:) = r_temp_dependence*visc(i)*gradT(:)/T(i)
     
  
     return
  end subroutine evaluate_viscosity_gradient  
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
end module thermodynamics
