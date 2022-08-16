module thermodynamics
  !! This module contains routines to calculate various thermodynamic
  !! and temperature dependent transport properties

  use kind_parameters
  use common_parameter
  use common_vars
  implicit none


contains
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_temperature_and_cp
     !! Evaluate temperature, either as constant (isothermal), or analytically (assuming cp=const)
     !! or numerically (assuming cp=cp(T)).
     integer(ikind) :: i
     real(rkind) :: tmp_kinetic,c

     allocate(T(np),cp(np))
     T(:)=zero
     !! Calculating Cp is hard-coded for now
     cp(:) = 1.0045d3
#ifndef isotT

#ifdef pgl       
     !! Perfect gas law. cp is constant, analytic expression.
     !$omp parallel do private(tmp_kinetic)
     do i=1,np
        tmp_kinetic = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        T(i) = (roE(i)/exp(lnro(i)) - tmp_kinetic)*0.4d0/Rgas_mix(i)!/(cp(i)-Rgas_mix(i))
     end do
     !$omp end parallel do
#else
     !! Semi-perfect gas. cp=poly(T). Numerical solution. TBC.
     !$omp parallel do private(tmp_kinetic)
     do i=1,np
        tmp_kinetic = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        T(i) = (roE(i)/exp(lnro(i)) - tmp_kinetic)/(cp(i)-Rgas_mix(i))       
     end do
     !$omp end parallel do
#endif    

     !! Find the maximum sound speed (useful later for time-step constraint)
     cmax = zero
     !$omp parallel do private(c) reduction(max:cmax)
     do i=1,npfb
        c = sqrt(cp(i)*Rgas_mix(i)*T(i)/(cp(i)-Rgas_mix(i)))
        cmax = max(c,cmax)
     end do
     !$omp end parallel do
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
    
     return
  end subroutine evaluate_mixture_gas_constant
!! ------------------------------------------------------------------------------------------------  
  subroutine evaluate_transport_properties
     !! Uses temperature, cp and density to evaluate thermal conductivity, viscosity and 
     !! molecular diffusivity
     integer(ikind) :: ispec,i
     real(rkind) :: tmp
     
     allocate(lambda_th(npfb),visc(npfb),Mdiff(npfb,nspec))
     
     !$omp parallel do private(ispec,tmp)
     do i=1,npfb
     
        !! Thermal conductivity
        lambda_th(i) = cp(i)*Alambda*(T(i)/T0)**rlambda
     
        !! Viscosity
        visc(i) = lambda_th(i)*Pr/cp(i)

  
        !! Molecular diffusivity
        tmp = lambda_th(i)/(exp(lnro(i))*cp(i))
        do ispec=1,nspec
           Mdiff(i,ispec) = tmp/Lewis_number(ispec)
        end do        

visc(i) = visc0
lambda_th(i) = lambda_th0  
Mdiff(i,:) = Mdiff0

!write(6,*) visc(i),lambda_th(i),Mdiff(i,1)        
     end do
     !$omp end parallel do
  
     return
  end subroutine evaluate_transport_properties
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
