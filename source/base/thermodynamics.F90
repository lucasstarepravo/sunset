module thermodynamics
  !! This module contains routines to calculate various thermodynamic
  !! and temperature dependent transport properties

  use kind_parameters
  use common_parameter
  use common_2d
  implicit none


contains
!! ------------------------------------------------------------------------------------------------
  subroutine pressure_from_primary_vars
     integer(ikind) :: i
     real(rkind) :: tmp_kinetic
     
     allocate(p(np))

#ifdef isoT    
     !! Isothermal, calculate pressure from density
     !$omp parallel do 
     do i=1,np    !! N.B. this is over ALL particles incl. ghosts
        p(i) = csq*exp(lnro(i))
     end do
     !$omp end parallel do      
#else
     !! Thermally perfect gas, calculate pressure from energy, density, velocity     
     !$omp parallel do private(tmp_kinetic)
     do i=1,np         !! N.B. this is over ALL particles incl. ghosts
        tmp_kinetic = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))*exp(lnro(i))
        p(i) = (roE(i) - tmp_kinetic)*gammagasm1
     end do
     !$omp end parallel do  
#endif     
  end subroutine pressure_from_primary_vars
!! ------------------------------------------------------------------------------------------------
  subroutine temp_from_primary_vars
     !! Temperature from energy and velocity, assuming perfect gas
     integer(ikind) :: i
     real(rkind) :: tmp_kinetic
     
     allocate(T(np))  
     
     !$omp parallel do private(tmp_kinetic)
     do i=1,np
#ifdef isoT
        T(i) = zero
#else        
        tmp_kinetic = half*(u(i)*u(i) + v(i)*v(i) + w(i)*w(i))
        T(i) = (roE(i)/exp(lnro(i)) - tmp_kinetic)*gammagasm1/Rs0
#endif        
     end do
     !$omp end parallel do
  
  end subroutine temp_from_primary_vars
!! ------------------------------------------------------------------------------------------------
  subroutine visc_from_temp
     !! Viscosity from temperature, according to ... (TBC)
     integer(ikind) :: i
     
     allocate(visc(np))
     
     !$omp parallel do
     do i=1,np
        visc(i) = visc0
     end do
     !$omp end parallel do
  
  
  end subroutine visc_from_temp
!! ------------------------------------------------------------------------------------------------
  function calc_sound_speed(p_local,ro_local) result(c)
     !! Sound speed from p,ro, or prescribed if isoThermal
     real(rkind), intent(in) :: p_local,ro_local
     real(rkind) ::  c
#ifdef isoT
     c = sqrt(csq)
#else  
     c = dsqrt(gammagas*p_local/ro_local)
#endif      
  end function calc_sound_speed
end module thermodynamics
