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
     end do
     !$omp end parallel do
     
  
     return
  end subroutine calculate_chemical_production_rates
!! ------------------------------------------------------------------------------------------------  
#endif
end module chemistry

