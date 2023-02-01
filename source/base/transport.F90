module transport
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !! JRCK               |Jan 2023         |New module to implement mix av transport
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines to calculate transport properties.
  
  !! Various options ::
  !!        1) isoT         - isothermal flows, viscosity and molecular diffusivity are set 
  !!                              as constant based on reference values.
  !!        2) not(isoT)    - thermal flows, viscosity, molecular diffusivity and thermal
  !!                          conductivity are dependent on temperature and composition, and
  !!                 a) 
  !!
  !!
  !!
  !!

  use kind_parameters
  use common_parameter
  use common_vars
  implicit none


contains
!! ------------------------------------------------------------------------------------------------
  subroutine load_transport_file
     !! Subroutine loads the transport file transport.in and prepares the coefficients for mixture
     !! average transport properties
     integer(ikind) :: i,j,ispec,jspec
     real(rkind) :: store1
     
     open(unit=13,file='transport.in')     
     
     !! Read header
     read(13,*)
     read(13,*)
     
     !! Reference temperature
     read(13,*)
     read(13,*) T_ref
     read(13,*)
     
     !! Reference pressure
     read(13,*)
     read(13,*) p_ref
     read(13,*)
     
     !! Polynomial coefficients for viscosity
     read(13,*)
     allocate(mxav_coef_visc(nspec,4))
     do ispec=1,nspec
        read(13,*) i,mxav_coef_visc(ispec,1),mxav_coef_visc(ispec,2),mxav_coef_visc(ispec,3),mxav_coef_visc(ispec,4)
     end do
     read(13,*)
     
     !! Precompute elements of combination rule for viscosity
     allocate(mxav_visc_combo1(nspec,nspec),mxav_visc_combo2(nspec,nspec))
     do ispec = 1,nspec
        do jspec = 1,nspec
           mxav_visc_combo1(ispec,jspec) = (molar_mass(jspec)/molar_mass(ispec))**(one/4.0d0)
           mxav_visc_combo2(ispec,jspec) = (one/sqrt(8.0d0))*(one + molar_mass(ispec)/molar_mass(jspec))**(-half)           
        end do
     end do     
          
     !! Polynomial coefficients for thermal conductivity
     read(13,*)
     allocate(mxav_coef_lambda(nspec,4))
     do ispec=1,nspec
        read(13,*) i,mxav_coef_lambda(ispec,1),mxav_coef_lambda(ispec,2), &
                     mxav_coef_lambda(ispec,3),mxav_coef_lambda(ispec,4)
     end do
     read(13,*)
     
     !! Polynomial coefficients for molecular diffusivity
     read(13,*)
     allocate(mxav_coef_Diff(nspec,nspec,4));mxav_coef_Diff=zero
     do ispec=1,nspec
        do jspec=1,ispec
           read(13,*) i,j,mxav_coef_Diff(ispec,jspec,1),mxav_coef_Diff(ispec,jspec,2), &
                          mxav_coef_Diff(ispec,jspec,3),mxav_coef_Diff(ispec,jspec,4)
        end do
     end do
     read(13,*)
     
     !! Copy triangle matrix to symmetric matrix
     do ispec=1,nspec
        do jspec=1,ispec
           mxav_coef_Diff(jspec,ispec,1) = mxav_coef_Diff(ispec,jspec,1)
           mxav_coef_Diff(jspec,ispec,2) = mxav_coef_Diff(ispec,jspec,2)
           mxav_coef_Diff(jspec,ispec,3) = mxav_coef_Diff(ispec,jspec,3)
           mxav_coef_Diff(jspec,ispec,4) = mxav_coef_Diff(ispec,jspec,4)                                 
        end do
     end do
           
     
    
  
     return
  end subroutine load_transport_file
!! ------------------------------------------------------------------------------------------------
  subroutine evaluate_transport_properties
     !! Uses temperature, cp and density to evaluate thermal conductivity, viscosity and 
     !! molecular diffusivity. For isothermal flows, or if not(tdtp), use reference values.
     integer(ikind) :: ispec,i
     real(rkind) :: tmp   
       
#ifndef isoT     
     if(mix_av_flag.eq.0) then
        !! Constant Lewis numbers, with option for temp dependence
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
           !! Molecular diffusivity - actually returning ro*Mdiff
           tmp = visc(i)/Pr
           do ispec=1,nspec
              roMdiff(i,ispec) = tmp*one_over_Lewis_number(ispec)
           end do        
#endif        
   
        end do
        !$omp end parallel do
     else
        !! For mixture averaged transport, there is no (or very complex) analytic relation
        !! between transport gradients and temperature gradients, so we evaluate the transport
        !! properties at halos and mirrors too!
        !$omp parallel do 
        do i=1,np
     
           call evaluate_mixture_average_transport_at_node(i)
 
   
        end do
        !$omp end parallel do     
     end if
#else
     visc(:) = visc_ref
#ifdef ms
     roMdiff(:,:) = ro(i)*Mdiff_ref
#endif     
#endif     
  
     return
  end subroutine evaluate_transport_properties
!! ------------------------------------------------------------------------------------------------  
  subroutine evaluate_mixture_average_transport_at_node(inode)
     !! Evaluates the transport properties at an individual node based on mixture averaged
     !! rules for individual species transport properties.
     integer(ikind),intent(in) :: inode
     integer(ikind) :: ispec,icoef,jspec  
     real(rkind) :: logT,Tmp,logvisc,tmpro,loglambda,molar_mass_mix,logdiff,p_ratio
     real(rkind),dimension(:),allocatable :: species_visc,Xspec,species_lambda
     real(rkind),dimension(:,:),allocatable :: species_diff
     real(rkind) :: store1,store2,store3,store4
     
     !! Allocate space for viscosity and other properties across all species
     allocate(species_visc(nspec),species_lambda(nspec),species_diff(nspec,nspec))
     allocate(Xspec(nspec))
     
     !! Density
     tmpro = ro(inode)
              
     !! Temperature and its logarithm
     Tmp = T(inode)
     logT = log(Tmp/T_ref)
     
     !! Pressure ratio (relative to reference)
     p_ratio = p_ref/p(inode)
     
     !! Loop over all species and evaluate viscosity, thermal conductivity, mixture mol mass, and Xspec
     molar_mass_mix = zero
     do ispec=1,nspec
         
        !! Evaluate species viscosity
        logvisc = mxav_coef_visc(ispec,4)
        logvisc = logvisc*logT + mxav_coef_visc(ispec,3)
        logvisc = logvisc*logT + mxav_coef_visc(ispec,2)
        logvisc = logvisc*logT + mxav_coef_visc(ispec,1)        
        species_visc(ispec) = exp(logvisc)
        
        !! Evaluate species thermal conductivity
        loglambda = mxav_coef_lambda(ispec,4)
        loglambda = loglambda*logT + mxav_coef_lambda(ispec,3)
        loglambda = loglambda*logT + mxav_coef_lambda(ispec,2)
        loglambda = loglambda*logT + mxav_coef_lambda(ispec,1)                
        species_lambda(ispec) = exp(loglambda)

        !! Build Xspec and mixture molar mass
        Xspec(ispec) = Yspec(inode,ispec)*one_over_molar_mass(ispec)/tmpro
        molar_mass_mix = molar_mass_mix + Xspec(ispec)
        
        !! Evaluate pairwise species diffusion (evaluate lower triangle, and copy to upper)
        do jspec = 1,ispec
           logdiff = mxav_coef_diff(ispec,jspec,4)
           logdiff = logdiff*logT + mxav_coef_diff(ispec,jspec,3)
           logdiff = logdiff*logT + mxav_coef_diff(ispec,jspec,2)
           logdiff = logdiff*logT + mxav_coef_diff(ispec,jspec,1)           
           !! Include pressure scaling here
           species_diff(ispec,jspec) = exp(logdiff)*p_ratio                                      
           species_diff(jspec,ispec) = species_diff(ispec,jspec) !! Copy to upper triangle
        end do
        
     end do
     
     !! Finalise mixture molar mass and Xspec
     molar_mass_mix = one/molar_mass_mix
     Xspec(:) = Xspec(:)*molar_mass_mix
         
     
     !! Combination rule for viscosity
if(.false.)then     
     !! Combination rule used by Senga+, Hamish, based on Wilkes 1968
     store1=zero
     do ispec=1,nspec
        store2=zero
        do jspec=1,nspec
        
           !! Evaluate PHI_ispec,jspec
           store3 = sqrt(species_visc(ispec)/species_visc(jspec))
           store3 = one + store3*mxav_visc_combo1(ispec,jspec)
           store3 = mxav_visc_combo2(ispec,jspec)*store3*store3 
        
           !! Sum molar concentration times PHI_ispec,jspec
           store2 = store2 + Xspec(jspec)*store3                   
        end do
        
        !! Species_visc divided by lower sum(over jspec)  
        store3 = species_visc(ispec)/store2
        
        !! Upper sum
        store1 = store1 + Xspec(ispec)*store3    
        
     end do
     visc(inode) = store1    
else
     !! Combination rule used by some versions of PeleC/PeleLM 
     store1=zero 
     do ispec=1,nspec
        store1 = store1 + Xspec(ispec)*species_visc(ispec)**6.0d0
     end do
     visc(inode) = store1**(one/6.0d0)
endif     
     
     !! Combination rule for thermal conductivity
     store1=zero 
     do ispec=1,nspec
        store1 = store1 + Xspec(ispec)*species_lambda(ispec)**(one/4.0d0)
     end do
     lambda_th(inode) = store1**(4.0d0)     

     !! Combination rule for molecular diffusivity
     do ispec = 1,nspec
        !! Initialise with the ispec-ispec contribution subtracted
        store3 = Xspec(ispec) + quitesmall
        store1 = -store3*molar_mass(ispec) 
        store2 = -store3/species_diff(ispec,ispec)
        do jspec = 1,nspec
           store3 = Xspec(jspec) + quitesmall
           store1 = store1 + store3*molar_mass(jspec)
           store2 = store2 + store3/species_diff(ispec,jspec)
        end do
        store1 = store1/(molar_mass_mix*store2)        
        roMdiff(inode,ispec) = store1*tmpro     
     end do
     
     deallocate(Xspec,species_visc,species_lambda,species_diff)
         
     
     !! Temporary non-mix-averaged forms =================================================
     
!     visc(inode) = visc_ref*(T(inode)/T_ref)**r_temp_dependence
    
!     write(6,*) inode,Tmp,lambda_th(inode),cp(inode)*visc(inode)/Pr
    
     !! Thermal conductivity
!     lambda_th(inode) = cp(inode)*visc(inode)/Pr

     !! Molecular diffusivity - actually returning ro*Mdiff
!     store1 = visc(inode)/Pr
!     do ispec=1,nspec
!!        write(6,*) Tmp,ispec,roMdiff(inode,ispec),store1*one_over_Lewis_number(ispec)
!        roMdiff(inode,ispec) = store1*one_over_Lewis_number(ispec)
!     end do        
  
  
     return
  end subroutine evaluate_mixture_average_transport_at_node
end module transport
