module setup_flow
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines read in control data, chemistry, and initial conditions
  use kind_parameters
  use common_parameter
  use common_vars
#ifdef mp
  use mpi_transfers
#endif    
  implicit none
  
contains
!! ------------------------------------------------------------------------------------------------
  subroutine initial_solution
     use mirror_boundaries
     use derivatives
     use thermodynamics
     !! Temporary subroutine whilst developing. Initialises all fields
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro     
     
     !! Allocate arrays for properties - primary
     allocate(u(np),v(np),w(np),lnro(np),roE(np),divvel(np))
     allocate(Yspec(np,nspec))
     u=zero;v=zero;w=zero;lnro=zero;roE=one;Yspec=one;divvel=zero

     !! Secondary properties
     allocate(T(np));T=T_ref
     allocate(p(np));p=zero
     allocate(alpha_out(np));alpha_out = zero     
     
     !! Transport properties
     allocate(visc(npfb));visc = visc_ref
#ifndef isoT
     allocate(lambda_th(npfb))
#endif
#ifdef ms     
     allocate(Mdiff(npfb,nspec))
#endif     
     allocate(cp(np),Rgas_mix(np))
     
     !! Allocate the boundary temperatures
     if(nb.ne.0) then
        allocate(T_bound(nb));T_bound = T_ref                 
     end if
     
     !! =======================================================================
     !! Choose initial conditions
#ifndef restart     
#ifdef react
     if(nsteps.eq.1) then
!        call make_1d_1step_flame
     else if(nsteps.eq.21) then
!        call make_1d_21step_flame
     else if(nsteps.eq.35) then
!        call make_1d_25step_flame
     end if
     call load_flame_file    
#else
     call hardcode_initial_conditions     
#endif
#else    
     !! RESTART OPTION. Ask for input number (just hard-coded for now...)
     call load_restart_file(2)
#endif
     !! =======================================================================
            
     !! Set energy from lnro,u,Y,T
     call initialise_energy         
   
     !! Mirrors and halos                        
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchanges_all
#endif          
     
     !! Set the initial velocity divergence
     call calc_divergence(u,v,w,divvel(1:npfb))
     
     !! Mirrors and halos for divvel                   
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchange_divvel
#endif         
     
     !! Set the initial forcing to zero. It will only be changed if PID controller in velcheck is used.
     driving_force(:) = zero

     !! Profiling - re-zero time accumulators
     segment_time_local = zero
     cputimecheck = zero         
     
     !! Initialise the time-stepping (necessary for PID controlled stepping)
     dt = 1.0d-9             
         
     return
  end subroutine initial_solution   
!! ------------------------------------------------------------------------------------------------
  subroutine load_control_data_LUonly
     integer(ikind) :: dummy_int
     real(rkind) :: dummy_real
     
     !! Load data from the control file
     open(unit=12,file='control.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)

     !! Length-scale
     read(12,*)
     read(12,*) L_char
     read(12,*)
     
     !! Velocity-scale
     read(12,*) 
     read(12,*) U_char
     read(12,*)
     
     !! Set inflow velocity, Z-length-scale and characteristic time-scale
     u_inflow = u_char
     Lz = L_char
     Time_char = L_char/u_char
     
     close(12)
          
     return
  end subroutine load_control_data_LUonly
!! ------------------------------------------------------------------------------------------------  
  subroutine load_control_data_all
     integer(ikind) :: dummy_int
     real(rkind) :: dummy_real
     
     !! Load data from the control file
     open(unit=12,file='control.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)

     !! Length-scale
     read(12,*)
     read(12,*) dummy_real
     read(12,*)
     
     !! Velocity-scale
     read(12,*) 
     read(12,*) dummy_real
     read(12,*)
     
     !! Start and end time (in multiples of characteristic time)
     read(12,*)
     read(12,*) time,time_end
     read(12,*)
     time = time*Time_char;time_end = time_end*Time_char
     itime = 0
     
     !! Output frequency (in multiples of characteristic time)
     read(12,*)
     read(12,*) dt_out
     read(12,*)
     dt_out = dt_out*Time_char
        
     !! Gravity?
     read(12,*)
     read(12,*) grav(:)
     read(12,*)
     
     !! Characteristic density
     read(12,*)
     read(12,*) rho_char
     read(12,*)
     
     !! Reference temp
     read(12,*)
     read(12,*) T_ref
     read(12,*) 
     
     !! Reference viscosity
     read(12,*)
     read(12,*) visc_ref
     read(12,*) 
     
     !! Reference pressure
     read(12,*)
     read(12,*) p_ref
     read(12,*) 
         
     !! Prandtl number
     read(12,*)
     read(12,*) Pr
     read(12,*) 
     
     !! Mach number (only used for isothermal flows)
     read(12,*)
     read(12,*) Ma
     read(12,*) 
     
     !! Set the Reynolds number (never used, just nice to calculate it)
     Re = rho_char*u_char*L_char/visc_ref

     !! Set a reference molecular diffusivity
     Mdiff_ref = visc_ref/rho_char/Pr/one   !! the 1 represents Lewis number
     
     !! set the sound speed squared
#ifdef isoT
     csq = (u_char/Ma)**two
#endif      

     !! Read in T-exponent for TDTP if required
     read(12,*)
     read(12,*) r_temp_dependence
     read(12,*)     
#ifndef tdtp
     r_temp_dependence = zero  !! zero it if not required
#endif 
     
     close(12)
          
     return
  end subroutine load_control_data_all
!! ------------------------------------------------------------------------------------------------  
  subroutine load_chemistry_data
     integer(ikind) :: ispec,iorder,istep,dummy_int,jspec
     real(rkind) :: dummy_real
     
     !! Load data from the thermochemistry control file
     open(unit=12,file='thermochem.in')
     read(12,*)                       !! Ignore header and blank line
     read(12,*)
     
     !! Number of chemical species
     read(12,*)
     read(12,*) nspec
     read(12,*)
#ifndef ms
     nspec = 1  !! If not multispecies, nspec must be 1
#endif          

     !! This section reads in transport data for all species.
     
     !! Number of coefs and order of polynomial for cp(T)
     read(12,*)
     read(12,*) ncoefs_cp
     read(12,*)
     !! ncoefs = (polyorder + 1) + 1 + 1: terms in polynomial,coef for h, coef for s
     polyorder_cp = ncoefs_cp - 3
     
     !! Temperature limits
     read(12,*)
     read(12,*) T_low,T_high
     read(12,*)

     !! Allocate space for molar mass, Lewis number, and polynomial fitting for cp(T)   
     allocate(molar_mass(nspec),one_over_Lewis_number(nspec),one_over_molar_mass(nspec))
     allocate(coef_cp(nspec,ncoefs_cp),coef_h(nspec,ncoefs_cp),coef_dcpdT(nspec,ncoefs_cp))
     allocate(coef_gibbs(nspec,ncoefs_cp))
          
     !! Load molar mass, Lewis, and polynomial fits   
     read(12,*) !! Read comment line  
     do ispec = 1,nspec
        read(12,*) !! Read species identifier comment line
        read(12,*) molar_mass(ispec)
        read(12,*) one_over_Lewis_number(ispec)
        one_over_Lewis_number(ispec) = one/one_over_Lewis_number(ispec)
        one_over_molar_mass(ispec) = one/molar_mass(ispec)

        read(12,*)  !! Comment line
        do iorder=1,ncoefs_cp
           read(12,*) coef_cp(ispec,iorder)
        end do
        read(12,*) !! Blank line                            
                           
        !! Pre-divide gibbs coefficients (and leave in molar form). Also include log(P_ref/R0)
        !! and a log(T) term, so result of creating gibbs includes all terms required for
        !! backwards rate calculation
        coef_gibbs(ispec,:) = coef_cp(ispec,:)
        do iorder = 2,polyorder_cp+1
           coef_gibbs(ispec,iorder) = coef_gibbs(ispec,iorder)/dble(iorder*(iorder-1))
        end do                           
        coef_gibbs(ispec,1) = coef_cp(ispec,polyorder_cp+3) &
                            - coef_cp(ispec,1) &
                            + log(P_ref/Rgas_universal)
        coef_gibbs(ispec,polyorder_cp+3) = coef_cp(ispec,1) - one
                     
                           
        !! Convert cp coefficients from molar to mass based
        coef_cp(ispec,:) = coef_cp(ispec,:)*Rgas_universal*one_over_molar_mass(ispec)

        !! Pre-divide coefs by iorder for h.
        do iorder = 1,polyorder_cp + 1
           coef_h(ispec,iorder) = coef_cp(ispec,iorder)/dble(iorder)
        end do
        coef_h(ispec,polyorder_cp+2) = coef_cp(ispec,polyorder_cp+2)
        
        !! Pre-multiply coefs by iorder-1 for dcp/dT
        do iorder = 1,polyorder_cp + 1
           coef_dcpdT(ispec,iorder) = coef_cp(ispec,iorder)*dble(iorder-1)
        end do
                
     end do     

     !! Next section reads in reaction mechanism. 
     read(12,*) !! Read comment line    
     
     !! Number of steps
     read(12,*)
     read(12,*) nsteps
     read(12,*)
     
     !! Number of different third body efficiencies
     read(12,*)
     read(12,*) nthirdbodies
     read(12,*)
     
     
     !! Space for rate constants and coefficients etc
     allocate(Arrhenius_coefs(nsteps,3))
     
     !! Numbers of r and p, and lists
     allocate(num_reactants(nsteps),num_products(nsteps))
     allocate(reactant_list(nsteps,3),product_list(nsteps,3))  !! Limit of 3 reactants and 3 products in any given step
     allocate(stepspecies_list(nsteps,6)) !! List of all species in reaction (reactants then products)
     
     !! Stoichiometric coefficients
     allocate(nu_dash(nsteps,nspec),nu_ddash(nsteps,nspec),delta_nu(nsteps,nspec)) 
     nu_dash = zero;nu_ddash = zero;delta_nu = zero
     
     !! Flags, efficiencies and coefficients
     allocate(gibbs_rate_flag(nsteps),lindemann_form_flag(nsteps)) !! Flags for backwards and lindemann
     allocate(third_body_flag(nsteps))
     if(nthirdbodies.ne.0) allocate(third_body_efficiencies(nthirdbodies,nspec))
     allocate(lindemann_coefs(nsteps,4));lindemann_coefs = zero
    
     
     read(12,*) !! Comment line
     do istep = 1,nsteps
        read(12,*) dummy_int  !! Reaction number
        if(dummy_int.ne.istep) write(6,*) "Warning, error reading reaction mech. Expect seg fault."

        read(12,*) num_reactants(istep)            !! Number of reactants
        do ispec = 1,num_reactants(istep)  !! Loop over all reactants
           read(12,*) dummy_int,dummy_real
           reactant_list(istep,ispec) = dummy_int !! Identity of reactant
           stepspecies_list(istep,ispec) = dummy_int
           nu_dash(istep,dummy_int) = dummy_real  !! Stoichiometric coefficient
        end do

        read(12,*) num_products(istep)            !! Number of products
        do ispec = 1,num_products(istep)  !! Loop over all products
           read(12,*) dummy_int,dummy_real
           product_list(istep,ispec) = dummy_int !! Identity of product
           stepspecies_list(istep,num_reactants(istep) + ispec) = dummy_int
           nu_ddash(istep,dummy_int) = dummy_real  !! Stoichiometric coefficient
        end do
        
        !! Calculate deltas
        delta_nu(istep,:) = nu_ddash(istep,:) - nu_dash(istep,:)
        
        !! Coefficients for arrhenius rate constant
        read(12,*) Arrhenius_coefs(istep,1:3)

        !! Take logarithm of pre-exponential factor
        Arrhenius_coefs(istep,1) = log(Arrhenius_coefs(istep,1))
        Arrhenius_coefs(istep,3) = Arrhenius_coefs(istep,3)/(Rgas_universal)                     
        
        !! Gibbs based backwards rate?
        read(12,*) gibbs_rate_flag(istep)
        
        !! Lindemann form?
        read(12,*) lindemann_form_flag(istep)
        
        !! Third bodies?
        read(12,*) third_body_flag(istep)        
              
        read(12,*) !! Blank line      
     end do

     !! Lists of third body efficiencies
     if(nthirdbodies.ne.0) then
        
        do istep = 1,nthirdbodies
           read(12,*) !! Comment line
           do ispec = 1,nspec
              read(12,*) dummy_int,third_body_efficiencies(istep,ispec)     !! List of efficiencies
           end do
        end do
     end if
     read(12,*) !! Blank line
     
     !! Lists of Lindemann coefficients
     nlindemann = maxval(lindemann_form_flag(1:nsteps))
     if(nlindemann.ne.0) then
        read(12,*) !! Comment line
        do istep = 1,nlindemann
           read(12,*) dummy_int,lindemann_coefs(istep,1:4)
           !! Take logarithm of pre-exponential factor etc
           lindemann_coefs(istep,1) = log(lindemann_coefs(istep,1))
           lindemann_coefs(istep,3) = lindemann_coefs(istep,3)/(Rgas_universal)
           lindemann_coefs(istep,4) = log(lindemann_coefs(istep,4))
        end do         
     end if
     
     !! Build list of species which require gibbs evaluation
     allocate(gibbs_flag_species(nspec));gibbs_flag_species=0     
     do istep = 1,nsteps
        if(gibbs_rate_flag(istep).ne.0) then          !! For gibbs steps
           do jspec = 1, num_reactants(istep) + num_products(istep)  !! Loop over products and reactants
              ispec = stepspecies_list(istep,jspec)              
              gibbs_flag_species(ispec) = 1              !! Flag this species as requiring gibbs evaluation
           end do
        end if
     end do
     !! Count the number of species which require gibbs evaluation
     num_gibbs_species=0
     do ispec=1,nspec
        if(gibbs_flag_species(ispec).ne.0)then
           num_gibbs_species = num_gibbs_species  + 1           
           gibbs_flag_species(ispec) = num_gibbs_species  !! Modify flag to become a counter
        end if
     end do
     
     close(12)               

     
     return
  end subroutine load_chemistry_data    
!! ------------------------------------------------------------------------------------------------
  subroutine make_1d_1step_flame
     !! Make a single-step fixed stoichiometry flame
     integer(ikind) :: i,ispec,j
     real(rkind) :: flame_location,flame_thickness
     real(rkind) :: T_reactants,T_products
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z

     !! Position and scale     
     flame_location = zero
     flame_thickness = 5.0d-4/L_char !! Scale thickness because position vectors are scaled...

     !! Temperatures
     T_reactants = 3.0d2
     T_products = 2.3d3
     
     !! Pressure through flame
     P_flame = rho_char*Rgas_universal*T_reactants*one_over_molar_mass(1) 
     
     !! Inflow speed
     u_reactants = u_inflow
     
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Error function based progress variable
        c = half*(one + erf((x-flame_location)/flame_thickness))
        
        !! Temperature profile
        T(i) = T_reactants + (T_products - T_reactants)*c
        
        !! Composition
        Yspec(i,1) = one - c
        Yspec(i,2) = c
        
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal
        
        !! Density
        lnro(i) = log(P_flame) - log(Rmix_local) - log(T(i))
        
        !! Velocity
        u(i) = u_reactants*rho_char/exp(lnro(i))
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
!              T(i) = T_reactants + half*half*(T_products-T_reactants)
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char                
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              !! Do nothing
           end if
           T_bound(j) = T(i)  !! set T_bound to T
        end do
     end if
  
     return
  end subroutine make_1d_1step_flame
!! ------------------------------------------------------------------------------------------------
  subroutine make_1d_21step_flame
     !! Make a 21 step (9 species) H2-AIR flame. Initial profiles are erf(x').
     integer(ikind) :: i,ispec,j
     real(rkind) :: flame_location,flame_thickness
     real(rkind) :: T_reactants,T_products
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z
     real(rkind) :: Yin_H2,Yin_O2,Yin_N2,Yout_H2O

     !! Position and scale     
     flame_location = zero
     flame_thickness = 5.0d-4/L_char !! Scale thickness because position vectors are scaled...

     !! Inlet composition
     Yin_H2 = 0.0283126
     Yin_O2 = 0.226501
     Yin_N2 = 0.745187
     
     !! Outlet composition
     Yout_H2O = one - Yin_N2     

     !! Temperatures
     T_reactants = 3.0d2
     T_products = 2.366d3
     
     !! Pressure through flame
     P_flame = rho_char*Rgas_universal*T_reactants* &
             (Yin_H2*one_over_molar_mass(1) + Yin_O2*one_over_molar_mass(2) + Yin_N2*one_over_molar_mass(9))
     
     !! Inflow speed
     u_reactants = u_inflow
     
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Error function based progress variable
        c = half*(one + erf((x-flame_location)/flame_thickness))
        
        !! Temperature profile
        T(i) = T_reactants + (T_products - T_reactants)*c
        
        !! Composition
        Yspec(i,1) = (one - c)*Yin_H2
        Yspec(i,2) = (one - c)*Yin_O2
        Yspec(i,3) = c*Yout_H2O
        Yspec(i,4:8) = zero
        Yspec(i,9) = Yin_N2
        
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal
        
        !! Density
        lnro(i) = log(P_flame) - log(Rmix_local) - log(T(i))
        
        !! Velocity
        u(i) = u_reactants*rho_char/exp(lnro(i))
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              T_bound(j) = T(i)
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
              T_bound(j) = T(i) !! Inflow temperature is T_cold
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              T_bound(j) = T(i)  !! Outflow temperature is T_hot
           end if
        end do
     end if
 
  
  
     return
  end subroutine make_1d_21step_flame
!! ------------------------------------------------------------------------------------------------  
  subroutine make_1d_25step_flame
     !! Make a 25 step (16 species) CH4-AIR flame. Initial profiles are erf(x').
     integer(ikind) :: i,ispec,j
     real(rkind) :: flame_location,flame_thickness
     real(rkind) :: T_reactants,T_products
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z
     real(rkind) :: Yin_CH4,Yin_O2,Yin_N2,Yout_H2O,Yout_CO2

     !! Position and scale     
     flame_location = zero
     flame_thickness = 5.0d-4/L_char !! Scale thickness because position vectors are scaled...

     !! Inlet composition
     Yin_CH4 = 0.055046
     Yin_O2 = 0.22018
     Yin_N2 = one - Yin_CH4 - Yin_O2
     
     !! Outlet composition
     Yout_CO2 = (one - Yin_N2)*44.0d0/80.0d0
     Yout_H2O = Yout_CO2*36.0d0/44.0d0

     !! Temperatures
     T_reactants = 3.0d2
     T_products = 2.3d3
     
     !! Pressure through flame
     P_flame = rho_char*Rgas_universal*T_reactants* &
             (Yin_CH4*one_over_molar_mass(1) + Yin_O2*one_over_molar_mass(2) + Yin_N2*one_over_molar_mass(16))
     
     !! Inflow speed
     u_reactants = u_inflow
     
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Error function based progress variable
        c = half*(one + erf((x-flame_location)/flame_thickness))
        
        !! Temperature profile
        T(i) = T_reactants + (T_products - T_reactants)*c
        
        !! Composition
        Yspec(i,1) = (one - c)*Yin_CH4
        Yspec(i,2) = (one - c)*Yin_O2
        Yspec(i,3) = c*Yout_CO2
        Yspec(i,4) = c*Yout_H2O
        Yspec(i,5:15) = zero
        Yspec(i,16) = Yin_N2
        
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal
        
        !! Density
        lnro(i) = log(P_flame) - log(Rmix_local) - log(T(i))
        
        !! Velocity
        u(i) = u_reactants*rho_char/exp(lnro(i))
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              T_bound(j) = T(i)
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
              T_bound(j) = T(i) !! Inflow temperature is T_cold
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              T_bound(j) = T(i)  !! Outflow temperature is T_hot
           end if
        end do
     end if
 
  
  
     return
  end subroutine make_1d_25step_flame  
!! ------------------------------------------------------------------------------------------------
  subroutine load_flame_file
     use thermodynamics
     integer(ikind) :: i,ispec,j,nflamein
     real(rkind) :: flame_location,flame_thickness,cell_pos
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z
     real(rkind),dimension(:),allocatable :: flamein_ro,flamein_u,flamein_v,flamein_w,flamein_roE
     real(rkind),dimension(:,:),allocatable :: flamein_Y
     real(rkind) :: dx_flamein

     !! Open a file containing flame profile
     open(unit=19,file='init_flame.in')
     
     !! Specify the flame pressure
     P_flame = 1.0d5

     !! Read file, then close it.     
     !! Format of the file is ro,ro*u,ro*v,ro*w,roE,ro*Y (all ispec) over 1001 evenly spaced steps covering 
     !! a 0.01m long domain.
     nflamein = 1002
     dx_flamein = 0.01/(nflamein-2)
     allocate(flamein_ro(nflamein),flamein_u(nflamein),flamein_v(nflamein),flamein_w(nflamein))
     allocate(flamein_Y(nflamein,nspec))
     allocate(flamein_roE(nflamein))
     do i=1,nflamein-1
        !! Load conservative data
        read(19,*) flamein_ro(i),flamein_u(i),flamein_v(i),flamein_w(i), &
                   flamein_roE(i),flamein_Y(i,1:nspec)

        !! Divide by rho as required
        flamein_u(i) = flamein_u(i)/flamein_ro(i)
        flamein_v(i) = flamein_v(i)/flamein_ro(i)
        flamein_w(i) = flamein_w(i)/flamein_ro(i)        
        flamein_Y(i,:) = flamein_Y(i,:)/flamein_ro(i)                
     end do
     flamein_ro(nflamein) = flamein_ro(nflamein-1)   !! Data for extra point at end
     flamein_u(nflamein) = flamein_u(nflamein-1)
     flamein_v(nflamein) = flamein_v(nflamein-1)
     flamein_w(nflamein) = flamein_w(nflamein-1)
     flamein_roE(nflamein) = flamein_roE(nflamein-1)
     flamein_Y(nflamein,:) = flamein_Y(nflamein-1,:)   
     close(19)
     
       
     !! Loop through all particles. Find the "cell" the particle resides in. Copy data.     
     !$omp parallel do private(j,x,y,z,c,Rmix_local,ispec,cell_pos)
     do i=1,npfb
        x = (rp(i,1)+half)*L_char  !! Scale x - this requires L_char to match the length
        !! also needs shifting depending on set up.
        
        !! Nearest index in flame-in data (to left of x)
        j = floor(x/dx_flamein) + 1      
        
        !! proportion along cell
        cell_pos= x/dx_flamein - dble(j-1)
        
        !! Copy data
        lnro(i) = log(flamein_ro(j)*(one - cell_pos) + flamein_ro(j+1)*cell_pos)
        u(i) = flamein_u(j)*(one - cell_pos) + flamein_u(j+1)*cell_pos
        v(i) = flamein_v(j)*(one - cell_pos) + flamein_v(j+1)*cell_pos
        w(i) = flamein_w(j)*(one - cell_pos) + flamein_w(j+1)*cell_pos
        roE(i) = flamein_roE(j)*(one - cell_pos) + flamein_roE(j+1)*cell_pos   
        Yspec(i,:) = flamein_Y(j,:)*(one - cell_pos) + flamein_Y(j+1,:)*cell_pos

        !! Temporary assign p and T
        T(i) = T_ref
        p(i) = P_flame
        
     end do
     !$omp end parallel do
     
     !! Free up space
     deallocate(flamein_ro,flamein_u,flamein_v,flamein_w,flamein_roE,flamein_Y)
     
     !! Temporarily copy some energy data to halos and mirrors (it will be later overwritten, but
     !! just prevents the NR solver from crashing at set-up)
     !$omp parallel do
     do i=npfb+1,np
        roE(i) = roE(1)
        lnro(i) = lnro(1)
        u(i) = u(1);v(i) = v(1);w(i) = w(1)
        Yspec(i,:) = Yspec(1,:)
     end do
     !$omp end parallel do
     
     !! Re-evaluate temperature from energy.  
     call evaluate_temperature_and_pressure
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
              T_bound(j) = T(i)
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
              T_bound(j) = T(i) !! Inflow temperature is T_cold
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              T_bound(j) = T(i)  !! Outflow temperature is T_hot
           end if
        end do
     end if   
  
     return
  end subroutine load_flame_file
!! ------------------------------------------------------------------------------------------------
  subroutine hardcode_initial_conditions
     !! Temporary routine to generate initial conditions from some hard-coded functions.
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro
     
     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(x,y,z,tmp,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        u(i) = u_char!-cos(two*pi*x)*sin(two*pi*y)!*cos(two*pi*z/Lz)!*oosqrt2
        v(i) = zero!sin(two*pi*x)*cos(two*pi*y)!*cos(two*pi*z/Lz)    !!c c
        w(i) = zero!u(i);u(i)=zero
!        tmp = -half*half*(cos(4.0d0*pi*x) + cos(4.0d0*pi*y))/csq  !! Modify for not(isoT)
        lnro(i) = log(rho_char)!log(rho_char + tmp)
!if(x.le.zero) lnro(i) = log(1.2*rho_char)             

        tmp = half*(one+erf(5.0d0*x))
        T(i) = T_ref    
        tmp = rho_char*T_ref/T(i)
        lnro(i) = log(tmp)                                            
              
#ifdef ms    
!        do ispec=1,nspec      
           tmp = one - half*(one + erf(5.0d0*x))
           Yspec(i,1) = tmp
!           Yspec(i,2) = one - tmp

!        end do
#endif         
                    
     end do
     !$OMP END PARALLEL DO
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero

              tmp = T_ref*(one + 0.01*sin(two*pi*rp(i,3)/Lz))
              T_bound(j) = tmp !! Might want changing                
              T(i) = T_bound(j)
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
              u(i)=u_char
           end if
        end do
     end if   
  
     return
  end subroutine hardcode_initial_conditions  
!! ------------------------------------------------------------------------------------------------
  subroutine load_restart_file(n_restart)
     !! Load initial conditions from a dump file
     integer(ikind),intent(in) :: n_restart
     integer(ikind) :: k,i,j
     real(rkind) :: tmp,tmpro
     character(70) :: fname  

#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     
     if( n_restart .lt. 10 ) then 
        write(fname,'(A17,I5,A1,I1)') './data_out/layer_',k,'_',n_restart
     else if( n_restart .lt. 100 ) then 
        write(fname,'(A17,I5,A1,I2)') './data_out/layer_',k,'_',n_restart        
     else if( n_restart .lt. 1000 ) then
        write(fname,'(A17,I5,A1,I3)') './data_out/layer_',k,'_',n_restart        
     else
        write(fname,'(A17,I5,A1,I4)') './data_out/layer_',k,'_',n_restart        
     end if 
     !! Open the file
     open(14,file=fname)
     read(14,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart"
     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(14,*) tmp,tmp,tmp,tmp,h(i),k,tmpro,u(i),v(i),w(i),tmp,T(i),Yspec(i,1:nspec)
#else
        read(14,*) tmp,tmp,tmp,h(i),k,tmpro,u(i),v(i),tmp,T(i),Yspec(i,1:nspec)
#endif        
        if(k.ne.node_type(i)) then
           write(6,*) "ERROR: Problem in restart file. STOPPING."
#ifdef mp
           call MPI_Abort(MPI_COMM_WORLD, k, ierror)
#else
           stop
#endif
        end if
        lnro(i) = log(tmpro)
     end do
     
     !! Re-specify the boundary temperatures
     if(nb.ne.0) then
        do j=1,nb
           i=internal_list(j)
           T_bound(j) = T(i)
        end do
     end if
      
     close(14)
  
     return
  end subroutine load_restart_file
!! ------------------------------------------------------------------------------------------------
end module setup_flow
