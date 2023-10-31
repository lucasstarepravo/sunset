module setup_flow
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains routines generate initial conditions, either from restart files or
  !! analytic forms.
  use kind_parameters
  use common_parameter
  use common_vars
#ifdef mp
  use mpi_transfers
#endif    
  use turbulence
  implicit none
  
  real(rkind),dimension(:),allocatable :: Yspec_reactants,Yspec_products
  
contains
!! ------------------------------------------------------------------------------------------------
  subroutine initial_solution
     !! Allocates arrays for primary and secondary properties, and calls routines to populate arrays
     !! for initial conditions.
     use mirror_boundaries
     use derivatives
     use thermodynamics
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro     
     
     !! Allocate arrays for properties - primary
     allocate(rou(np),rov(np),row(np),ro(np),roE(np),divvel(np))
     allocate(Yspec(np,nspec))
     rou=zero;rov=zero;row=zero;ro=one;roE=one;Yspec=one;divvel=zero

     !! Secondary properties
     allocate(T(np));T=T_ref
     allocate(p(np));p=zero
     allocate(u(np),v(np),w(np));u=zero;v=zero;w=zero
     allocate(hrr(npfb));hrr = zero !! Array for heat release rate  
     
     !! Transport properties
     allocate(visc(np));visc = visc_ref
#ifndef isoT
     allocate(lambda_th(np))
#endif
     allocate(roMdiff(np,nspec))
     allocate(cp(np),Rgas_mix(np))
     
     !! Allocate the boundary temperatures
     if(nb.ne.0) then
        allocate(T_bound(nb));T_bound = T_ref  
        allocate(u_inflow_local(nb));u_inflow_local = u_char   
        allocate(dudt_inflow_local(nb));dudt_inflow_local=zero            
     end if
     
     !! =======================================================================
     !! Choose initial conditions
#ifndef restart     
     !! Evaluate mass fractions of reactants and products assuming complete combustion
     call initialise_composition        

     !! Make the base flow
     call make_baseflow   
     
     !! HARD-CODED-CHOICE =====================================
     !! Un-comment ONE of the following routines ==============
     !! The routines themselves can be modified a bit to change
     !! things like the shape of a hotspot
        
        
     !! Make a 1D flame: pass X-position,flame_thickness and T_hot
     call make_1d_flame(-0.0d0,2.0d-4,2.366d3)  !-0.275d0
     
     !! Make a gaussian hotspot: pass X,Y-positions, hotspot size and T_hot
!     call make_gaussian_hotspot(-0.23d0,zero,2.0d-4,2.5d3)   !-0.23d0 !0.045    

     !! Load an existing 1D flame file
!     call load_flame_file

     !! Make the initial conditions for a Rayleigh-Taylor instability
!     call make_RT_initial_conditions(-0.0d0,2.0d-4,two*T_ref)
     
     !! A messy routine to play with for other initial conditions
!     call hardcode_initial_conditions     

     !! END HARD-CODED-CHOICE =================================

     deallocate(Yspec_reactants,Yspec_products)
#else    
     !! RESTART OPTION. 
     call load_restart_file

     !! Un-comment this to (re-)ignite a restarting simulation
!     call make_gaussian_hotspot(-2.0d0,zero,2.0d-4,2.5d3)   !-0.23d0 !0.045    

#endif

     !! Add some turbulence to the velocity field
!     call make_turbulent_velocity_field(1.0d-3,5.0d0*u_char)
     !! =======================================================================
     
     !! Convert from velocity to momentum and Y to roY
     !$omp parallel do private(ispec)
     do i=1,np
        rou(i) = ro(i)*u(i)
        rov(i) = ro(i)*v(i)
        row(i) = ro(i)*w(i)                
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec(i,ispec)*ro(i)
        end do
     end do
     !$omp end parallel do   
              
     !! Set energy from ro,u,Y,T
     call initialise_energy     
     
     !! Mirrors and halos                        
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchanges_all
#endif          

     !! Obtain velocity from momentum for mirrors
     !$omp parallel do
     do i=npfb+1,np
        u(i)=rou(i)/ro(i)
        v(i)=rov(i)/ro(i)
        w(i)=row(i)/ro(i)                
     end do
     !$omp end parallel do
     
     !! Set the initial velocity divergence
#ifdef dim3        
     call calc_divergence(u,v,w,divvel(1:npfb))
#else
     call calc_divergence(u,v,divvel(1:npfb))     
#endif              
     
     !! Mirrors and halos for divvel                   
     call reapply_mirror_bcs
#ifdef mp
     call halo_exchange_divvel
#endif         
         
     !! Initialise the variable which holds inflow velocity locally
     if(nb.ne.0)then
        !$omp parallel do
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then !! Inflow node
              u_inflow_local(j) = u(i)
           endif
        end do
        !$omp end parallel do        
     end if
     
     !! Initialise the variable holding the inflow species     
     if(nb.ne.0) then
        allocate(Yspec_inflow(nspec))
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.1) then !! Inflow node
              do ispec=1,nspec
                 Yspec_inflow(ispec) = Yspec(i,ispec)/ro(i)
              end do
           end if
        end do
     endif
           
     !! SOME ADDITIONAL INITIALISATION STUFF ------------------------------------------------------
     !! Profiling - re-zero time accumulators
     segment_time_local = zero
     cputimecheck = zero         
     
#ifndef restart
     !! Pre-set the time-stepping (necessary for PID controlled stepping)
     dt = 1.0d-10            
     !! Initialise PID controller variables
     emax_np1=pid_tol;emax_n=pid_tol;emax_nm1=pid_tol
#endif     
     
     !! Initialise time-stepping error normalisation based on expected magnitudes
     ero_norm = one/one !! Unity density
     erou_norm = one/(one*u_char)  !! Characteristic velocity
     eroE_norm = one/p_ref         !! Energy of the order of pressure
     eroY_norm = one/(one*one)     !! Unity mass fraction
     

     !! Initialise PID controller variables for <|u|>
#ifndef restart
     eflow_nm1 = one
     sum_eflow = zero   
     driving_force(:) = zero !! Will only be changed if using PID controller
#endif     
                   
     return
  end subroutine initial_solution   
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
!! N.B. In the routines below here, we are loading or generating initial conditions on the
!! primitive variables (ro,u,v,w,T,Y), and hence Yspec holds Y. Everywhere else in the code, 
!! Yspec holds roY.
!! ------------------------------------------------------------------------------------------------
  subroutine initialise_composition
     !! This subroutine evaluates the reactant and product mass fractions for the given
     !! stoichiometry.
     real(rkind) :: Yin_H2,Yin_O2,Yin_N2,Yout_H2O
     real(rkind) :: o2n2_ratio,h2o2_stoichiometric,h2o2_ratio
     real(rkind) :: Yin_CH4,Yout_CO2
     real(rkind) :: ch4o2_ratio,ch4o2_stoichiometric,co2h2o_ratio     
  
     !! Allocate two arrays
     allocate(Yspec_reactants(nspec),Yspec_products(nspec))
  
     !! Single-phase
     if(nspec.eq.1) then
        Yspec_reactants(1)=one
        Yspec_products(1)=one
     end if
  
     !! Single-step chemistry
     if(nspec.eq.2) then 
        Yspec_reactants(1) = one
        Yspec_reactants(2) = zero
        Yspec_products(1) = zero
        Yspec_products(2) = one
     end if
    
     !! 9-species H2 chemistry
     if(nspec.eq.9) then
        o2n2_ratio = one*molar_mass(2)/(3.76d0*molar_mass(9))
        h2o2_stoichiometric = two*molar_mass(1)/(one*molar_mass(2))
        h2o2_ratio = h2o2_stoichiometric*phi_in

        Yin_O2 = one/(one + h2o2_ratio + one/o2n2_ratio)
        Yin_H2 = Yin_O2*h2o2_ratio
        Yin_N2 = one - Yin_H2 - Yin_O2            
        Yout_H2O = one - Yin_N2      

        Yspec_reactants(1) = Yin_H2
        Yspec_reactants(2) = Yin_O2
        Yspec_reactants(3:8) = zero
        Yspec_reactants(9) = Yin_N2
        
        Yspec_products(1:2) = zero
        Yspec_products(3) = Yout_H2O
        Yspec_products(4:8) = zero
        Yspec_products(9) = Yin_N2      
     end if
     
     !! 15 species CH4 chemistry
     if(nspec.eq.16) then
        o2n2_ratio = one*molar_mass(2)/(3.76d0*molar_mass(16))
        ch4o2_stoichiometric = one*molar_mass(1)/(two*molar_mass(2))
        ch4o2_ratio = ch4o2_stoichiometric*phi_in
          
        Yin_O2 = one/(one + ch4o2_ratio + one/o2n2_ratio)
        Yin_CH4 = Yin_O2*ch4o2_ratio
        Yin_N2 = one - Yin_CH4 - Yin_O2
   
        co2h2o_ratio = one*molar_mass(3)/(two*molar_mass(4))  !! Assume complete combustion??
        Yout_CO2 = (one - Yin_N2)/(one + one/co2h2o_ratio)
        Yout_H2O = Yout_CO2/co2h2o_ratio        
     
        Yspec_reactants(1) = Yin_CH4
        Yspec_reactants(2) = Yin_O2
        Yspec_reactants(3:15) = zero
        Yspec_reactants(16) = Yin_N2
        
        Yspec_products(1:2) = zero
        Yspec_products(3) = Yout_CO2
        Yspec_products(4) = Yout_H2O
        Yspec_products(5:15) = zero
        Yspec_products(16) = Yin_N2   
     end if
  
     return
  end subroutine initialise_composition
!! ------------------------------------------------------------------------------------------------ 
  subroutine make_1d_flame(flame_location,flame_thickness,T_products)
     integer(ikind) :: i,ispec,j
     real(rkind),intent(in) :: flame_location,flame_thickness,T_products
     real(rkind) :: fl_thck
     real(rkind) :: c,Rmix_local,x,y,z,ro_inflow

     !! Scale thickness because position vectors are scaled...
     fl_thck = flame_thickness/L_char 

     !! Inflow mixture gas constant
     Rmix_local = zero
     do ispec=1,nspec
        Rmix_local = Rmix_local + Yspec_reactants(ispec)*one_over_molar_mass(ispec)
     end do
     Rmix_local = Rmix_local*Rgas_universal     

     !! Inflow density
     ro_inflow = p_ref/(Rmix_local*T_ref)
    
     !! Loop over all nodes and impose values
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
                
        !! Error function based progress variable
        z=0.01*sin(two*pi*y/(ymax-ymin))
        c = half*(one + erf((x-flame_location)/fl_thck))
        
        !! Temperature profile
        T(i) = T(i) + (T_products - T(i))*c
        
        !! Composition
        do ispec = 1,nspec
           Yspec(i,ispec) = (one-c)*Yspec_reactants(ispec) + c*Yspec_products(ispec)
        end do
        
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal
        
        !! Density
        ro(i) = p_ref/(Rmix_local*T(i))
        
        !! Adjust the velocity
        u(i) = u(i)*ro_inflow/ro(i)
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Special apply values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
           end if                 
           T_bound(j) = T(i)
        end do
     end if
 
  
  
     return
  end subroutine make_1d_flame
!! ------------------------------------------------------------------------------------------------  
  subroutine make_RT_initial_conditions(flame_location,flame_thickness,T_products)
     integer(ikind) :: i,ispec,j
     real(rkind),intent(in) :: flame_location,flame_thickness,T_products
     real(rkind) :: fl_thck
     real(rkind) :: c,Rmix_local,x,y,z,ro_ref_hot,ro_ref_cold,ro_ref_local
     real(rkind) :: T_y0,dcdy,p_local

     !! Scale thickness because position vectors are scaled...
     fl_thck = flame_thickness/L_char 

     !! Inflow mixture gas constant
     Rmix_local = zero
     do ispec=1,nspec
        Rmix_local = Rmix_local + Yspec_reactants(ispec)*one_over_molar_mass(ispec)
     end do
     Rmix_local = Rmix_local*Rgas_universal     

     !! Inflow density
     ro_ref_cold = p_ref/(Rmix_local*T_ref)
     ro_ref_hot = p_ref/(Rmix_local*T_products)
         
     !! Loop over all nodes and impose values
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec,dcdy,ro_ref_local)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
                
        !! Error function based progress variable
        c = half*(one + erf((flame_location-x)/fl_thck))        

        !! Temperature profile
         
        !! Composition
        do ispec = 1,nspec
!           Yspec(i,ispec) = (one-c)*Yspec_reactants(ispec) + c*Yspec_products(ispec)
        end do
        
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal
        
        
        if(x.lt.zero) then !! Cold region
           T(i) = T_ref
           ro_ref_local = ro_ref_cold
        else               !! Hot region
           T(i) = T_products
           ro_ref_local = ro_ref_hot
        end if
        
        !! Density
!        ro(i) = ro_ref_local*exp(grav(2)*x/(Rmix_local*T(i)))
        ro(i) = p_ref/(Rmix_local*T(i))        
        
        !! Local pressure
        p(i) = ro(i)*Rmix_local*T(i)        
        
        !! Adjust the velocity
        u(i) = u(i)!*ro_ref_cold/ro(i)        
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Special apply values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
           end if                 
           T_bound(j) = T(i)
        end do
     end if
 
  
  
     return
  end subroutine make_RT_initial_conditions   
!! ------------------------------------------------------------------------------------------------ 
  subroutine make_gaussian_hotspot(f_loc_x,f_loc_y,flame_thickness,T_hot)
     !! Add a Gaussian hotspot to the flow. The composition and velocity are unchanged.
     !! Temperature follows a prescribed Gaussian, whilst density is modified to ensure pressure
     !! is unchanged.
     integer(ikind) :: i,ispec,j
     real(rkind),intent(in) :: f_loc_x,f_loc_y,flame_thickness,T_hot
     real(rkind) :: fl_thck,pR_local
     real(rkind) :: c,u_reactants,Rmix_local,x,y,z

     !! Scale thickness because position vectors are scaled...
     fl_thck = flame_thickness/L_char 

   
     !! Loop over all nodes and impose values
     !$omp parallel do private(x,y,z,c,pR_local)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Get local pressure
        pR_local = ro(i)*T(i) 
               
        !! Gaussian progress variable
!        c = exp(-((x-f_loc_x)/fl_thck)**two - ((y-f_loc_y)/fl_thck)**two) !! 2D
        c = exp(-((x-f_loc_x)/fl_thck)**two)                               !! 1D 
                
        !! Adjust temperature profile
        T(i) = T(i)*(one-c) + T_hot*c
        
        !! Composition - unchanged
       
        !! Density
        ro(i) = pR_local/(T(i))
        
        !! Velocity
        !! ...is unchanged.
                        
     end do
     !$omp end parallel do
     
     !! Special impose values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
           end if                 
           T_bound(j) = T(i)
        end do
     end if
    
     return
  end subroutine make_gaussian_hotspot
!! ------------------------------------------------------------------------------------------------
  subroutine make_baseflow
     !! Initialise the base flow: uniform ro, T, p, Y=Y_reactants, and spatially varying u,v,w
     integer(ikind) :: i,ispec,j
     real(rkind) :: Rmix_local,x,y,z,ro_inflow


#ifndef isoT
     !! Inflow mixture gas constant
     Rmix_local = zero
     do ispec=1,nspec
        Rmix_local = Rmix_local + Yspec_reactants(ispec)*one_over_molar_mass(ispec)
     end do
     Rmix_local = Rmix_local*Rgas_universal     

     !! Inflow density based on reference pressure and temperature
     ro_inflow = p_ref/(Rmix_local*T_ref)
#else
     ro_inflow = rho_char
#endif     
    
     !! Loop over all nodes and impose values
     !$omp parallel do private(x,y,z,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
              
        !! Temperature profile
        T(i) = T_ref
        
        !! Composition - reactants everywhere
        do ispec = 1,nspec
           Yspec(i,ispec) = Yspec_reactants(ispec)
        end do
              
        !! Density
        ro(i) = ro_inflow
        
        !! Velocity
        if(flag_base_flow_profile.eq.0) then !! Zero flow
           u(i) = zero
        else if(flag_base_flow_profile.eq.1) then !! Uniform flow
           u(i) = u_inflow_start
        else if(flag_base_flow_profile.eq.2) then !! Parabolic flow
           y=y/(ymax-ymin)
           u(i) = u_inflow_start*six*(half-y)*(half+y)         
        end if
        v(i) = zero
        w(i) = zero
                        
     end do
     !$omp end parallel do
     
     !! Special impose values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero  !! Will impose an initial shock!!
           end if                 
           T_bound(j) = T(i)
        end do
     end if
    
     return
  end subroutine make_baseflow  
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
     P_flame = p_ref

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
        ro(i) = flamein_ro(j)*(one - cell_pos) + flamein_ro(j+1)*cell_pos
        u(i) = flamein_u(j)*(one - cell_pos) + flamein_u(j+1)*cell_pos
        v(i) = flamein_v(j)*(one - cell_pos) + flamein_v(j+1)*cell_pos
        w(i) = flamein_w(j)*(one - cell_pos) + flamein_w(j+1)*cell_pos
        roE(i) = flamein_roE(j)*(one - cell_pos) + flamein_roE(j+1)*cell_pos   
        Yspec(i,:) = flamein_Y(j,:)*(one - cell_pos) + flamein_Y(j+1,:)*cell_pos
        
        !! Convert to conservative variables
        rou(i) = u(i)*ro(i)
        rov(i) = v(i)*ro(i)
        row(i) = w(i)*ro(i)
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec(i,ispec)*ro(i)
        end do

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
        ro(i) = ro(1)
        u(i) = u(1);v(i) = v(1);w(i) = w(1)
        rou(i) = rou(1);rov(i)=rov(1);row(i) = row(1)
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
  
     !! Return Yspec to primitive form
     !$omp parallel do private(ispec)
     do i=1,npfb
        do ispec=1,nspec
           Yspec(i,ispec) = Yspec(i,ispec)/ro(i)
        end do
     end do
     !$omp end parallel do
     
  
     return
  end subroutine load_flame_file
!! ------------------------------------------------------------------------------------------------
  subroutine hardcode_initial_conditions
     !! Temporary routine to generate initial conditions from some hard-coded functions.
     integer(ikind) :: i,j,k,ispec
     real(rkind) :: x,y,z,tmp,tmpro,Rmix_local
     
     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(x,y,z,tmp,ispec,Rmix_local)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)

        !! TG 3D Re1600 as in Cant 2022, Sandam 2017 etc (ish)
        u(i) = u_char!-cos(x)*sin(y)*cos(z)!*oosqrt2
        v(i) = zero!sin(x)*cos(y)*cos(z)    !!c c
        w(i) = zero!u(i);u(i)=zero
                
        
!        tmp = Rgas_universal*one_over_molar_mass(1)*T_ref  !! RT0        
!        tmp = rho_char*U_char*U_char/tmp   !! roUU/RT0
!        tmp = -tmp*(one/16.0d0)*(cos(two*x)+cos(two*y))*(two+cos(two*z))        
!        T(i) = T_ref

        ro(i) = rho_char! + tmp

!        tmp = one - half*(one + erf(5.0d0*x))
        Yspec(i,:) = one

!        u(i) = zero
!        v(i) = zero
 !       T(i) = T_ref
 !       tmp = exp(-y*y/0.01d0)
 !       ro(i) = one + 0.1d0*tmp
        
#ifndef isoT
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal       
        
        !! Density
        p(i) = ro(i)*Rmix_local*T(i)
#else
        p(i) = ro(i)*csq
#endif        
        
                   
     end do
     !$OMP END PARALLEL DO
     
    
     
     !! Values on boundaries
     if(nb.ne.0)then
        do j=1,nb
           i=boundary_list(j)
           if(node_type(i).eq.0) then !! wall initial conditions
              u(i)=zero;v(i)=zero;w(i)=zero

              tmp = T_ref*(one + 0.01*sin(two*pi*rp(i,3)/L_domain_z))
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
  subroutine load_restart_file
     !! Load initial conditions from a dump file
     integer(ikind) :: k,i,j
     real(rkind) :: tmp,tmpro
     character(70) :: fname,fname2  

#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     

     !! Construct the file name:
     write(fname,'(A17,I5)') './restart/fields_',k
     write(fname2,'(A16,I5)') './restart/nodes_',k

     !! Load the "smoothing length" from nodes file
     open(15,file=fname2)
     read(15,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart. NODES FILE."
     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(15,*) tmp,tmp,tmp,tmp,h(i),k
#else
        read(15,*) tmp,tmp,tmp,h(i),k
#endif        
        if(k.ne.node_type(i)) then
           write(6,*) "ERROR: Problem in restart file. STOPPING."
#ifdef mp
           call MPI_Abort(MPI_COMM_WORLD, k, ierror)
#else
           stop
#endif
        end if
     end do
     close(15)

     !! Open the field files
     open(14,file=fname)
     read(14,*) !! Skip line
     read(14,*) k
     read(14,*) emax_np1,emax_n,emax_nm1,dt
     read(14,*) eflow_nm1,sum_eflow,driving_force
     read(14,*) !! Skip line
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart. FIELDS FILE."

     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(14,*) tmpro,u(i),v(i),w(i),tmp,T(i),p(i),hrr(i),Yspec(i,1:nspec)
#else
        read(14,*) tmpro,u(i),v(i),tmp,T(i),p(i),hrr(i),Yspec(i,1:nspec)
#endif        
        ro(i) = tmpro
        p(i) = p(i) + p_ref !! The output files hold p-p_ref
     end do      
     close(14)
     
     !! Re-specify the boundary temperatures
     if(nb.ne.0) then
        do j=1,nb
           i=internal_list(j)
           T_bound(j) = T_ref!T(i)
        end do
     end if

  
     return
  end subroutine load_restart_file
!! ------------------------------------------------------------------------------------------------
end module setup_flow
