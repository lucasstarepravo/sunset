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
  implicit none
  
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
     rou=zero;rov=zero;row=zero;ro=zero;roE=one;Yspec=one;divvel=zero

     !! Secondary properties
     allocate(T(np));T=T_ref
     allocate(p(np));p=zero
     allocate(u(np),v(np),w(np));u=zero;v=zero;w=zero
     allocate(alpha_out(np));alpha_out = zero   
     allocate(hrr(npfb));hrr = zero !! Array for heat release rate  
     
     !! Transport properties
     allocate(visc(np));visc = visc_ref
#ifndef isoT
     allocate(lambda_th(np))
#endif
#ifdef ms     
     allocate(roMdiff(np,nspec))
#endif     
     allocate(cp(np),Rgas_mix(np))
     
     !! Allocate the boundary temperatures
     if(nb.ne.0) then
        allocate(T_bound(nb));T_bound = T_ref  
        allocate(u_inflow_local(nb));u_inflow_local = u_char               
     end if
     
     !! =======================================================================
     !! Choose initial conditions
#ifndef restart     
#ifdef react
     if(.true.)then
        if(nsteps.eq.1) then
           call make_1d_1step_flame
        else if(nsteps.eq.21) then
           call make_1d_21step_flame
        else if(nsteps.eq.35) then
           call make_1d_25step_flame
        end if
      else
         call load_flame_file
      end if
#else
     call hardcode_initial_conditions     
#endif

#else    
     !! RESTART OPTION. 
     call load_restart_file
!     call make_ignition_hotspot         
#endif
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
     dt = 1.0d-10
     
     !! Initialise the variable which holds inflow velocity - this will need modifying in due course
     !! as we introduce corners and non-uniform inflows
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
              
     return
  end subroutine initial_solution   
!! ------------------------------------------------------------------------------------------------
!! ------------------------------------------------------------------------------------------------
!! N.B. In the routines below here, we are loading or generating initial conditions on the
!! primitive variables (ro,u,v,w,T,Y), and hence Yspec holds Y. Everywhere else in the code, 
!! Yspec holds roY.
!! ------------------------------------------------------------------------------------------------
  subroutine make_1d_1step_flame
     !! Make a single-step fixed stoichiometry flame
     integer(ikind) :: i,ispec,j
     real(rkind) :: flame_location,flame_thickness
     real(rkind) :: T_reactants,T_products
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z,ro_inflow

     !! Position and scale     
     flame_location = -0.0d0
     flame_thickness = 2.0d-4/L_char !! Scale thickness because position vectors are scaled...

     !! Temperatures
     T_reactants = T_ref
     T_products = 2.3d3
         
     !! Pressure through flame     
     P_flame = p_ref
     
     !! Inflow density
     ro_inflow = p_ref/(Rgas_universal*T_reactants*one_over_molar_mass(1))    
     
     !! Inflow speed
     u_reactants = u_char
     
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Error function based progress variable
        c = half*(one + erf((x-flame_location)/flame_thickness))
!        c = exp(-((x-flame_location)/flame_thickness)**two - (y/flame_thickness)**two)
                
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
        ro(i) = P_flame/(Rmix_local*T(i))   
               
        !! Velocity
        u(i) = u_reactants*ro_inflow/ro(i)
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
              T(i) = T_reactants !+ half*half*(T_products-T_reactants)              
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
!              u(i)=u_char                
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
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z,ro_inflow
     real(rkind) :: Yin_H2,Yin_O2,Yin_N2,Yout_H2O
     real(rkind) :: o2n2_ratio,h2o2_stoichiometric,h2o2_ratio

     !! Position and scale     
     flame_location = -0.25d0!zero!-0.27d0
     flame_thickness = 2.0d-4/L_char !! Scale thickness because position vectors are scaled...
   
     !! Determine inlet composition
     o2n2_ratio = one*molar_mass(2)/(3.76d0*molar_mass(9))
     h2o2_stoichiometric = two*molar_mass(1)/(one*molar_mass(2))
     h2o2_ratio = h2o2_stoichiometric*phi_in

     Yin_O2 = one/(one + h2o2_ratio + one/o2n2_ratio)
     Yin_H2 = Yin_O2*h2o2_ratio
     Yin_N2 = one - Yin_H2 - Yin_O2     
     
     !! Outlet composition
     Yout_H2O = one - Yin_N2     

     !! Temperatures
     T_reactants = T_ref
     T_products = 2.366d3
     
     !! Pressure through flame is reference pressure
     P_flame = p_ref

     !! Inflow density
     ro_inflow = p_flame/(Rgas_universal*T_reactants* &
             (Yin_H2*one_over_molar_mass(1) + Yin_O2*one_over_molar_mass(2) + Yin_N2*one_over_molar_mass(9)))      
    
     !! Inflow speed
     u_reactants = u_char
     
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Error function based progress variable
        c = half*(one + erf((x-flame_location)/flame_thickness))
!        c = exp(-((x-flame_location)/flame_thickness)**two - (y/flame_thickness)**two)        
!        c = exp(-((x-flame_location)/flame_thickness)**two - ((y-0.05d0)/flame_thickness)**two)    &
!          + exp(-((x-flame_location)/flame_thickness)**two - ((y+0.05d0)/flame_thickness)**two)        
!        c = exp(-((x-flame_location)/flame_thickness)**two)
        
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
        ro(i) = P_flame/(Rmix_local*T(i))
        
        !! Velocity
        u(i) = u_reactants*ro_inflow/ro(i)
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
              T(i) = T_reactants
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
!              u(i)=u_char
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
           end if
           T_bound(j) = T(i)
        end do
     end if
 
  
  
     return
  end subroutine make_1d_21step_flame
!! ------------------------------------------------------------------------------------------------  
  subroutine make_ignition_hotspot
     !! Make a hot spot for ignition (hard-coded to fit 9spec H2 combustion at present
     integer(ikind) :: i,ispec,j
     real(rkind) :: flame_location,flame_thickness
     real(rkind) :: T_hot
     real(rkind) :: P_local,c,Rmix_local,x,y,z

     !! Position and scale     
     flame_location = -0.27d0!zero!-0.27d0
     flame_thickness = 2.0d-4/L_char !! Scale thickness because position vectors are scaled...
   
     !! Temperature peak
     T_hot = 2.5d3
     
     !$omp parallel do private(x,y,z,c,Rmix_local,ispec,P_local)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
        
        !! Gaussian hotspot
        c = exp(-((x-flame_location)/flame_thickness)**two - (y/flame_thickness)**two)        
!        c = exp(-((x-flame_location)/flame_thickness)**two)
        
               
        !! Local mixture gas constant
        Rmix_local = zero
        do ispec=1,nspec
           Rmix_local = Rmix_local + Yspec(i,ispec)*one_over_molar_mass(ispec)
        end do
        Rmix_local = Rmix_local*Rgas_universal
        
        !! Local pressure
        P_local = ro(i)*Rmix_local*T(i)
        
        !! Modify temperature
        T(i) = T(i)*(one-c) + c*T_hot
        
        !! Modify density so pressure is unaffected
        ro(i) = P_local/(Rmix_local*T(i))
        
        !! composition and velocity are not modified
                        
     end do
     !$omp end parallel do 
  
  
     return
  end subroutine make_ignition_hotspot  
!! ------------------------------------------------------------------------------------------------  
  subroutine make_1d_25step_flame
     !! Make a 25 step (16 species) CH4-AIR flame. Initial profiles are erf(x').
     integer(ikind) :: i,ispec,j
     real(rkind) :: flame_location,flame_thickness
     real(rkind) :: T_reactants,T_products
     real(rkind) :: P_flame,c,u_reactants,Rmix_local,x,y,z,ro_inflow
     real(rkind) :: Yin_CH4,Yin_O2,Yin_N2,Yout_H2O,Yout_CO2
     real(rkind) :: o2n2_ratio,ch4o2_ratio,ch4o2_stoichiometric,co2h2o_ratio

     !! Position and scale     
     flame_location = zero
     flame_thickness = 5.0d-4/L_char !! Scale thickness because position vectors are scaled...

     !! Determine inlet composition
     o2n2_ratio = one*molar_mass(2)/(3.76d0*molar_mass(16))
     ch4o2_stoichiometric = one*molar_mass(1)/(two*molar_mass(2))
     ch4o2_ratio = ch4o2_stoichiometric*phi_in
          
     Yin_O2 = one/(one + ch4o2_ratio + one/o2n2_ratio)
     Yin_CH4 = Yin_O2*ch4o2_ratio
     Yin_N2 = one - Yin_CH4 - Yin_O2

     !! Outlet composition 
     co2h2o_ratio = one*molar_mass(3)/(two*molar_mass(4))  !! Assume complete combustion??
     Yout_CO2 = (one - Yin_N2)/(one + one/co2h2o_ratio)
     Yout_H2O = Yout_CO2/co2h2o_ratio            

     !! Temperatures
     T_reactants = T_ref
     T_products = 2.3d3     
     
     !! Pressure through flame 
     P_flame = p_ref
 
     !! inflow density
     ro_inflow = p_flame/(Rgas_universal*T_reactants* &
             (Yin_CH4*one_over_molar_mass(1) + Yin_O2*one_over_molar_mass(2) + Yin_N2*one_over_molar_mass(16)))

     
     !! Inflow speed
     u_reactants = u_char
     
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
        ro(i) = P_flame/(Rmix_local*T(i))
        
        !! Velocity
        u(i) = u_reactants*ro_inflow/ro(i)
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
           end if                 
           if(node_type(i).eq.1) then !! inflow initial conditions
              u(i)=u_char
           end if
           if(node_type(i).eq.2) then !! outflow initial conditions
           end if
           T_bound(j) = T(i)           
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
     real(rkind) :: x,y,z,tmp,tmpro
     
     !! Values within domain
     !$OMP PARALLEL DO PRIVATE(x,y,z,tmp,ispec)
     do i=1,npfb
        x = rp(i,1);y=rp(i,2);z=rp(i,3)
!        u(i) = -cos(two*pi*x)*sin(two*pi*y)*cos(two*pi*z/Lz)!*oosqrt2
!        v(i) = sin(two*pi*x)*cos(two*pi*y)*cos(two*pi*z/Lz)    !!c c
        !! old
!        tmp = -half*half*(cos(two*x) + cos(two*y))/csq  !! Modify for not(isoT)

        !! TG 3D Re1600 as in Cant 2022, Sandam 2017 etc (ish)
        u(i) = -cos(x)*sin(y)*cos(z)!*oosqrt2
        v(i) = sin(x)*cos(y)*cos(z)    !!c c
        w(i) = zero!u(i);u(i)=zero
        tmp = Rgas_universal*one_over_molar_mass(1)*T_ref  !! RT0
        tmp = rho_char*U_char*U_char/tmp   !! roUU/RT0
        tmp = -tmp*(one/16.0d0)*(cos(two*x)+cos(two*y))*(two+cos(two*z))        
        T(i) = T_ref

        ro(i) = rho_char + tmp

#ifdef ms    
        tmp = one - half*(one + erf(5.0d0*x))
        Yspec(i,1) = tmp
#else
        Yspec(i,1) = one
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
  subroutine load_restart_file
     !! Load initial conditions from a dump file
     integer(ikind) :: k,i,j
     real(rkind) :: tmp,tmpro
     character(70) :: fname  

#ifdef mp
     k=10000+iproc
#else
     k=10000
#endif     

     !! Construct the file name:
     write(fname,'(A16,I5)') './restart/layer_',k

     !! Open the file
     open(14,file=fname)
     read(14,*) k
     if(k.ne.npfb) write(6,*) "WARNING, expecting problem in restart"

     !! Load the initial conditions
     do i=1,npfb
#ifdef dim3
        read(14,*) tmp,tmp,tmp,tmp,h(i),k,tmpro,u(i),v(i),w(i),tmp,T(i),hrr(i),Yspec(i,1:nspec)
#else
        read(14,*) tmp,tmp,tmp,h(i),k,tmpro,u(i),v(i),tmp,T(i),hrr(i),Yspec(i,1:nspec)
#endif        
        if(k.ne.node_type(i)) then
           write(6,*) "ERROR: Problem in restart file. STOPPING."
#ifdef mp
           call MPI_Abort(MPI_COMM_WORLD, k, ierror)
#else
           stop
#endif
        end if
        ro(i) = tmpro
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
