module common_parameter
  !! ----------------------------------------------------------------------------------------------
  !! SUNSET CODE: Scalable Unstructured Node-SET code for DNS.
  !! 
  !! Author             |Date             |Contributions
  !! --------------------------------------------------------------------------
  !! JRCK               |2019 onwards     |Main developer                     
  !!
  !! ----------------------------------------------------------------------------------------------
  !! This module contains parameters for the sunset code
  use kind_parameters
  implicit none 

  !! Numbers --------------------------------------------------------------------------------------
  real(rkind), parameter :: pi=3.141592653589793238462643383279502884197d0
  real(rkind), parameter :: pi4 = pi**4.0
  real(rkind), parameter :: oosqrt2 = 1.0d0/dsqrt(2.0d0)
  real(rkind), parameter :: fourthirds = 4.0d0/3.0d0
  real(rkind), parameter :: onethird = 1.0d0/3.0d0
  real(rkind), parameter :: twothirds = 2.0d0/3.0d0
  real(rkind), parameter :: zero = 0.0d0
  real(rkind), parameter :: one = 1.0d0
  real(rkind), parameter :: two = 2.0d0
  real(rkind), parameter :: half = 0.5d0
  real(rkind), parameter :: oosix = 1.0d0/6.0d0
  real(rkind), parameter :: verysmall = 1.0d-30
  
  !! Physical constants ---------------------------------------------------------------------------
  real(rkind), parameter :: Rgas_universal = 8.3142d3         !! Universal gas constant  

  !! (NODE-SET) Discretisation related parameters
  integer(ikind) ,parameter :: dims = 3
  real(rkind), parameter :: hovs = 2.7d0   !! stencil scale over discretisation scale (h/s)
  real(rkind), parameter :: hovs_bound = 2.4d0 !! as above, reduced near bounds for stability
  real(rkind), parameter :: ss = 2.0d0       !! Stencil size (radius, in multiples of h)

  !! Maximum allowable number of species
  integer(ikind), parameter :: nspec_max = 20     
  
  !! Runge Kutta coefficients ---------------------------------------------------------------------
  !! RK3(2)4[2R+]C Kennedy (2000) Appl. Num. Math. 35:177-219
  real(rkind),parameter :: rk3_4s_2r_a21 = 11847461282814.0d0/36547543011857.0d0
  real(rkind),parameter :: rk3_4s_2r_a32 = 3943225443063.0d0/7078155732230.0d0
  real(rkind),parameter :: rk3_4s_2r_a43 = -346793006927.0d0/4029903576067.0d0
  real(rkind),parameter :: rk3_4s_2r_b1 = 1017324711453.0d0/9774461848756.0d0
  real(rkind),parameter :: rk3_4s_2r_b2 = 8237718856693.0d0/13685301971492.0d0
  real(rkind),parameter :: rk3_4s_2r_b3 = 57731312506979.0d0/19404895981398.0d0
  real(rkind),parameter :: rk3_4s_2r_b4 = -101169746363290.0d0/37734290219643.0d0
  real(rkind),parameter :: rk3_4s_2r_bh1 = 15763415370699.0d0/46270243929542.0d0
  real(rkind),parameter :: rk3_4s_2r_bh2 = 514528521746.0d0/5659431552419.0d0
  real(rkind),parameter :: rk3_4s_2r_bh3 = 27030193851939.0d0/9429696342944.0d0
  real(rkind),parameter :: rk3_4s_2r_bh4 = -69544964788955.0d0/30262026368149.0d0
  real(rkind),dimension(3),parameter :: rk3_4s_2r_a=(/rk3_4s_2r_a21,rk3_4s_2r_a32,rk3_4s_2r_a43/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_b=(/rk3_4s_2r_b1,rk3_4s_2r_b2,rk3_4s_2r_b3,rk3_4s_2r_b4/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_bh=(/rk3_4s_2r_bh1,rk3_4s_2r_bh2,rk3_4s_2r_bh3,rk3_4s_2r_bh4/)
  real(rkind),dimension(4),parameter :: rk3_4s_2r_bmbh = rk3_4s_2r_b - rk3_4s_2r_bh

  !! Normalisation constants and parameters for PID error estimators (OK for combustion at standard P,T)
  real(rkind), parameter :: elnro_norm = 1.0d-10
  real(rkind), parameter :: eu_norm = 1.0d-2 
  real(rkind), parameter :: ev_norm = 1.0d-2 !! 1d-6 might be necessary for steadier flows.
  real(rkind), parameter :: ew_norm = 1.0d-2 !! 1d-2 is better for dealing with initial shocks
  real(rkind), parameter :: eroE_norm = 1.0d-2 
  real(rkind), parameter :: eY_norm = 1.0d-10       
  real(rkind), parameter :: pid_tol = 1.0d-4        !! Error tolerance
  real(rkind), parameter :: pid_a=0.7d0/two  !! P-coefficient   ! 0.7
  real(rkind), parameter :: pid_b=0.4d0/two  !! I-coefficient   ! 0.4
  real(rkind), parameter :: pid_c=0.1d0/two  !! D-coefficient   ! 0.1
  real(rkind), parameter :: pid_k=0.9d0      !! kappa coefficient !0.9
  
  !! Finite difference stencil sizes
#define FDORDER 8              
  !! Size of Stencil
#if FDORDER==4
  integer(ikind),parameter :: ij_count_fd = 5
#elif FDORDER==6
  integer(ikind),parameter :: ij_count_fd = 7
#elif FDORDER==8
  integer(ikind),parameter :: ij_count_fd = 9
#elif FDORDER==10
  integer(ikind),parameter :: ij_count_fd = 11
#elif FDORDER==12
  integer(ikind),parameter :: ij_count_fd = 13
#endif      
  
end module common_parameter
