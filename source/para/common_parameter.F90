module common_parameter
  use kind_parameters
  implicit none 

  !! Discretisation related parameters
  integer(ikind) ,parameter :: dims = 3

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
  
  !! Runge Kutta coefficients ---------------------------------------------------------------------
  !! classic RK4 (from Wikipedia...)  
  real(rkind),dimension(3),parameter :: rk3_2N_a = (/ 0.0d0, -2.0d0/3.0d0, -1.0d0 /)
  real(rkind),dimension(3),parameter :: rk3_2N_b = (/ 1.0d0/3.0d0, 1.0d0, 0.5d0 /)
  real(rkind),dimension(3),parameter :: rk3_2N_c = (/ 0.0d0, 1.0d0/3.0d0, 2.0d0/3.0d0 /)  
  
  !! RK3[2N] Symmetric (from PENCIL CODE)
  real(rkind),dimension(4),parameter :: rk4_a = (/ 1.0d0,0.5d0,0.5d0,1.0d0 /)
  real(rkind),dimension(4),parameter :: rk4_b = (/ 1.0d0/6.0d0, 1.0d0/3.0d0, 1.0d0/3.0d0, 1.0d0/6.0d0 /)
  real(rkind),dimension(4),parameter :: rk4_c = (/ 0.0d0,0.5d0,0.5d0,1.0d0 /)
  
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
   
  !! SIMULATION PARAMETERS ========================================================================
  !! Primary domain parameters (i.e. those we can specify) ----------------------------------------
  real(rkind), parameter :: L_char = half*half*half*half*half !! Characteristic lengthscale
  real(rkind), dimension(dims), parameter :: grav = (/zero,zero,zero/) !! Gravity  
  
  !! Primary physical fluid properties ------------------------------------------------------------
  real(rkind), parameter :: visc0 = 1.8d-5
  real(rkind), parameter :: rho_char = one
  real(rkind), parameter :: gammagas = 1.4d0
  real(rkind), parameter :: Rs0 = 287.058d0 !! Reference specific gas constant  
  
  !! Primary dimensionless groups -----------------------------------------------------------------
  real(rkind), parameter :: Re = 1000.0d0 !! Reynolds number
  real(rkind), parameter :: Sc = one !! Schmidt number (higher Sc reduces importance of MD)  
  real(rkind), parameter :: Pr = one !! Prandtl number
  real(rkind), parameter :: Ma = 0.0016d0 !! Mach number  
  
  !! Secondary domain parameters ------------------------------------------------------------------
  real(rkind), parameter :: Lz = half*L_char  !! Domain length-scale in third dimension
  real(rkind), parameter :: u_char = Re*visc0/L_char/rho_char  !! Characteristic velocity
  real(rkind), parameter :: u_inflow = u_char  !! Inflow velocity (occasionally differs from u_char)
  real(rkind), parameter :: Time_char= L_char/u_char            !! Characteristic time-scale
    
  !! Secondary fluid parameters -------------------------------------------------------------------
  real(rkind), parameter :: gammagasm1 = gammagas - one
  real(rkind), parameter :: csq = (u_char/Ma)**two  !! Reference sound speed   
  real(rkind), parameter :: T0 = csq/gammagas/Rs0   !! Reference temperature 
  real(rkind), parameter :: MD = visc0/rho_char/Sc !! Mass/molecular diffusivity   
  real(rkind), parameter :: lambda_cond = Rs0*gammagas*visc0/gammagasm1/Pr  
#ifdef isoT
  real(rkind), parameter :: p_infinity = csq    !! Reference pressure
#else
  real(rkind), parameter :: p_infinity = rho_char*Rs0*T0
#endif


  !! Multi-species parameters ---------------------------------------------------------------------
#ifdef ms
  integer(ikind), parameter :: nspec = 2
#else
  integer(ikind), parameter :: nspec = 1
#endif  
   

end module common_parameter
