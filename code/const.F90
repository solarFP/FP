module const
! Contains several useful constants
  implicit none
  double precision, parameter :: me       = 0.5109989461d6      ! (eV) electron mass
  double precision, parameter :: mp       = 938.2720881629d6    ! (eV) proton mass
  double precision, parameter :: clight   = 2.99792458d10       ! (cm/s) speed of light
  double precision, parameter :: kb       = 8.617330350d-5      ! (eV/K) Boltzmann constant
  double precision, parameter :: pi       = acos(-1d0)!3.1415926535897931  ! pi
  double precision, parameter :: sqrtpi2  = 2d0/sqrt(pi)        ! 2/sqrt(pi)
  double precision, parameter :: e2       = 1.4399764d-7        ! (eV cm) (e charge)^2
  double precision, parameter :: pi4e4    = 4d0*pi*e2*e2        ! (eV^2 cm^2) 4 pi e^4
  double precision, parameter :: ergperev = 1.60217653d-12      ! (erg / eV) ;convert eV to erg
  double precision, parameter :: alpha2   = 0.0072973525693**2  ! (fine stucture constant) **2
  double precision, parameter :: hbc2     = 6.1992097e-05       ! hbar c /2 in units of eV cm
  double precision, parameter :: sqrtpie2 = sqrt(pi/e2)         ! sqrt( pi / e2)
END module 
