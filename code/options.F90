module options
! This module stores input options for internal use by FP
 double precision Emin, Emax, implicit_theta, tolres, toldiff ! Convergence tolerances
 integer maxiter
 logical inc_relativity, writeoutput, reflecttop, reflectbottom ! BC conditions
 logical inc_CC, inc_synchro, inc_magmirror, inc_RC, oneD  !INCLUDE FORCES 
 double precision  resist_fact
 character(len=256) outfile

end module 
