MODULE fpinopts
! This module defines fpinputtype and fpoutputtype
! fpinputtype encapsulate the input parameters and options which are communicated to FP from the calling procedure 
  TYPE fpinputtype 
    integer nE, nmu, nz, nion, nneutral, maxiter, patype
    double precision Emin, Emax, tolres, toldiff, implicit_theta, mbeam, zbeam, Ecut, dlt, eflux, pasigma, resist_fact
    double precision, allocatable::  zin(:), tg(:), bfield(:), dni(:,:), dnn(:, :), mion(:), Zion(:), Zn(:), Enion(:)
    logical inc_relativity, inc_CC, inc_synchro, inc_magmirror, inc_RC, oneD, reflecttop, reflectbottom
  ENDTYPE
! fpoutputtype encapsulates the output produced by FP and communicated to the calling procedure
  TYPE fpoutputtype
    double precision, allocatable :: f(:,:,:), E(:), mu(:), z(:), heatrate(:), momrate(:),esvol(:,:)
    integer nE, nmu, nz
  ENDTYPE
ENDMODULE
