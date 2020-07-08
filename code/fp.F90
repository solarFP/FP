MODULE FP
! This module encapsulates the Fokker-Plank solver as described in Allred et al (2020). 
! Module Variables:
!   fpinput:	fpinputtype class which stores several input options that control how FP is run.
!   fpoutput:	fpoutputtype class which stores the distribution function, f, after FP has finished.
!   fp_restart:	True if FP has been restarted from a previous solution. False, otherwise
  use fpinopts
  implicit none
  public :: FP_Solve

  TYPE(fpinputtype) fpinput
  TYPE(fpoutputtype) fpoutput
  logical, private :: fp_restart = .false.
  contains
  SUBROUTINE FP_Solve()
    implicit none

    call FP_set_options()
    call FP_set_beam_particle()
    call FP_set_sizes()
    call FP_set_loop_atm()
    call FP_set_grids()
    call FP_InjectPowerLaw()
    call FP_SetupMatrix()
    call FP_InitSolution()
    call FP_Output_Init()
    call FP_DoSolve()
    call FP_Output()
    call FP_End()
  ENDSUBROUTINE

  SUBROUTINE FP_RestartFromFile(infile, maxiter, tolres, toldiff, implicit_theta)
    use writer, only : readout
    implicit none
    integer maxiter
    double precision tolres, toldiff, implicit_theta
    character(len = 256) infile

    fp_restart = .true.
    call readout(infile, fpinput,fpoutput)
    ! replace read-in options with those specified in this routine
    fpinput%maxiter = maxiter; fpinput%tolres = tolres; fpinput%toldiff = toldiff; fpinput%implicit_theta = implicit_theta
    call FP_Solve()
  ENDSUBROUTINE

  SUBROUTINE FP_set_options()
    use options
    implicit none

    Emin = fpinput%Emin; Emax = fpinput%Emax; inc_relativity = fpinput%inc_relativity; inc_CC= fpinput%inc_CC 
    inc_synchro = fpinput%inc_synchro; inc_magmirror = fpinput%inc_magmirror
    inc_RC = fpinput%inc_RC; reflecttop = fpinput%reflecttop; oneD = fpinput%oneD
    reflectbottom = fpinput%reflectbottom;  maxiter = fpinput%maxiter; writeoutput = fpinput%writeoutput; tolres =  fpinput%tolres 
    toldiff =  fpinput%toldiff; implicit_theta =  fpinput%implicit_theta; outfile = fpinput%outfile; 
    resist_fact =  fpinput%resist_fact
    if (oneD) then 
      fpinput%nmu = 1
      fpinput%patype = 0
    endif
  ENDSUBROUTINE
 
  SUBROUTINE FP_set_beam_particle()
    use beam
    implicit none
    mbeam =  fpinput%mbeam; zbeam =  fpinput%zbeam
  ENDSUBROUTINE

  SUBROUTINE FP_set_sizes()
    use grid
    use loop
    implicit none
    nE_all =  fpinput%nE; nmu_all =  fpinput%nmu; nz_all =  fpinput%nz
    Nion =  fpinput%Nion; Nneutral =  fpinput%Nneutral
  ENDSUBROUTINE

  SUBROUTINE FP_set_loop_atm()
    use grid
    use loop
    implicit none

    call init_loop()
    zin_all =  fpinput%zin; tg_all =  fpinput%tg; bfield_all =  fpinput%bfield; dni_all =  fpinput%dni; dnn_all =  fpinput%dnn; 
    mion =  fpinput%mion; Zion =  fpinput%Zion; Zn =  fpinput%Zn; Enion =  fpinput%Enion
  ENDSUBROUTINE
  
  SUBROUTINE FP_set_grids()
    use options
    use parmpi
    use grid
    use loop
    use beam, only : mbeam
    implicit none
    call init_par()
    call setup_par_grids()
    call getmyatm()
    ! convert Emax and Emin to eV
    call init_grid(mbeam,Emax*1d3,Emin*1d3, fpinput%Ecut*1d3, inc_relativity, zin_all,fulloffset,oneD)
  ENDSUBROUTINE

  SUBROUTINE FP_InjectPowerLaw()
    use beam
    implicit none
    call init_beam()
    call InjectPowerLaw(fpinput%Ecut,  fpinput%dlt,  fpinput%eflux,  fpinput%patype,  fpinput%pasigma)
    call updateBC()
  ENDSUBROUTINE
  
  SUBROUTINE FP_SetupMatrix()
    use matrix
    implicit none
    call SetupMatrix()
  ENDSUBROUTINE

  SUBROUTINE FP_InitSolution()
    use matrix
    implicit none
    if (fp_restart) then
      call FP_RestartSoln()
    else
      call Analytic_Soln()
    endif
  ENDSUBROUTINE

  SUBROUTINE FP_RestartSoln()
    use const
    use grid
    use beam
    use parmpi
    implicit none
    double precision, allocatable :: fE(:,:,:) ! f_E is stored in output. Need to convert to f_p
    integer i,j,k

    allocate(fE(nE,nmu,nz))
    call par_scatterall(fE, fpoutput%f)
    !convert from fE to fp. fp = fE * esvol/msvol
    forall (k = 1:nz, j = 1:nmu, i = 1:nE)
      f(i,j,k) = fe(i,j,k)*(E(i+1)-E(i))*(mbeam/1d3)*(mu(j)-mu(j+1))*2*pi / msvol(i,j) ! #particles/cm^3 /d^3p
    endforall
    deallocate(fE)
    call UpdateBC()
  ENDSUBROUTINE

  SUBROUTINE FP_Output_Init()
    use writer
    use options
    implicit none
    if (writeoutput) call init_write(fpinput,fpoutput) 
  ENDSUBROUTINE

  SUBROUTINE FP_DoSolve()
    use options
    use beam
    use matrix
    use parmpi
    implicit none
    double precision r, diff
    integer ts
    logical conv

    call CalcNflux()
    r= residual()
    do ts = 1, maxiter
      if (.not.oneD) call SweepMu()
      call SweepP()
      call CalcNflux()
      r = residual()
      diff= CalcDiff()
      if (myid.eq.0) write(*,'(I5,A,2ES12.3)') ts, ': (l2 norm, Ave. change)', r,diff
      conv = IsConverged(r,diff)
      if (conv) exit ! Converged so break out of loop
    enddo
    if (conv.and.myid.eq.0) then
      write(*,*) 'Converged in',ts,'iterations.'
    else if (myid.eq.0) then
      write(*,'(A,I0,A,2ES10.3)') 'After ',maxiter, ' iterations not yet reached requested tolerance of',tolres, toldiff  
    endif
  ENDSUBROUTINE
  SUBROUTINE FP_Output()
    use options
    use writer
    implicit none
    if (writeoutput) then
      call writeout(fpoutput)
      call closeout()
    endif
  ENDSUBROUTINE
  SUBROUTINE FP_End()
    use options
    use beam
    use matrix
    use parmpi
    use loop
    use grid
    call deallocate_matrix()
    call deallocate_beam()
    call deallocate_loop()
    call deallocate_grid()
    call par_end()
  ENDSUBROUTINE
ENDMODULE
