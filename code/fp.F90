MODULE FP
! This module encapsulates the Fokker-Plank solver as described in Allred et al (2020). 
  use fpinopts
  implicit none
  public :: FP_Solve

  contains
  SUBROUTINE FP_Solve(fpinput,fpoutput, restart)
    implicit none
    TYPE(fpinputtype) fpinput
    TYPE(fpoutputtype) fpoutput
    logical restart
    call FP_set_options(fpinput)
    call FP_set_beam_particle(fpinput)
    call FP_set_sizes(fpinput)
    call FP_set_loop_atm(fpinput)
    call FP_set_grids(fpinput)
    call FP_InjectPowerLaw(fpinput)
    call FP_SetupMatrix()
    call FP_InitSolution(fpoutput,restart)
    call FP_Output_Init(fpoutput)
    call FP_DoSolve()
    call FP_Output(fpoutput)
    call FP_End()
  ENDSUBROUTINE

  SUBROUTINE FP_set_options(fpinput)
    use options
    implicit none
    TYPE(fpinputtype) fpinput

    Emin = fpinput%Emin; Emax = fpinput%Emax; inc_relativity = fpinput%inc_relativity; inc_CC= fpinput%inc_CC 
    inc_synchro = fpinput%inc_synchro; inc_magmirror = fpinput%inc_magmirror
    inc_RC = fpinput%inc_RC; reflecttop = fpinput%reflecttop; oneD = fpinput%oneD
    reflectbottom = fpinput%reflectbottom;  maxiter = fpinput%maxiter; tolres =  fpinput%tolres 
    toldiff =  fpinput%toldiff; implicit_theta =  fpinput%implicit_theta 
    resist_fact =  fpinput%resist_fact
    if (oneD) then 
      fpinput%nmu = 1
      fpinput%patype = 0
    endif
  ENDSUBROUTINE
 
  SUBROUTINE FP_set_beam_particle(fpinput)
    use beam
    implicit none
    TYPE(fpinputtype) fpinput
    mbeam =  fpinput%mbeam; zbeam =  fpinput%zbeam
  ENDSUBROUTINE

  SUBROUTINE FP_set_sizes(fpinput)
    use grid
    use loop
    implicit none
    TYPE(fpinputtype) fpinput
    nE_all =  fpinput%nE; nmu_all =  fpinput%nmu; nz_all =  fpinput%nz
    Nion =  fpinput%Nion; Nneutral =  fpinput%Nneutral
  ENDSUBROUTINE

  SUBROUTINE FP_set_loop_atm(fpinput)
    use grid
    use loop
    implicit none
    TYPE(fpinputtype) fpinput
    call init_loop()
    zin_all =  fpinput%zin; tg_all =  fpinput%tg; bfield_all =  fpinput%bfield; dni_all =  fpinput%dni; dnn_all =  fpinput%dnn; 
    mion =  fpinput%mion; Zion =  fpinput%Zion; Zn =  fpinput%Zn; Enion =  fpinput%Enion
  ENDSUBROUTINE
  
  SUBROUTINE FP_set_grids(fpinput)
    use options
    use parmpi
    use grid
    use loop
    use beam, only : mbeam
    implicit none
    TYPE(fpinputtype) fpinput
    call init_par()
    call setup_par_grids()
    call getmyatm()
    ! convert Emax and Emin to eV
    call init_grid(mbeam,Emax*1d3,Emin*1d3, fpinput%Ecut*1d3, inc_relativity, zin_all,fulloffset,oneD)
  ENDSUBROUTINE

  SUBROUTINE FP_InjectPowerLaw(fpinput)
    use beam
    implicit none
    TYPE(fpinputtype) fpinput
    call init_beam()
    call InjectPowerLaw(fpinput%Ecut,  fpinput%dlt,  fpinput%eflux,  fpinput%patype,  fpinput%pasigma)
    call updateBC()
  ENDSUBROUTINE
  
  SUBROUTINE FP_SetupMatrix()
    use matrix
    implicit none
    call SetupMatrix()
  ENDSUBROUTINE

  SUBROUTINE FP_InitSolution(fpoutput,restart)
    use matrix
    implicit none
    logical restart
    TYPE(fpoutputtype) fpoutput
    if (restart) then
      call FP_RestartSoln(fpoutput)
    else
      call Analytic_Soln()
    endif
  ENDSUBROUTINE

  SUBROUTINE FP_RestartSoln(fpoutput)
    use const
    use grid
    use beam
    use parmpi
    implicit none
    TYPE(fpoutputtype) fpoutput
    double precision, allocatable :: fE(:,:,:) ! f_E is stored in output. Need to convert to f_p
    integer i,j,k

    allocate(fE(nE,nmu,nz))
    call par_scatterall(fE, fpoutput%f)
    !convert from fE to fp. fp = fE * esvol/msvol
    forall (k = 1:nz, j = 1:nmu, i = 1:nE)
      f(i,j,k) = fe(i,j,k)*esvol(i,j) / msvol(i,j) ! #particles/cm^3 /d^3p
    endforall
    deallocate(fE)
    call UpdateBC()
  ENDSUBROUTINE

  SUBROUTINE FP_Output_Init(fpoutput)
    use options
    use grid
    use beam
    use parmpi
    implicit none
    integer i,j
    TYPE(fpoutputtype) fpoutput

    fpoutput%nE = nE_all; fpoutput%nmu = nmu_all; fpoutput%nz = nz_all
    if (myid.eq.0) then
      if (.not.allocated(fpoutput%E)) allocate(fpoutput%E(nE_all+1))
      if (.not.allocated(fpoutput%mu)) allocate(fpoutput%mu(nmu_all+1))
      if (.not.allocated(fpoutput%z)) allocate(fpoutput%z(nz_all+1))
    else
      if (.not.allocated(fpoutput%E)) allocate(fpoutput%E(1))
      if (.not.allocated(fpoutput%mu)) allocate(fpoutput%mu(1))
      if (.not.allocated(fpoutput%z)) allocate(fpoutput%z(1))
    endif
    if (myplaneid(3).eq.0) then
      if (.not.allocated(fpoutput%esvol)) allocate(fpoutput%esvol(nE_all,nmu_all))
    else
      if (.not.allocated(fpoutput%esvol)) allocate(fpoutput%esvol(1,1))
    endif
    call par_GetGridAll(fpoutput%E,fpoutput%mu,fpoutput%z)
    call par_CollectMS(esvol,fpoutput%esvol)
    fpoutput%E = fpoutput%E * mbeam / 1d3 ! convert E to keV
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
  SUBROUTINE FP_Output(fpoutput)
    use options
    use parmpi
    use grid
    use matrix
    use beam
    implicit none
    TYPE(fpoutputtype) fpoutput
    double precision, allocatable :: fe(:,:,:)
    double precision qh(nz), momrate(nz)
    integer i,j,k

    if (.not.allocated(fpoutput%heatrate)) allocate(fpoutput%heatrate(nz_all))
    if (.not.allocated(fpoutput%momrate)) allocate(fpoutput%momrate(nz_all))

    call CalcHeatRate(qh,momrate)
    call par_Gather1D(qh,fpoutput%heatrate,nz,3)
    call par_Gather1D(momrate,fpoutput%momrate,nz,3)

    if (myid.eq.0) then
      if (.not.allocated(fpoutput%f)) allocate(fpoutput%f(nE_all,nmu_all,nz_all))
    else
      if (.not.allocated(fpoutput%f)) allocate(fpoutput%f(1,1,1))
    endif
    !convert from fp to fE. fE = fp * msvol / esvol
    allocate(fe(nE,nmu,nz))
    forall (k = 1:nz, j = 1:nmu, i = 1:nE)
      fe(i,j,k) = f(i,j,k)*msvol(i,j)/esvol(i,j) ! #particles/cm^3 /keV /sr
    endforall
    call par_collectall(fe,fpoutput%f)
    deallocate(fe)
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
