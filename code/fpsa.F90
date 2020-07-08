PROGRAM FPSA
! This is the standalone driver that calls FP.
! It does the following steps:
!   1. Initialize MPI
!   2. Read in input parameters from param.cnt
!   3. Read in the loop atmospheric stratification
!   4. Copies input parameters into the fpinput structure
!   5. Calls FP_solve() or if restarting then calls FP_RestartFromFile()
!   6. Deallocate variables and closes MPI 
  use mpi
  use fp
  implicit none
  integer nE, nmu, maxiter, patype, ierr, myid
  double precision Emin, Emax, tolres, toldiff, implicit_theta, mbeam, zbeam, Ecut, dlt, eflux, pasigma, resist_fact
  !double precision, allocatable :: zin(:), tg(:), bfield(:), dni(:,:), dnn(:,:), mion(:), Zion(:), Zn(:), Enion(:)
  logical inc_CC, inc_synchro, inc_magmirror, inc_RC, oneD, reflecttop, reflectbottom, &
     writeoutput, inc_relativity, restart
  character(len=256) outfile, paramfile, atmfile
  
  namelist /control/ nE, nmu, Emin, Emax, inc_relativity, inc_CC, inc_synchro,inc_magmirror, inc_RC, oneD, &
                     reflecttop, reflectbottom, maxiter, writeoutput,tolres, toldiff, implicit_theta, atmfile, outfile, mbeam, &
                     zbeam, Ecut, dlt, eflux,patype, pasigma, resist_fact, restart
  call MPI_Init(ierr)
  call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)
  if (command_argument_count() .gt. 0) then
    call get_command_argument(1,paramfile)
  else
    paramfile = 'param.cnt'
  endif
  open(unit=80, file=paramfile, status='old')
  read(80,control)
  close(unit=80)
  
  if (.not.restart) then
    call read_atm(atmfile, fpinput%nz, fpinput%Nion, fpinput%Nneutral, fpinput%zin,fpinput%tg, fpinput%bfield, fpinput%dni, &
                  fpinput%dnn, fpinput%mion, fpinput%Zion, fpinput%Zn, fpinput%Enion, myid)
    fpinput%nE = nE; fpinput%nmu = nmu; fpinput%Emin = Emin; fpinput%Emax = Emax; fpinput%mbeam = mbeam; fpinput%Zbeam = Zbeam
    fpinput%inc_relativity = inc_relativity; fpinput%inc_CC = inc_CC 
    fpinput%inc_synchro = inc_synchro; fpinput%inc_magmirror = inc_magmirror; fpinput%inc_RC = inc_RC
    fpinput%oneD = oneD; fpinput%reflecttop = reflecttop; fpinput%reflectbottom = reflectbottom; fpinput%maxiter = maxiter
    fpinput%writeoutput = writeoutput; fpinput%tolres = tolres; fpinput%toldiff = toldiff; fpinput%implicit_theta = implicit_theta
    fpinput%outfile = outfile; fpinput%Ecut = Ecut; fpinput%dlt = dlt; fpinput%eflux = eflux; fpinput%patype = patype
    fpinput%pasigma = pasigma; fpinput%resist_fact = resist_fact
    call FP_Solve()
  else
    call FP_RestartFromFile(outfile, maxiter, tolres, toldiff, implicit_theta)
  endif
  deallocate(fpinput%zin,fpinput%tg,fpinput%bfield,fpinput%dni,fpinput%dnn,fpinput%mion,fpinput%Zion, fpinput%Zn, fpinput%Enion)
  call MPI_Finalize(ierr)

contains
SUBROUTINE read_atm(atmfile, nz, Nion, Nneutral, zin,tg,bfield,dni, dnn, mion, Zion, Zn, Enion, myid)
  use mpi
  implicit none
  double precision,allocatable :: zin(:), tg(:), bfield(:), dni(:,:), dnn(:,:), mion(:), Zion(:), Zn(:), Enion(:)
  integer nz, Nion, Nneutral, ns(3),ierr, myid
  character (len=256) atmfile
  if (myid.eq.0) then
    open(unit=20, file = atmfile,form = 'unformatted',status='old')
    read(20) nz, Nion, Nneutral
    ns = (/nz, Nion, Nneutral/)
  endif
  call MPI_Bcast(ns,3,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
  nz = ns(1); Nion = ns(2); Nneutral = ns(3)
  allocate(zin(nz), tg(nz), bfield(nz), dni(Nion,nz), dnn(Nneutral,nz), &
           mion(Nion), Zion(Nion), Zn(Nneutral), Enion(Nneutral))
  if (myid.eq.0) then
    read(20) zin,tg, bfield,dni, dnn, mion, Zion, Zn, Enion
    close(unit=20)
  endif
  call MPI_Bcast(zin,nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(tg,nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(bfield,nz,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(dni,nz*Nion,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(dnn,nz*Nneutral,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(mion,Nion,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Zion,Nion,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Zn,Nneutral,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
  call MPI_Bcast(Enion,Nneutral,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
ENDSUBROUTINE
END
