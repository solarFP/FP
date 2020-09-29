MODULE writer
! This module handles the writing (and reading in the case of restarting) to a binary HDF5 file.
! Module Variables:
!  file_id:	The HDF5 file id

  use hdf5readwrite
  implicit none
  integer(HID_T), private :: file_id       ! File identifier

  contains 
  SUBROUTINE init_write(outfile,fpin, fpout)
    use fpinopts
    use parmpi
    implicit none
    TYPE(fpinputtype) fpin
    TYPE(fpoutputtype) fpout
    integer i,j, error
    integer(HID_T) :: group_id
    character(len=256) outfile

    !if (myid.eq.0) then
      call h5open_f(error)
      call h5fcreate_f(outfile, H5F_ACC_TRUNC_F, file_id, error)
      call h5gcreate_f(file_id, 'inputparams', group_id, error)
      call writesattr(fpin%nE,'nE',group_id)
      call writesattr(fpin%nmu,'nmu',group_id)
      call writesattr(fpin%nz,'nz',group_id)
      call writesattr(fpin%maxiter,'maxiter',group_id)
      call writesattr(fpin%patype,'patype',group_id)
      call writesattr(fpin%inc_relativity,'inc_relativity',group_id)
      call writesattr(fpin%inc_CC,'inc_CC',group_id)
      call writesattr(fpin%inc_synchro,'inc_synchro',group_id)
      call writesattr(fpin%inc_magmirror,'inc_magmirror',group_id)
      call writesattr(fpin%inc_RC,'inc_RC',group_id)
      call writesattr(fpin%oneD,'oneD',group_id)
      call writesattr(fpin%reflecttop,'reflecttop',group_id)
      call writesattr(fpin%reflectbottom,'reflectbottom',group_id)
      call writesattr(fpin%Emin,'Emin',group_id)
      call writesattr(fpin%Emax,'Emax',group_id)
      call writesattr(fpin%tolres,'tolres',group_id)
      call writesattr(fpin%toldiff,'toldiff',group_id)
      call writesattr(fpin%implicit_theta,'implicit_theta',group_id)
      call writesattr(fpin%mbeam,'mbeam',group_id)
      call writesattr(fpin%zbeam,'zbeam',group_id)
      call writesattr(fpin%Ecut,'Ecut',group_id)
      call writesattr(fpin%dlt,'dlt',group_id)
      call writesattr(fpin%eflux,'eflux',group_id)
      call writesattr(fpin%pasigma,'pasigma',group_id)
      call writesattr(fpin%resist_fact,'resist_fact',group_id)
      call h5gclose_f(group_id, error)

      call h5gcreate_f(file_id, 'atm', group_id, error)
      call writesattr(fpin%nneutral,'nNeutral',group_id)
      call writesattr(fpin%nion,'nIon',group_id)
      call writearray(fpin%zin,'zin',group_id)
      call writearray(fpin%tg,'tg',group_id)
      call writearray(fpin%bfield,'bfield',group_id)
      call writearray(fpin%dni,'dni',group_id)
      call writearray(fpin%dnn,'dnn',group_id)
      call writearray(fpin%mion,'mion',group_id)
      call writearray(fpin%zion,'Zion',group_id)     
      call writearray(fpin%zn,'Zn',group_id)
      call writearray(fpin%Enion,'Enion',group_id)
      call h5gclose_f(group_id, error)

      call h5gopen_f(file_id,'/',group_id,error)
      call writearray(fpout%E,'E',group_id)
      call writearray(fpout%mu,'mu',group_id)
      call writearray(fpout%z,'z',group_id)
      call writearray(fpout%esvol,'esvol',group_id)
      call h5gclose_f(group_id, error)
    !endif
  ENDSUBROUTINE

  SUBROUTINE writeout(fpout)
    use fpinopts
    use hdf5readwrite
    implicit none
    integer ierr
    integer(HID_T) group_id
    TYPE(fpoutputtype) fpout

    !if (myid.eq.0) then
      call h5gopen_f(file_id,'/',group_id,ierr)
      call writearray(fpout%heatrate,'heatrate',group_id)
      call writearray(fpout%momrate,'momrate',group_id)
      call writearray(fpout%f,'f',group_id)
      call h5gclose_f(group_id,ierr)
      call h5fclose_f(file_id, ierr)
      call h5close_f(ierr)
    !endif
    !call mpi_barrier(cart_comm,ierr)
  ENDSUBROUTINE

  SUBROUTINE readout(infile, fpin, fpout,myid)
    use fpinopts
    !use mpi
    use hdf5readwrite
    implicit none
    TYPE(fpinputtype) fpin
    TYPE(fpoutputtype) fpout
    !double precision, allocatable :: fe(:,:,:)
    integer myid, ierr
    integer(HID_T) :: group_id
    character(len = 256) infile
 
!    call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)
    !Either have every processor read the file or have one do it and Bcast to the rest.
    !right now it seems easiest to have every processor do it
!    if (myid.eq.0) then 
      call h5open_f(ierr)
      call h5fopen_f(infile, H5F_ACC_RDONLY_F, file_id, ierr)
      call h5gopen_f(file_id, 'inputparams', group_id, ierr)
      call readsattr(fpin%nE,'nE',group_id)
      call readsattr(fpin%nmu,'nmu',group_id)
      call readsattr(fpin%nz,'nz',group_id)
      call readsattr(fpin%maxiter,'maxiter',group_id)
      call readsattr(fpin%patype,'patype',group_id)
      call readsattr(fpin%inc_relativity,'inc_relativity',group_id)
      call readsattr(fpin%inc_CC,'inc_CC',group_id)
      call readsattr(fpin%inc_synchro,'inc_synchro',group_id)
      call readsattr(fpin%inc_magmirror,'inc_magmirror',group_id)
      call readsattr(fpin%inc_RC,'inc_RC',group_id)
      call readsattr(fpin%oneD,'oneD',group_id)
      call readsattr(fpin%reflecttop,'reflecttop',group_id)
      call readsattr(fpin%reflectbottom,'reflectbottom',group_id)
      call readsattr(fpin%Emin,'Emin',group_id)
      call readsattr(fpin%Emax,'Emax',group_id)
      call readsattr(fpin%tolres,'tolres',group_id)
      call readsattr(fpin%toldiff,'toldiff',group_id)
      call readsattr(fpin%implicit_theta,'implicit_theta',group_id)
      call readsattr(fpin%mbeam,'mbeam',group_id)
      call readsattr(fpin%zbeam,'zbeam',group_id)
      call readsattr(fpin%Ecut,'Ecut',group_id)
      call readsattr(fpin%dlt,'dlt',group_id)
      call readsattr(fpin%eflux,'eflux',group_id)
      call readsattr(fpin%pasigma,'pasigma',group_id)
      call readsattr(fpin%resist_fact,'resist_fact',group_id)
      call h5gclose_f(group_id, ierr)  

      call h5gopen_f(file_id, 'atm', group_id, ierr)
      call readsattr(fpin%nneutral,'nNeutral',group_id)
      call readsattr(fpin%nion,'nIon',group_id)

      allocate(fpin%zin(fpin%nz), fpin%tg(fpin%nz), fpin%bfield(fpin%nz), fpin%dni(fpin%Nion,fpin%nz), & 
           fpin%dnn(fpin%Nneutral,fpin%nz), fpin%mion(fpin%Nion), fpin%Zion(fpin%Nion), fpin%Zn(fpin%Nneutral), &
           fpin%Enion(fpin%Nneutral))

      call readarray(fpin%zin,'zin',group_id)
      call readarray(fpin%tg,'tg',group_id)
      call readarray(fpin%bfield,'bfield',group_id)
      call readarray(fpin%dni,'dni',group_id)
      call readarray(fpin%dnn,'dnn',group_id)
      call readarray(fpin%mion,'mion',group_id)
      call readarray(fpin%zion,'Zion',group_id)
      call readarray(fpin%zn,'Zn',group_id)
      call readarray(fpin%Enion,'Enion',group_id)
      call h5gclose_f(group_id, ierr)
    !endif
    if (myid.eq.0) then
      allocate(fpout%f(fpin%nE,fpin%nmu,fpin%nz))
      call h5gopen_f(file_id, '/', group_id, ierr)
      call readarray(fpout%f,'f',group_id)
      call h5gclose_f(group_id, ierr)
    else
      allocate(fpout%f(1,1,1))
    endif
    call h5fclose_f(file_id,ierr)
    call h5close_f(ierr)
  ENDSUBROUTINE

ENDMODULE
