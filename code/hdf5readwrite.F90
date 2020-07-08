!Setup several subroutines to manage reading and writing HDF5 files.
SUBROUTINE writesattr_int(val,nme,id)
  use hdf5
  implicit none
  integer, intent(in) :: val
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) aspace_id, attr_id
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  integer error

  dims(1) = 1
  call h5screate_f(H5S_SCALAR_F, aspace_id, error)
  call h5acreate_f(id, nme, H5T_NATIVE_INTEGER, aspace_id, attr_id, error)
  call h5awrite_f(attr_id, H5T_NATIVE_INTEGER, (/val/), dims, error)
  call h5aclose_f(attr_id, error)
  call h5sclose_f(aspace_id, error)
ENDSUBROUTINE
SUBROUTINE writesattr_logical(val,nme,id)
  use hdf5
  implicit none
  logical, intent(in) :: val
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  call writesattr_int(merge(1,0,val),nme,id)
ENDSUBROUTINE
SUBROUTINE writesattr_double(val,nme,id)
  use hdf5
  implicit none
  double precision, intent(in) :: val
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) aspace_id, attr_id
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  integer error

  dims(1) = 1
  call h5screate_f(H5S_SCALAR_F, aspace_id, error)
  call h5acreate_f(id, nme, H5T_NATIVE_DOUBLE, aspace_id, attr_id, error)
  call h5awrite_f(attr_id, H5T_NATIVE_DOUBLE, (/val/), dims, error)
  call h5aclose_f(attr_id, error)
  call h5sclose_f(aspace_id, error)
ENDSUBROUTINE
SUBROUTINE readsattr_int(val,nme,id)
  use hdf5
  implicit none
  integer, intent(inout) :: val
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) attr_id
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  integer error, aval(1)

  dims(1) = 1
  call h5aopen_f(id, nme, attr_id, error)
  call h5aread_f(attr_id, H5T_NATIVE_INTEGER, aval, dims,error)
  call h5aclose_f(attr_id, error)
  val = aval(1)
ENDSUBROUTINE
SUBROUTINE readsattr_logical(val,nme,id)
  use hdf5
  implicit none
  logical, intent(inout) :: val
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer ival
  call readsattr_int(ival,nme,id)
  if (ival .eq.1) then 
    val = .true.
  else
    val = .false.
  endif
ENDSUBROUTINE
SUBROUTINE readsattr_double(val,nme,id)
  use hdf5
  implicit none
  double precision, intent(inout) :: val
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) attr_id
  INTEGER(HSIZE_T), DIMENSION(1) :: dims
  integer error
  double precision aval(1)

  dims(1) = 1
  call h5aopen_f(id, nme, attr_id, error)
  call h5aread_f(attr_id, H5T_NATIVE_DOUBLE, aval, dims,error)
  call h5aclose_f(attr_id, error)
  val = aval(1)
ENDSUBROUTINE

SUBROUTINE writearray_1d_double(arr, nme, id)
  use hdf5
  implicit none
  double precision, intent(in) :: arr(:)
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) dspace_id, dataset_id
  INTEGER(HSIZE_T) :: dims(1)
  integer error

  dims(1) = size(arr)
  call h5screate_simple_f(1, dims, dspace_id, error)
  CALL h5dcreate_f(id, nme, H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
  CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, arr, dims, error)
  CALL h5sclose_f(dspace_id, error)
  CALL h5dclose_f(dataset_id, error)
ENDSUBROUTINE
SUBROUTINE writearray_2d_double(arr, nme, id)
  use hdf5
  implicit none
  double precision, intent(in) :: arr(:,:)
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) dspace_id, dataset_id
  INTEGER(HSIZE_T) :: dims(2)
  integer error

  dims = shape(arr)
  call h5screate_simple_f(2, dims, dspace_id, error)
  CALL h5dcreate_f(id, nme, H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
  CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, arr, dims, error)
  CALL h5sclose_f(dspace_id, error)
  CALL h5dclose_f(dataset_id, error)
ENDSUBROUTINE
SUBROUTINE writearray_3d_double(arr, nme, id)
  use hdf5
  implicit none
  double precision, intent(in) :: arr(:,:,:)
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) dspace_id, dataset_id
  INTEGER(HSIZE_T) :: dims(3)
  integer error

  dims = shape(arr)
  call h5screate_simple_f(3, dims, dspace_id, error)
  CALL h5dcreate_f(id, nme, H5T_NATIVE_DOUBLE, dspace_id, dataset_id, error)
  CALL h5dwrite_f(dataset_id, H5T_NATIVE_DOUBLE, arr, dims, error)
  CALL h5sclose_f(dspace_id, error)
  CALL h5dclose_f(dataset_id, error)
ENDSUBROUTINE

SUBROUTINE readarray_1d_double(arr, nme, id)
  use hdf5
  implicit none
  double precision, intent(inout) :: arr(:)
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) dataset_id
  INTEGER(HSIZE_T) :: dims(1)
  integer error

  dims(1) = size(arr)
  CALL h5dopen_f(id, nme, dataset_id, error)
  CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, arr, dims, error)
  CALL h5dclose_f(dataset_id, error)
ENDSUBROUTINE
SUBROUTINE readarray_2d_double(arr, nme, id)
  use hdf5
  implicit none
  double precision, intent(inout) :: arr(:,:)
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) dataset_id
  INTEGER(HSIZE_T) :: dims(2)
  integer error

  dims = shape(arr)
  CALL h5dopen_f(id, nme, dataset_id, error)
  CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, arr, dims, error)
  CALL h5dclose_f(dataset_id, error)
ENDSUBROUTINE
SUBROUTINE readarray_3d_double(arr, nme, id)
  use hdf5
  implicit none
  double precision, intent(inout) :: arr(:,:,:)
  character(*), intent(in) :: nme
  integer(HID_T), intent(in) :: id
  integer(HID_T) dataset_id
  INTEGER(HSIZE_T) :: dims(3)
  integer error

  dims = shape(arr)
  CALL h5dopen_f(id, nme, dataset_id, error)
  CALL h5dread_f(dataset_id, H5T_NATIVE_DOUBLE, arr, dims, error)
  CALL h5dclose_f(dataset_id, error)
ENDSUBROUTINE

MODULE hdf5readwrite
  use hdf5
  implicit none
  INTERFACE writesattr
    SUBROUTINE writesattr_int(val,nme,id)
      use hdf5
      integer, intent(in) :: val
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
    SUBROUTINE writesattr_logical(val,nme,id)
      use hdf5
      logical, intent(in) :: val
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
    SUBROUTINE writesattr_double(val,nme,id)
      use hdf5
      double precision, intent(in) :: val
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
  ENDINTERFACE

  INTERFACE readsattr
    SUBROUTINE readsattr_int(val,nme,id)
      use hdf5
      integer, intent(inout) :: val
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
    SUBROUTINE readsattr_logical(val,nme,id)
      use hdf5
      logical, intent(inout) :: val
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
    SUBROUTINE readsattr_double(val,nme,id)
      use hdf5
      double precision, intent(inout) :: val
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
  ENDINTERFACE  

  INTERFACE writearray
    SUBROUTINE writearray_1d_double(arr, nme, id)
      use hdf5
      double precision, intent(in) :: arr(:)
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
    SUBROUTINE writearray_2d_double(arr, nme, id)
      use hdf5
      double precision, intent(in) :: arr(:,:)
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
    SUBROUTINE writearray_3d_double(arr, nme, id)
      use hdf5
      double precision, intent(in) :: arr(:,:,:)
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
  ENDINTERFACE
  INTERFACE readarray
    SUBROUTINE readarray_1d_double(arr, nme, id)
      use hdf5
      double precision, intent(inout) :: arr(:)
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
    SUBROUTINE readarray_2d_double(arr, nme, id)
      use hdf5
      double precision, intent(inout) :: arr(:,:)
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
    SUBROUTINE readarray_3d_double(arr, nme, id)
      use hdf5
      double precision, intent(inout) :: arr(:,:,:)
      character(*), intent(in) :: nme
      integer(HID_T), intent(in) :: id
    ENDSUBROUTINE
  ENDINTERFACE
ENDMODULE
