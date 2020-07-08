MODULE parmpi
! This module manages the MPI interface for parallel computing
! The parallel scheme decomposes the full E;Mu;z grid into logically cartesian 
! subdomains owned by the local processors
!
! Module Variables:
!  nproc:	Total number of processors
!  myid: 	Local processor ID in MPI_COMM_WORLD
!  cart_comm:	A cartesian domain decomposition communicator
!  mycoords:	Local processor position in cartesian domain
!  neighbors:	The IDs of the local processors's cartesian neighbors
!  dimsizes:	The number of processors along a grid dimension
!  line_comm:   An array of communicators collecting the processors along cartesian lines in each grid dimension  
!  plane_comm:  An array of communicators collecting the processors in cartesian planes in each grid dimension
!  mylineid:	This processors ID in the line_comm communicators
!  myplaneid:	This processors ID in the plane_comm communicators
!  lineids:	The IDs of all processors in this processors line_comm communicators
!  planeids:    The IDs of all processors in this processors plane_comm communicators
!  nall:        The dimensions of the global grid (nE_all,nmu_all, nz_all)
!  nsub:	The dimensions of this processors local grid (nE,nmu, nz)
!  ncells:	The dimensions of local grids collected across all processors
!  fulloffset:	The offset in the global grids where this processors local grid begins
!  offset:      The fulloffset for each processor in the cartesian communicator
!  par_leftedge True if this processor owns the i = 1 (i.e., the left) side of the global grid 
!  par_rightedge True if this processor owns the i = n (i.e., the right) side of the global grid
  use mpi
  implicit none
  integer, parameter :: ndim = 3
  integer, parameter :: minsize = 5
  integer mycoords(ndim), neighbors(-1:1,-1:1,-1:1), myid, nproc, dimsizes(ndim), mylineid(ndim), myplaneid(ndim), fulloffset(ndim)! fulloffset is the offset in the full grid owned by this process. 
  logical par_leftedge(ndim), par_rightedge(ndim) ! par_leftedge (par_rightedge) is true when this process manages the i=1 (i=n) edge of the full grid
  integer, allocatable :: ncells(:,:), offset(:,:), lineids(:,:), planeids(:,:)
  integer nall(ndim), nsub(ndim)
  integer line_comm(ndim), cart_comm, plane_comm(ndim)
  integer mysubarr_type

  contains
  SUBROUTINE init_par()
    implicit none
    integer ierr

    call MPI_Comm_rank(MPI_COMM_WORLD,myid,ierr)
    call MPI_Comm_size(MPI_COMM_WORLD,nproc,ierr)
  ENDSUBROUTINE

  SUBROUTINE setup_par_grids()
    use grid
    implicit none
    integer ierr, nz1, nmu1, nE1, reorder, periods(ndim), npE, npz, npmu, i,j,k, coord(ndim), idx

    !NEED TO ADD CHECKS FOR WHEN NPROC DOES NOT EVENLY DIVIDE GRID SIZES 
    nz1 = min(max(nz_all / nproc, minsize),nz_all) ; npz = min(nz_all / nz1, nproc)
    nmu1 = min(max((nmu_all * npz) / nproc, minsize),nmu_all); npmu = nmu_all / nmu1
    nE1 = min(max( (nE_all * npz*npmu)/nproc, minsize),nE_all); npE = nE_all / nE1
    
    dimsizes = (/npE,npmu,npz/)
    nall = (/nE_all,nmu_all,nz_all/)
    periods = (/0,0,0/)
    reorder = 1

    call MPI_Barrier(MPI_COMM_WORLD,ierr)
    call MPI_Cart_create(MPI_COMM_WORLD, ndim, dimsizes, periods, reorder, cart_comm,ierr)
    call MPI_Comm_rank(cart_comm,myid,ierr)
    call MPI_Cart_coords(cart_comm,myid,ndim, mycoords, ierr)

    nE = nE_all / npE; if (mod(nE_all,npE) .gt. mycoords(1)) nE = nE + 1
    nmu = nmu_all / npmu; if (mod(nmu_all,npmu) .gt. mycoords(2)) nmu = nmu + 1
    nz = nz_all / npz; if (mod(nz_all,npz) .gt. mycoords(3)) nz = nz + 1 
    nsub = (/nE,nmu,nz/)

    do i = 1, ndim
      fulloffset(i) = (nall(i) / dimsizes(i))*mycoords(i) + min(mycoords(i),mod(nall(i),dimsizes(i)))
    enddo

    allocate(ncells(ndim,nproc))
    call MPI_Allgather(nsub, 3, MPI_INTEGER, ncells, 3, MPI_INTEGER, cart_comm, ierr) ! SEND size of array to all processes
    allocate(offset(ndim,nproc))
    call MPI_Allgather(fulloffset,3,MPI_INTEGER, offset, 3, MPI_INTEGER, cart_comm,ierr)

    neighbors = -1
    do k=-1,1; do j=-1,1; do i=-1,1
      if (mycoords(1)+i.lt.0.or.mycoords(1)+i.gt.dimsizes(1)-1 .or. &
          mycoords(2)+j.lt.0.or.mycoords(2)+j.gt.dimsizes(2)-1 .or. &
          mycoords(3)+k.lt.0.or.mycoords(3)+k.gt.dimsizes(3)-1 ) cycle
      call MPI_Cart_rank(cart_comm, mycoords + (/i,j,k/), neighbors(i,j,k),ierr)
    enddo;enddo; enddo
    par_leftedge = .false.
    par_rightedge = .false.
    do i = 1, ndim
      if (mycoords(i).eq.0) par_leftedge(i) = .true.
      if (mycoords(i).ge.dimsizes(i)-1) par_rightedge(i) = .true.
    enddo
    allocate(lineids(maxval(dimsizes),ndim), planeids(maxval(dimsizes),ndim))
    do i = 1, ndim
      coord = mycoords
      coord(i) = 0
      idx = coord(3)*dimsizes(2)*dimsizes(1) + coord(2)*dimsizes(1) + coord(1)
      call MPI_Comm_split(cart_comm, idx, mycoords(i), line_comm(i),ierr) 
      call MPI_Comm_rank(line_comm(i),mylineid(i), ierr)
      call MPI_AllGather(myid,1,MPI_INTEGER, lineids(:,i), 1, MPI_INTEGER, line_comm(i), ierr)
    enddo
    do i = 1, ndim
      call MPI_Comm_split(cart_comm, mycoords(i), mycoords(i), plane_comm(i),ierr)
      call MPI_Comm_rank(plane_comm(i),myplaneid(i), ierr)
      call MPI_AllGather(myid,1,MPI_INTEGER, planeids(:,i), 1, MPI_INTEGER, plane_comm(i), ierr)
    enddo

    call MPI_Type_create_subarray(ndim,nall,nsub,fulloffset,MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, mysubarr_type,ierr)
    call MPI_Type_commit(mysubarr_type,ierr)
  ENDSUBROUTINE

  SUBROUTINE par_updateGC(array)
! Swap guard cell values with neighboring processors
    use grid
    implicit none 
    double precision array(-1:nE+2,-1:nmu+2,-1:nz+2)
    integer i,j,k,nsize,ii,jj,kk,ii1,jj1,kk1,rii,rjj,rkk,rii1,rjj1,rkk1,tag, status(MPI_STATUS_SIZE),ierr

    do k = -1,1; do j = -1,1; do i = -1,1
      if (k.eq.0.and.j.eq.0.and.i.eq.0) cycle
      ii = (i + 1)/2 * (nE-2) + 1; jj = (j+1)/2*(nmu-2) + 1; kk = (k+1)/2*(nz-2)+1
      ii1 = (i +2)/2 * (nE-2) + 2; jj1 = (j+2)/2*(nmu-2) + 2; kk1 = (k+2)/2*(nz-2)+2
      rii = (-i + 1)/2 * (nE-2)-2*i+1; rjj = (-j+1)/2*(nmu-2)-2*j+1; rkk = (-k+1)/2*(nz-2)-2*k+1
      rii1 = (-i +2)/2 * (nE-2)-2*i+2; rjj1 = (-j+2)/2*(nmu-2)-2*j+2; rkk1 = (-k+2)/2*(nz-2)-2*k+2
      nsize = (ii1-ii+1)*(jj1-jj+1)*(kk1-kk+1)

      tag = (k+1)*9 + (j+1)*3 + i+1 +1 
      if (neighbors(i,j,k) .ge.0) then
        call MPI_Send(array(ii:ii1,jj:jj1,kk:kk1),nsize, MPI_DOUBLE_PRECISION, neighbors(i,j,k), tag, cart_comm,ierr)
      endif
      if (neighbors(-i,-j,-k) .ge.0) then
        call MPI_Recv(array(rii:rii1,rjj:rjj1,rkk:rkk1), nsize, MPI_DOUBLE_PRECISION, neighbors(-i,-j,-k),tag, cart_comm, & 
                      status,ierr)
      endif
    enddo; enddo; enddo
    call MPI_Barrier(cart_comm,ierr)
  ENDSUBROUTINE

  SUBROUTINE par_collectall(array, array_all)
    use grid
    implicit none
    double precision array(nE,nmu,nz), array_all(:,:,:)
    integer i,ierr, status(MPI_STATUS_SIZE), recv_type
    if (myid.eq.0) then
      do i = 2, nproc
        call MPI_Type_create_subarray(ndim, nall, ncells(:,i), offset(:,i), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                      recv_type, ierr)
        call MPI_Type_commit(recv_type,ierr)
        call MPI_Recv(array_all,1,recv_type,i-1,i-1,cart_comm,status,ierr)
        call MPI_Type_free(recv_type,ierr)
      enddo
      array_all(1:nE,1:nmu,1:nz) = array
    else
      call MPI_Send(array,nE*nmu*nz, MPI_DOUBLE_PRECISION, 0, myid, cart_comm,ierr)
    endif
  ENDSUBROUTINE

  SUBROUTINE par_scatterall(array, array_all)
    use grid
    implicit none
    double precision array(nE,nmu,nz), array_all(:,:,:)
    integer i,ierr, status(MPI_STATUS_SIZE), send_type
    if (myid.eq.0) then
      do i = 2, nproc
        call MPI_Type_create_subarray(ndim, nall, ncells(:,i), offset(:,i), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                      send_type, ierr)
        call MPI_Type_commit(send_type,ierr)
        call MPI_Send(array_all,1,send_type,i-1,i-1,cart_comm,ierr)
        call MPI_Type_free(send_type,ierr)
      enddo
      array = array_all(1:nE,1:nmu,1:nz)
    else
      call MPI_Recv(array,nE*nmu*nz, MPI_DOUBLE_PRECISION, 0, myid, cart_comm,status, ierr)
    endif
  ENDSUBROUTINE

  SUBROUTINE par_GetGridAll(Eall,muall,zall)
    use grid
    implicit none
    double precision Eall(:), muall(:), zall(:)
    integer n1(dimsizes(1)), n2(dimsizes(2)), n3(dimsizes(3)), disp1(dimsizes(1)), disp2(dimsizes(2)), disp3(dimsizes(3))
    integer ie, ierr
 
    if (mycoords(2) .eq.0.and.mycoords(3).eq.0) then
       ie = nE
       if (par_rightedge(1)) ie = ie+1
       n1 = ncells(1,lineids(1:dimsizes(1),1)+1)
       n1(dimsizes(1)) = n1(dimsizes(1))+1
       disp1 = offset(1,lineids(1:dimsizes(1),1)+1)
       call MPI_Gatherv(E(1:ie),ie,MPI_DOUBLE_PRECISION, Eall, n1, disp1, MPI_DOUBLE_PRECISION, 0, line_comm(1), ierr)
    endif
    if (mycoords(1) .eq.0.and.mycoords(3).eq.0) then
       ie = nmu
       if (par_rightedge(2)) ie = ie+1
       n2 = ncells(2,lineids(1:dimsizes(2),1)+1)
       n2(dimsizes(2)) = n2(dimsizes(2))+1
       disp2 = offset(2,lineids(1:dimsizes(2),2)+1)
       call MPI_Gatherv(mu(1:ie),ie,MPI_DOUBLE_PRECISION, muall, n2, disp2, MPI_DOUBLE_PRECISION, 0, line_comm(2), ierr)
    endif
    if (mycoords(1) .eq.0.and.mycoords(2).eq.0) then
       ie = nz
       if (par_rightedge(3)) ie = ie+1
       n3 = ncells(3,lineids(1:dimsizes(3),3)+1)
       n3(dimsizes(3)) = n3(dimsizes(3))+1
       disp3 = offset(3,lineids(1:dimsizes(3),3)+1)
       call MPI_Gatherv(z(1:ie),ie,MPI_DOUBLE_PRECISION, zall, n3, disp3, MPI_DOUBLE_PRECISION, 0, line_comm(3), ierr)
    endif
    call MPI_Barrier(cart_comm,ierr)
  ENDSUBROUTINE
  
  SUBROUTINE par_Gather1D(array, array_all, n, dir)
    implicit none
    double precision array(:), array_all(:)
    integer n, dir, ierr
    integer, allocatable :: n1(:), disp(:)

    allocate(n1(dimsizes(dir)), disp(dimsizes(dir)))
    disp = offset(dir,lineids(1:dimsizes(dir),dir)+1) 
    n1 = ncells(dir, lineids(1:dimsizes(dir),dir)+1)
    call MPI_Allgatherv(array, n, MPI_DOUBLE_PRECISION, array_all, n1, disp, MPI_DOUBLE_PRECISION, line_comm(dir), ierr)
    deallocate(n1,disp)
  ENDSUBROUTINE

  SUBROUTINE par_CollectMS(array, array_all)
    use grid
    implicit none
    double precision array(:,:), array_all(:,:)
    integer np,i,ierr, status(MPI_STATUS_SIZE), recv_type

    call MPI_Comm_size(plane_comm(3),np,ierr)
    if (myplaneid(3).eq.0) then
      do i = 2, np
        call MPI_Type_create_subarray(2, nall(1:2), ncells(1:2,i), offset(1:2,i), MPI_ORDER_FORTRAN, MPI_DOUBLE_PRECISION, &
                                      recv_type, ierr)
        call MPI_Type_commit(recv_type,ierr)
        call MPI_Recv(array_all,1,recv_type,i-1,i-1,plane_comm(3),status,ierr)
        call MPI_Type_free(recv_type,ierr)
      enddo
      array_all(1:nE,1:nmu) = array(1:nE,1:nmu)
    else
      call MPI_Send(array(1:nE,1:nmu),nE*nmu, MPI_DOUBLE_PRECISION, 0, myplaneid(3), plane_comm(3),ierr)
    endif
  ENDSUBROUTINE

  SUBROUTINE par_CumSum3(cumsum, n1, n2,n3)
    implicit none
    double precision cumsum(n1,n2, n3), cumsum_all(n1,n2,dimsizes(3)), csum(n1,n2)
    integer n1,n2, n3,ierr, i, j,k

    do k = 2, n3; do j = 1, n2; do i= 1, n1
      cumsum(i,j,k) = cumsum(i,j,k) + cumsum(i,j,k-1)
    enddo;enddo;enddo
    if (dimsizes(3).gt.1) then
      call MPI_Allgather(cumsum(:,:,n3), n1*n2, MPI_DOUBLE_PRECISION, cumsum_all, n1*n2,MPI_DOUBLE_PRECISION, line_comm(3),ierr)
      if (mylineid(3).gt.0) then
        csum = sum(cumsum_all(:,:,1:mylineid(3)),3)
        do k = 1, n3
          cumsum(:,:,k) = cumsum(:,:,k) + csum
        enddo
      endif
    endif
    call MPI_BARRIER(line_comm(3),ierr)
  ENDSUBROUTINE
  
  SUBROUTINE par_LineSpaceSum(lsum,nl, dir)
    use grid
    implicit none
    double precision lsum(nl),allsum(nl)
    integer nl, ierr, dir
    
    call MPI_Allreduce(lsum, allsum, nl, MPI_DOUBLE_PRECISION, MPI_SUM, line_comm(dir), ierr)
  ENDSUBROUTINE

  SUBROUTINE par_MomSpaceIntegrate(array, allsum)
    use grid
    implicit none
    double precision array(-1:nE+2,-1:nmu+2,-1:nz+2), mysum(nz), allsum(nz)
    integer ierr, k
    
    do k = 1, nz
      mysum(k) = sum(array(1:nE,1:nmu,k)*msvol)
    enddo
    call MPI_Allreduce(mysum, allsum, nz, MPI_DOUBLE_PRECISION, MPI_SUM, plane_comm(3),ierr)
  ENDSUBROUTINE 
  
  SUBROUTINE par_MomSpaceSum(lsum, nl)
    use grid
    implicit none
    double precision lsum(nl),allsum(nl)
    integer nl, ierr

    call MPI_Allreduce(lsum, allsum, nl, MPI_DOUBLE_PRECISION, MPI_SUM, plane_comm(3), ierr)
    lsum = allsum
  ENDSUBROUTINE 

  SUBROUTINE par_end()
    implicit none
    integer ierr,i
    deallocate (ncells,offset,lineids,planeids)
    do i = 1, ndim
      call MPI_Comm_free(line_comm(i),ierr)
      call MPI_Comm_free(plane_comm(i),ierr)
    enddo
    call MPI_Type_free(mysubarr_type,ierr)
    call MPI_Comm_free(cart_comm,ierr)
  ENDSUBROUTINE
ENDMODULE
