MODULE beam
! This module encapsulates the distribution of injected (beamed) particles.
! Module Variables:
!   f:	 	The current iteration estimate of the distribution function (particles / cm^3 / d3p) 
!   fo:	 	The previous iteration's value of f (particles / cm^3/ d3p).
!   injectedf:  The distribution function injected at the z=0 boundary (particles /cm^3 / d3p). 
!   nflux:	The number flux as a function of z (particles / cm^2 /s.
!   mbeam:	Mass of the beam particle (eV).
!   zbeam:	Charge of the beam particle in (elemenentary charge units)
  implicit none
  double precision, allocatable :: f(:,:,:), fo(:,:,:)
  double precision, allocatable:: injectedf(:,:), nflux(:)
  double precision mbeam, zbeam
  
  contains
  SUBROUTINE init_beam()
    use grid
    implicit none

    allocate(f(-1:nE+2,-1:nmu+2, -1:nz+2), fo(nE,nmu, nz), injectedf(nE,nmu), nflux(nz))
    f =0 ; fo = 0
  ENDSUBROUTINE

  SUBROUTINE deallocate_beam()
    implicit none
    deallocate(f,fo,injectedf,nflux)
  ENDSUBROUTINE

  SUBROUTINE InjectPowerLaw(Ec, dlt, eflx,patype, pasigma)
    ! Inputs:
    !  Ec is cutoff energy is keV
    !  dlt is power law index
    !  eflx is injected energy flux erg/cm^2/s
    !  patype determines the pitch angle distribution of the beam
    !    patype = 0 fully beamed (all injected particles have mu =1)
    !    patype = 1 isotropic in foward hemisphere
    !    patype = 2 Gaussian in forward hemisphere
    !  pasigma is the sigma of the pitch angle distribution for a Gaussian
    use grid
    use parmpi
    implicit none
    double precision Ec, dlt, eflx, pasigma, nrm, Ecut, eflux, padist(nmu), injectedflux, nn(1)
    integer patype, i,j

    padist = 0
    if (patype.eq.0) then
      if (par_leftedge(2)) padist(1) = 1
    else if (patype.eq.1) then
      padist = 1
    else if (patype.eq.2) then
      padist = exp(-.5*(thetam(1:nmu)**2/pasigma**2))/pasigma/sqrt(2*pi)
    endif
    where (mum(1:nmu) .le.0)
      padist = 0
    endwhere
    nrm = sum(padist*(mu(1:nmu)-mu(2:nmu+1)))
    nn(1) = nrm
    call par_LineSpaceSum(nn,1,2)
    padist = padist / nn(1) /( 2*pi)

    Ecut = Ec * 1d3 ! in eV
    eflux = eflx/ergperev ! in eV
    nrm = eflux / Ecut**2*(dlt-2)
    do j=1, nmu; do i=1, nE
      if (Em(i)*mbeam .ge. Ecut) then
        injectedflux = nrm * (Em(i)*mbeam/Ecut)**(-dlt)*padist(j)
      else
!        injectedflux = nrm*(Em(i)*mbeam/Ecut)**(30)*padist(j)
        injectedflux=0
      endif
      injectedf(i,j) = injectedflux / (clight*btam(i)) * (E(i+1)-E(i))*mbeam*2*pi*(mu(j)-mu(j+1))/msvol(i,j)
    enddo; enddo
    nrm = 0
    do j=1, nmu; do i=1, nE
      nrm = nrm + injectedf(i,j)*mum(j)*btam(i)*clight*Em(i)*mbeam * msvol(i,j)    ! normalize so that z component of energy flux = eflux
    enddo; enddo
    nn = nrm
    call par_MomSpaceSum(nn,1)
    nrm = eflux / nn(1)
    injectedf = injectedf * nrm
  ENDSUBROUTINE
  
  SUBROUTINE CalcNflux()
    use grid
    use parmpi
    implicit none
    integer i,j,k
    nflux = 0
    do k = 1, nz; do j = 1, nmu; do i = 1, nE
      nflux(k) = nflux(k) + f(i,j,k)*mum(j)*btam(i)*clight*msvol(i,j)
    enddo;enddo; enddo
    call par_MomSpaceSum(nflux,nz)
  ENDSUBROUTINE

  FUNCTION CalcDiff() result(diff)
    use helpers
    use options
    use grid
    use parmpi
    implicit none
    double precision diff,d(2)
    integer ierr
  
    d = (/sum((f(1:nE,1:nmu,1:nz)-fo)**2), sum((f(1:nE,1:nmu,1:nz)+fo)**2)/4 /)
    call MPI_Allreduce(MPI_IN_PLACE, d, 2, MPI_DOUBLE_PRECISION, MPI_SUM,  cart_comm, ierr)
    diff = sqrt(safedivide(d(1),d(2),0d0)) / implicit_theta
    fo = f(1:nE,1:nmu,1:nz)
  ENDFUNCTION

  SUBROUTINE UpdateBC()
    use options
    use helpers
    use grid
    use loop
    use parmpi
    implicit none
    integer i,j

    !REFLECT AT BOUNDARIES
    !NEED TO MAKE THIS WORK FOR WHEN MU IS SPLIT ACROSS PROCESSORS
    if (reflecttop.and.par_leftedge(3)) then
      do j = 1, nmu; do i = 1, nE
        f(i,nmu+1-j,0) = f(i,j,1)*safedivide(msvol(i,j),msvol(i,nmu+1-j),0d0) + injectedf(i,nmu+1-j)
        f(i,nmu+1-j,-1) = f(i,j,2)*safedivide(msvol(i,j),msvol(i,nmu+1-j),0d0) + injectedf(i,nmu+1-j)
      enddo;enddo
    else if (par_leftedge(3)) then
      f(1:nE,1:nmu,0) = injectedf
      f(1:nE,1:nmu,-1) = injectedf
    endif
    if (reflectbottom.and.par_rightedge(3)) then
      do j = 1, nmu; do i = 1, nE
        f(i,nmu+1-j,nz+1) = f(i,j,nz)*safedivide(msvol(i,j),msvol(i,nmu+1-j),0d0)
        f(i,nmu+1-j,nz+2) = f(i,j,nz-1)*safedivide(msvol(i,j),msvol(i,nmu+1-j),0d0)
      enddo;enddo
    else if (par_rightedge(3)) then
      f(1:nE,1:nmu,nz+1) = f(1:nE,1:nmu,nz)
      f(1:nE,1:nmu,nz+2) = f(1:nE,1:nmu,nz) ! or should this also be set to 0 ?
    endif

    if (par_leftedge(1)) then
      f(0,1:nmu,:) = f(1,1:nmu,:); f(-1,1:nmu,:) = f(1,1:nmu,:)
!      f(0,1:nmu,:) = 0 ; f(-1,1:nmu,:) =0
    endif
    if (par_rightedge(1)) then
      !f(nE+1,1:nmu,:) = f(nE,1:nmu,:); f(nE+2,1:nmu,:) = f(nE,1:nmu,:)
      f(nE+1,1:nmu,:) = 0; f(nE+2,1:nmu,:) = 0
    endif

    if (par_leftedge(2)) then
      f(:,0,:) = f(:,1,:); f(:,-1,:) = f(:,2,:)
    endif
    if (par_rightedge(2)) then
      f(:,nmu+1,:) = f(:,nmu,:); f(:,nmu+2,:) = f(:,nmu-1,:)
    endif

    call par_updateGC(f)
  ENDSUBROUTINE

ENDMODULE
