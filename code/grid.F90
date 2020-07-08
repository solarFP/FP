module grid
! This module encapsulates the phase space grid. 
! It is 3D (Energy, pitch-angle, z position)
! Module Variables:
!  nE_all:	Number of energy grid cells in global domain
!  nmu_all:	Number of pitch angle grid cells in global domain
!  nz_all:	Number of position grid cells in global domain
!  nE:		Number of energy grid cells on local processor
!  nmu:		Number of pitch angle grid cells on local processor
!  nz:		Number of position grid cells on local processor
!
! The following grid variables are defined on the local processor
!  E:		Energy grid cell boundaries (normalized in units of mbeam)
!  Em:		Energy grid cell centers (same units as E)
!  gma:		Lorentz relativistic gamma at grid cell boundaries. gma = E + 1
!  gmam:	Lorentz relativistic gamma at grid cell centers. gmam = Em + 1
!  bta:		speed at grid cell boundaries (normalized in units of clight)
!  btam:	speed at cell centers (same units as bta)
!  p: 		momentum at cell boundaries (normalized in units of m clight so p = gma * bta)
!  pm: 		momentum at cell centers (same units as p)
!  dpm:		1st derivative, d/dp, left stencil (uses i-2,i-1,i)
!  dpp:		1st derivative, d/dp, right stencil (uses i,i+1,i+2)
!  dpc:         1st derivative, d/dp, centered stencil (uses i-1,i,i+1)
!  dp2:		2nd derivative, d^2/dp^2, centered stencil (uses i-1,i,i+1)
!
!  theta:	pitch angle grid cell boundaries (radians)
!  thetam:	pitch angle grid cell centers (radians)
!  mu:		cos(theta)
!  mum:		cos(thetam)
!  sint:	sin(theta)
!  sintm:	sin(thetam)
!  dmum:	1st derivative, d/dmu, left stencil (j-2,j-1,j)
!  dmup:	1st derivative, d/dmu, right stencil (j,j+1,j+2)
!  dmuc:	1st derivative, d/dmu, centered stencil (j-1,j,j+1)
!  dmu2:	2nd derivative, d^2/dmu^2, centered stencil (j-1,j,j+1)
!
!  z:		position from loop top grid cell boundaries (cm)
!  zm:		position from loop top grid cell centeres (cm)
!  dzm:		1st derivative, d/dz, left stencil (k-2,k-1,k)
!  dzp:		1st derivative, d/dz, right stencil (k,k+1,k+2)
!  dzc:		1st derivative, d/dz, centered stencil (k-1,k,k+1)
!  dz2:		2nd derivative, d^2/dz^2, centered stencil (k-1,k,k+1)
!
!  msvol:	momemtum space differential volume element, d3p. (units mbeam eV / clight)
!
  use const
  implicit none

  integer nE, nmu, nz ! The sizes of the sub-portion of the full grid owned by this process
  integer nE_all,nmu_all, nz_all ! The sizes of the full grid
  logical doRel

  double precision, allocatable:: E(:),p(:), Em(:), pm(:), gma(:), bta(:), gmam(:), btam(:), &
                                  dpm(:,:), dpp(:,:), dp2(:,:), dpc(:,:)  !Energy related variables
  double precision, allocatable:: theta(:), mu(:), thetam(:), mum(:), sint(:), sintm(:), dmum(:,:), dmup(:,:) 
  double precision, allocatable:: dmu2(:,:), dmuc(:,:), dtc(:,:) ! pitch-angle related variables
  double precision, allocatable:: z(:), zm(:), dzm(:,:), dzp(:,:), dz2(:,:), dzc(:,:), msvol(:,:)
  
  contains 
    subroutine allocate_grid(oneD)
      implicit none
      logical oneD

      if (nE.lt.1) nE = 200 ! Default nE
      if (nmu.lt.1) nmu = 60 ! default nmu
      if (nz.lt.1) nz = 5

      if ((.not.oneD).and.(mod(nmu_all,2) .ne.0)) nmu_all = nmu_all+1 ! ensure nmu_all is even
      allocate(E(-1:nE+3), p(-1:nE+3), Em(-1:nE+2), pm(-1:nE+2)) ! Energy and momentum grids
      allocate(gma(-1:nE+3), bta(-1:nE+3), gmam(-1:nE+2), btam(-1:nE+2)) ! gamma and beta 
      allocate(dpm(-2:2,nE), dpp(-2:2,nE), dp2(-2:2,nE), dpc(-1:1,nE)) ! 3-point 1st and 2nd derivative stencils
   
      allocate(theta(-1:nmu+3), mu(-1:nmu+3), thetam(-1:nmu+2), mum(-1:nmu+2))! theta and mu = cos(theta) grids
      allocate(sint(-1:nmu+3), sintm(-1:nmu+2)); ! sin(theta) 
      allocate(dmum(-2:2,nmu), dmup(-2:2,nmu), dmu2(-2:2,nmu), dtc(-1:1,nmu),dmuc(-1:1,nmu))! 3-point 1st and 2nd derivative stencils
      allocate(z(-1:nz+3), zm(-1:nz+2), dzm(-2:0,nz), dzp(0:2,nz), dz2(-1:1,nz),dzc(-1:1,nz))
      allocate(msvol(nE,nmu))
    end subroutine 

    subroutine deallocate_grid()
      implicit none
      deallocate(E, p, Em, pm)
      deallocate(gma, bta, gmam, btam)
      deallocate(dpm, dpp, dp2, dpc)
    
      deallocate(theta,mu,thetam,mum)
      deallocate(sint,sintm)
      deallocate(dmum,dmup,dmu2,dmuc, dtc)

      deallocate(z, zm,dzm,dzp,dz2,dzc, msvol)
    end subroutine

    subroutine init_grid(mbeam, Emax, Emin, Epiv, rel, zin_all,fulloffset,oneD)
      use derivs
      implicit none
      integer i,j,k,ib, jall, fulloffset(3), nEbp
      double precision Emax, Emx, mbeam, Emin, Emn, Epiv, Epv,z0, zma(-1:nz_all+2), za(-1:nz_all+3), zin_all(nz_all), zn,gs
      logical rel, oneD

      doRel = rel
      call allocate_grid(oneD)
      
      Emx = Emax / mbeam
      Emn = Emin / mbeam
      Epv = Epiv / mbeam

      gs=.5d0
      ! E grid in units of /mc^2
      nEbp = nE_all / 2 ! # of grid cells below Epiv
      
      ib = -1
      do i = ib, nE+3
!        E(i) = (real(i-1 + fulloffset(1))/real(nE_all+1) * (Emx**gs - Emn**gs) + Emn**gs)**(1d0/gs)
        E(i) = exp(real(i-1 + fulloffset(1))/real(nE_all+1) * (log(Emx) - log(Emn)) + log(Emn))
!        if (i+fulloffset(1).lt.nEbp) then
!          E(i) = (real(i-1 + fulloffset(1))/real(nEbp+1) * (Epv**gs - Emn**gs) + Emn**gs)**(1d0/gs) !power-law (gs) below Epiv
!        else
!          E(i) = exp(real(i-1-nEbp + fulloffset(1))/real(nE_all-nEbp+1) * (log(Emx) - log(Epv)) + log(Epv)) !log spaced above Epiv
!        endif
      enddo
      Em = .5*(E(-1:nE+2) + E(0:nE+3))
      if (doRel) then
        gma = E + 1
        gmam = Em + 1
        bta = sqrt(1 - 1/gma**2)
        btam = sqrt(1 - 1/gmam**2)
        where(E.gt.0.and.E.lt.1d-12) ! catch very small energies to avoid roundoff error
          bta = sqrt(2*E)
        endwhere
        where(Em.gt.0.and.Em.lt.1d-12)
          btam = sqrt(2*Em)
        endwhere
      else 
        gma = 1d0
        gmam = 1d0
        bta = sqrt(2*E)
        btam = sqrt(2*Em)
      endif
      p = gma*bta
      pm = gmam*btam
      call deriv2_3coef(pm,dp2(-1:1,:),nE); dp2(-2,:) = 0; dp2(2,:) = 0
      call deriv_3coef(pm(0:nE+1), pm(1:nE), dpc, nE)
      call deriv_3coef(pm(-1:nE), pm(1:nE), dpm(-2:0,:), nE) ; dpm(1:2,:) = 0! 3pt deriv stencil evalautes at pm(i) using pm(i-2), pm(i-1), pm(i)
      call deriv_3coef(pm(1:nE+2), pm(1:nE), dpp(0:2,:), nE); dpp(-2:-1,:) =0! 3pt deriv stencil evalautes at pm(i) using pm(i), pm(i+1), pm(i+2)
   
      dp2 = dp2 * (clight/mbeam)**2
      dpm = dpm * (clight/mbeam) 
      dpp = dpp * (clight/mbeam)
      dpc = dpc * (clight/mbeam)
      if (.not.oneD) then
        do j = -1, nmu+3
          ! setup mu grid
          jall = j + fulloffset(2)
          if (jall .le. nmu_all / 2 + 1) then 
            theta(j) = sign(((jall-1.)/(nmu_all/2))**1d0,jall-1d0) * pi/2 
          else 
            jall = nmu_all-jall+2
            theta(j) = pi - sign(((jall-1.)/(nmu_all/2))**1d0,jall-1d0) * pi/2
          endif
        enddo
        thetam = .5*(theta(-1:nmu+2) + theta(0:nmu+3))
        mu = cos(theta)
      
        mum = cos((theta(-1:nmu+2) + theta(0:nmu+3))*.5)
        sint = sin(theta)
        sintm = sin((theta(-1:nmu+2) + theta(0:nmu+3))*.5)
        call deriv2_3coef(thetam, dmu2(-1:1,:), nmu) ; dmu2(-2,:) = 0; dmu2(2,:) = 0
        call deriv_3coef(thetam(0:nmu+1),thetam(1:nmu), dmuc(-1:1,:), nmu) ; dtc = dmuc
        call deriv_3coef(thetam(-1:nmu), thetam(1:nmu), dmum(-2:0,:), nmu) ; dmum(1:2,:) = 0
        call deriv_3coef(thetam(1:nmu+2), thetam(1:nmu), dmup(0:2,:), nmu) ; dmup(-2:-1,:) = 0

        do j = 1, nmu
          dmum(:,j) = -dmum(:,j) / sintm(j)
          dmup(:,j) = -dmup(:,j) / sintm(j)
          dmuc(:,j) = -dmuc(:,j) / sintm(j)
          dmu2(-1:1,j) = (dmu2(-1:1,j) + dmuc(:,j)*mum(j))/sintm(j)**2
        enddo
      else !if oneD
        do j=-1,3
          theta(j) = pi*(j-1)
          thetam(j) = theta(j)
          mu(j) = nint(cos(theta(j)))
          mum(j) = mu(j)
          sint(j) = nint(sin(theta(j)))
          sintm(j) = sintm(j) 
        enddo
        dmum = 0; dmup=0; dmuc=0; dmu2 = 0
        dmum(0,:) = 1; dmup(0,:) = 1; dmuc(0,:) = 1; dmu2(0,:) = 1
      endif
      ! z grid
      z0 = zin_all(1) + (zin_all(1)-zin_all(2))/2. ! top boundary is 1/2 grid cell past the top cell center 
      zma(1:nz_all) = z0 - zin_all
      zma(nz_all+1) = zma(nz_all) + (zma(nz_all)-zma(nz_all-1))
      zma(nz_all+2) = zma(nz_all+1) + (zma(nz_all+1)-zma(nz_all))
      za(1) = 0
      do k = 1, nz_all+1
        zn = 2*zma(k) - za(k)
        if (zn.gt.zma(k+1)) zn = .5*(zma(k)+zma(k+1))
        za(k+1) = zn
      enddo
      za(nz_all+3) = 2*zma(nz_all+2) - za(nz_all+2)
      za(0) = -za(2)
      zma(0) = .5*(za(0)+za(1))
      za(-1) = -za(3)
      zma(-1) = .5*(za(-1)+za(0))
      z = za(fulloffset(3)-1:fulloffset(3)+nz+3)
      zm = zma(fulloffset(3)-1:fulloffset(3)+nz+2)
      call deriv_3coef(zm(-1:nz), zm(1:nz), dzm, nz)
      call deriv_3coef(zm(1:nz+2), zm(1:nz), dzp, nz)

      call deriv2_3coef(zm,dz2, nz)
      call deriv_3coef(zm(0:nz+1),zm(1:nz),dzc, nz)
      do k = 1, nz
        dzm(:,k) = 0; dzm(0,k) = 1d0/(zm(k)-zm(k-1)); dzm(-1,k) = -dzm(0,k)
        dzp(:,k) = 0; dzp(0,k) = -1d0/(zm(k+1)-zm(k)); dzp(1,k) = -dzp(0,k)
      enddo
      do j = 1, nmu; do i = 1, nE
        msvol(i,j) = (mbeam / clight)**3 * (p(i+1)**3 - p(i)**3)/3.*(mu(j)-mu(j+1)) *2*pi
      enddo; enddo
    end subroutine
end module
