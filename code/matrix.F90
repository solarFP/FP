Module matrix
! This module contains the matrix representing the Fokker-Planck operator.
! See Eqs 20-22 in Allred et al (2020)
! The variable "A" contains the linear components to the matrix and B contains 
! the nonlinear component prop to nflux 
! So matrix = A + B * nflux
! Module Variables:
! A:		The linear compoments of the matrix equation
! B:		The nonlinear components prop to nflux of the matrix equation
! rhs:		The right-hand side of the matrix equation
! res:		The residual of the current solution (Eq 23 of Allred et al 2020)
! Ap2:		The component of A prop to d^2 /dp^2
! Ap:		The component of A prop to d/dp
! Amu2:         The component of A prop to d^2/dmu^2
! Amu:          The component of A prop to d/dmu
! Az:           The component of A prop to d/dz
! Af:           The component of A prop to no derivatives
! Bp:           The component of B prop to d/dp
! Bmu:          The component of B prop to d/dmu
! Bf:           The component of B prop to no derivatives

  use grid
  implicit none
  
  double precision, allocatable :: A(:,:,:,:), rhs(:,:,:), res(:,:,:), B(:,:,:,:), flxpC(:,:,:,:), flxtC(:,:,:,:)
  double precision,allocatable :: Ap2(:,:), Ap(:,:,:), Amu2(:,:,:), Amu(:,:,:), Az(:,:), Af(:,:,:), Bmu(:,:,:), Bp(:,:,:), Bf(:,:,:)
  contains
  
  SUBROUTINE SetupMatrix()
    !Allocate memory for the matrix coefficients and populate them
    implicit none
    integer i,j,k
 
    allocate(Ap2(nE,nz), Ap(nE,nmu,nz), Amu2(nE, nmu,nz), Amu(nE, nmu, nz), Az(nE,nmu), Af(nE,nmu, nz)) 
    allocate(Bp(nE,nmu,nz), Bmu(nE,nmu,nz), Bf(nE,nmu,nz))
    allocate(A(0:6,0:nE+1,0:nmu+1, 0:nz+1), rhs(0:nE+1,0:nmu+1, 0:nz+1), res(nE,nmu,nz))
    allocate(B(0:4,nE,nmu,nz))
    allocate(flxpC(-1:1,nE,nmu,nz),flxtC(-1:1,nE,nmu,nz))
    call MatrixCoeff()
    do k=1, nz; do j=1, nmu; do i =1,nE
      A(0,i,j,k) = (Ap2(i,k)*dp2(0,i) + Ap(i,j,k)*dpc(0,i) + Amu2(i,j,k)*dmu2(0,j) + Amu(i,j,k)*dmuc(0,j) & 
                   +Az(i,j)*dzc(0,k) + Af(i,j,k))
      A(1,i,j,k) = (Ap2(i,k)*dp2(-1,i) + Ap(i,j,k)*dpc(-1,i))
      A(2,i,j,k) = (Ap2(i,k)*dp2( 1,i) + Ap(i,j,k)*dpc( 1,i))
      A(3,i,j,k) = (Amu2(i,j,k)*dmu2(-1,j) + Amu(i,j,k)*dmuc(-1,j))
      A(4,i,j,k) = (Amu2(i,j,k)*dmu2( 1,j) + Amu(i,j,k)*dmuc( 1,j))
      A(5,i,j,k) = Az(i,j)*dzc(-1,k)
      A(6,i,j,k) = Az(i,j)*dzc( 1,k)

      B(0,i,j,k) = (Bp(i,j,k) * dpc(0,i) + Bmu(i,j,k)*dmuc(0,j)) + Bf(i,j,k)
      B(1,i,j,k) = Bp(i,j,k) * dpc(-1,i)
      B(2,i,j,k) = Bp(i,j,k) * dpc( 1,i)
      B(3,i,j,k) = Bmu(i,j,k) * dmuc(-1,j)
      B(4,i,j,k) = Bmu(i,j,k) * dmuc( 1,j)
    enddo; enddo; enddo
    rhs = 0
    call ImposeBC()
  ENDSUBROUTINE

  SUBROUTINE deallocate_matrix()
    implicit none
    deallocate(Ap2,Ap,Amu2,Amu,Az,Af,Bp,Bmu,Bf,A,B,rhs,res,flxpC,flxtC)
  ENDSUBROUTINE

  SUBROUTINE ImposeBC()
    use beam
    implicit none
    
    A(0,0,:,:) = 1; A(2,0,:,:) = -1; A(1,0,:,:) = 0; A(3:,0,:,:) = 0 ! df/dp = 0 at lower p boundary
!    A(0,nE+1,:,:) = 1; A(1:,nE+1,:,:) = 0 ! at upper p boundary f = 0
!    rhs(nE+1,:,:) = 0
    A(0,nE+1,:,:) = 1; A(1,nE+1,:,:) = -1; A(2:,nE+1,:,:) = 0 ! df/dp = 0 at upper p boundary
    rhs(nE+1,:,:) = 0

    A(0,:,0,:) = 1; A(4,:,0,:) = -1   ! d/dmu = 0 at mu boundaries
    rhs(:,0,:) = 0
    A(0,:,nmu+1,:) = 1; A(3,:,nmu+1,:) = -1
    rhs(:,nmu+1,:) = 0

    A(0,:,:,0) = 1; A(1:,:,:,0) = 0; !impose injected flux at top boundary
    rhs(1:nE,1:nmu,0) = injectedf
    A(0,:,:,nz+1) = 1; A(1:,:,:,nz+1) = 0 ! at bottom of loop f is maxwellian
    rhs(:,:,nz+1) = 0!MaxwellDist(:,nz+1)
  ENDSUBROUTINE
  
  SUBROUTINE MatrixCoeff()
! Calculate the coefficients (Eq 20 of Allred et al 2020)
    use loop
    use beam
    use options
    implicit none
    double precision mr, Thetab, xb, xip, xi1, xi0, cl, kab, clp, kan, kpan, dlnBdz, S, eta,tmp,tmp2
    integer i,j,k,l
    logical doRelColl
    Ap2 = 0; Ap=0; Amu2=0; Amu=0; Az=0; Af=0; flxpC = 0; flxtC = 0
    Bp = 0; Bmu = 0; Bf = 0
    doRelColl = inc_relativity ! Do relativistic or non-relativistic Coulomb collisions?
    !doRelColl = .false. 
    do k = 1, nz; do i = 1, nE
      do l = 1, Nion
        mr = mbeam / mion(l)
        Thetab = (kb*tg(k)) / mion(l)
        call Calcxi(pm(i),Thetab, xb, xi0,xi1,xip,doRelColl)
        cl = CoulogCC(mbeam,Zbeam, mion(l), Zion(l), dni(l,k), xb, btam(i))
        kab = pi4e4*(Zion(l)**Zbeam)**2*cl*mbeam

!        if (inc_CC_el) then
        if (inc_CC.or.tg(k).lt.2d4) then
          if (doRelColl) then
            tmp2 = -dni(l,k)*kab*xi1/2/xb/pm(i)/mbeam/clight * gmam(i)
          else
            tmp2 = -dni(l,k)*kab*xi1/2/xb/pm(i)/mbeam/clight
          endif
          tmp = dni(l,k)*kab/pm(i)**2/mbeam**2
          Ap2(i,k) = Ap2(i,k) + tmp2
          do j = 1, nmu
            if (doRelColl) then
              Ap(i,j,k) = Ap(i,j,k) + tmp*(gmam(i)*xi1/2/xb - xip*Thetab - xi1*Thetab/gmam(i) - mr*xi1)
              Af(i,j,k) = Af(i,j,k) - dni(l,k)*kab*clight*xip/mbeam**2/mion(l)/gmam(i)/pm(i)
            else
              Ap(i,j,k) = Ap(i,j,k) + tmp*(xi1/2/xb - xip -mr*xi1)
              Af(i,j,k) = Af(i,j,k) - mr*dni(l,k)*kab*clight/pm(i)**3/mbeam**3*xip*2*xb
            endif
            flxpC(:,i,j,k) = flxpC(:,i,j,k) + tmp2*dpc(:,i)
            flxpC(0,i,j,k) = flxpC(0,i,j,k) - tmp*xi1*mr
          enddo
        endif
        if (inc_cc.and..not.oneD) then 
          if (doRelColl) then
            tmp = -dni(l,k)*kab*clight/2/pm(i)**3/mbeam**3/gmam(i)*(xi0 + gmam(i)*Thetab*xip - xi1/2/xb)
          else
            tmp = -dni(l,k)*kab*clight/2/pm(i)**3/mbeam**3*(xi1 + xip - xi1/2/xb)
          endif
          Amu2(i,:,k) = Amu2(i,:,k) + tmp*sintm(1:nmu)**2
          Amu(i,:,k) = Amu(i,:,k) -tmp*2*mum(1:nmu)
          do j =1,nmu
            flxtC(:,i,j,k) = flxtC(:,i,j,k) + tmp*dtc(:,j)*pm(i)*mbeam/clight
          enddo
        endif 
      enddo
      do l = 1, Nneutral
        if (abs(mbeam-me).lt.1d3) then ! if mbeam is within 1keV of electron mass assume its an electron
          cl = CoulogEN(btam(i),gmam(i),Em(i),Enion(l))
        else
          cl = CoulogHN(btam(i),gmam(i),Enion(l))
        endif
        clp =  CoulogND(btam(i),gmam(i),zn(l))
        mr = mbeam / me
        kan = pi4e4*Zbeam**2*Zn(l)*cl*mbeam
        kpan = pi4e4*(Zbeam*Zn(l))**2*clp*mbeam
!        if (inc_CC_el) then 
        if (inc_CC.or.tg(k).lt.2d4) then
          if (doRelColl) then
            tmp = - mr*dnn(l,k)*kan/pm(i)**2/mbeam**2*gmam(i)**2
          else
            tmp = - mr*dnn(l,k)*kan/pm(i)**2/mbeam**2
          endif
          Ap(i,:,k) = Ap(i,:,k) + tmp
          flxpC(0,i,:,k) = flxpC(0,i,:,k) + tmp
        endif
        if (inc_CC.and..not.oneD) then 
          if (doRelColl) then
            tmp = -dnn(l,k)*kpan*clight/2/pm(i)**3/mbeam**3/gmam(i)
          else
            tmp = -dnn(l,k)*kpan*clight/2/pm(i)**3/mbeam**3
          endif
          Amu2(i,:,k) = Amu2(i,:,k) + tmp * sintm(1:nmu)**2
          Amu(i,:,k) = Amu(i,:,k) - tmp*(2*mum(1:nmu))
          do j =1, nmu
            flxtC(:,i,j,k) = flxtC(:,i,j,k) + tmp*dtc(:,j)*pm(i)*mbeam/clight
          enddo
        endif
      enddo
    enddo; enddo
    S = 0 ; dlnBdz = 0 !Initialize Synchro and mag mirror scales 
    do k = 1, nz 
      if (inc_synchro) S = 2./3.*(e2*Zbeam**2/mbeam*Bfield(k))**2 / ergperev ! ev / cm ! Include synchrotron terms?
      if (inc_magmirror.and.(.not.oneD)) dlnBdz = sum(log(Bfield(k-1:k+1))*dzc(:,k)) ! Include Magmirror term?   
      do j = 1, nmu; do i = 1, nE
         tmp =  -S*btam(i)*gmam(i)**2*sintm(j)**2 
         Ap(i,j,k) = Ap(i,j,k) + tmp
         if (.not.oneD) then
           Amu(i,j,k) = Amu(i,j,k) + S*clight/gmam(i)**2/mbeam*sintm(j)**2*mum(j) &
                                   -.5/gmam(i)*btam(i)*clight*dlnbdz*sintm(j)**2
         endif
         Af(i,j,k) = Af(i,j,k) - 4*S*clight/gmam(i)/mbeam*(gmam(i)**2*sintm(j)+mum(j)**2-.5) &
                               + btam(i)*clight/gmam(i)*mum(j)*dlnbdz 
         if (oneD) Af(i,j,k) = Af(i,j,k) -S*btam(i)*(1-3*mum(j)**2) 
         !Do not add these forces to flxpC since they do not contribute to the local
         !heating rate. That is synchrotron radiation is thin and escapes so does
         !not heat the local plasma.
         !flxpC(0,i,j,k) = flxpC(0,i,j,k) + tmp
      enddo; enddo
    enddo
    do j = 1, nmu; do i = 1, nE
      Az(i,j) = mum(j)*btam(i)*clight
    enddo; enddo
    if (inc_RC) then
      do k = 1, nz 
        eta = 0
        do l = 2, Nion
        ! add eta for ambient drifting electron and ambient ion collisions
          eta = eta +  resist_fact * CalcEtaIon(mion(l), Zion(l), dni(l,k), dni(1,k), tg(k))
        enddo
        ! add eta for ambient drifting electron and ambient neutral hydrogen collisions
        eta = eta + resist_fact * CalcEtaNeut(dnn(1,k),dni(1,k),tg(k))
        do j = 1, nmu; do i = 1, nE
          if (.not.oneD) Bmu(i,j,k) = -Zbeam**2*e2*eta*sintm(j)**2/pm(i) *clight/mbeam
          Bp(i,j,k) = -Zbeam**2*e2*eta*mum(j)
          if (oneD) Bf(i,j,k) = Bp(i,j,k)*2/pm(i)*clight/mbeam
        enddo; enddo
      enddo
    endif
  ENDSUBROUTINE

  SUBROUTINE CalcHeatRate(qh, momrate)
    use helpers
    use beam
    use loop
    use parmpi
    implicit none
    double precision qh(nz), momrate(nz), flxp, flxt, ffz
    integer i,j,k
    qh = 0
    momrate = 0
    do k = 1, nz; do j = 1, nmu; do i = 1, nE
      flxp = sum(f(i-1:i+1,j,k)*flxpC(:,i,j,k)) + f(i,j,k) * Bp(i,j,k)*nflux(k) 
      flxt = sum(f(i,j-1:j+1,k)*flxtC(:,i,j,k)) - f(i,j,k) * safedivide(Bmu(i,j,k)*pm(i)*mbeam*nflux(k),clight*sintm(j),0d0)
      ffz = flxp*mum(j) - flxt*sintm(j)
      qh(k) = qh(k) - msvol(i,j) * btam(i)*clight*flxp*ergperev
      momrate(k) = momrate(k) - msvol(i,j)*ffz*ergperev
    enddo;enddo;enddo
    if (par_leftedge(1)) then ! add any energy flowing out of the grid boundaries
      do k = 1, nz; do j = 1, nmu
        flxp = min(sum(f(0:2,j,k)*flxpC(:,1,j,k)) + Bp(1,j,k)*nflux(k)*f(1,j,k),0d0)
        qh(k) = qh(k) - flxp*(mu(j)-mu(j+1))*2*pi*pm(1)**2*Em(1)*mbeam**3/clight**2*ergperev
        momrate(k) = momrate(k) -  flxp*(mu(j)-mu(j+1))*2*pi*(pm(1)*mum(j))*pm(1)**2*mbeam**3/clight**3*ergperev
      enddo; enddo
    endif
    if (par_rightedge(1)) then 
      do k = 1, nz; do j = 1, nmu
        flxp = max(sum(f(nE-1:nE+1,j,k)*flxpC(:,nE,j,k)) + Bp(nE,j,k)*nflux(k)*f(nE,j,k),0d0)
        qh(k) = qh(k) + flxp*(mu(j)-mu(j+1))*2*pi*pm(nE)**2*Em(nE)*mbeam**3/clight**2*ergperev
        momrate(k) = momrate(k) + flxp*(mu(j)-mu(j+1))*2*pi*(pm(nE)*mum(j))*pm(nE)**2*mbeam**3/clight**3*ergperev
      enddo; enddo
    endif

    call par_MomSpaceSum(qh,nz)
    call par_MomSpaceSum(momrate,nz)
  ENDSUBROUTINE

  SUBROUTINE SweepP()
    use helpers
    use options
    use parmpi
    use beam
    implicit none
    integer i,j,k
    double precision rhsd(0:nE+1), diag(0:nE+1,-2:2), muterm(-2:2), zterm(-2:2), forcep, forcet
    double precision newf(0:nE+1,nmu,nz), tt
 
    tt =implicit_theta
    rhsd(0) = 0
    diag(0,:) = 0
    diag(0,0) = 1
    diag(0,1) = -1

    rhsd(nE+1) = 0
    diag(nE+1,:) = 0d0
    diag(nE+1,0) = 1d0
    diag(nE+1,-1) = -1d0

    do k = 1, nz; do j = 1, nmu
      do i = 1,nE
        forcep = (Ap(i,j,k)+Bp(i,j,k)*nflux(k))*dpm(0,i)! Current df/dp component of the force in the p-hat direction
        forcet = (Amu(i,j,k)+Bmu(i,j,k)*nflux(k))*dmum(0,j)! Current df/dmu component of the force in the theta-hat direction

        if (forcep .gt.0) then ! acceleration so upwind derivative stencil is dpm
          diag(i,:) = (Ap(i,j,k) + Bp(i,j,k)*nflux(k)) * dpm(:,i) + Ap2(i,k)*dp2(:,i) 
          if (i.eq.1.and.par_leftedge(1)) diag(i,:) = 0
        else ! deceleration so upwind derivative stencil is dpp
          diag(i,:) = (Ap(i,j,k) + Bp(i,j,k)*nflux(k)) * dpp(:,i) + Ap2(i,k)*dp2(:,i) 
          if (i.eq.nE.and.par_rightedge(1)) diag(i,:) = 0
        endif
        if (forcet .gt.0) then !force is increasing theta so upwind derivative stencil is dmum
          muterm = (Amu(i,j,k) + Bmu(i,j,k)*nflux(k)) * dmum(:,j) + Amu2(i,j,k)*dmu2(:,j)
        else ! force is decreasing theta so upwind derivative stencil is dmup
          muterm = (Amu(i,j,k) + Bmu(i,j,k)*nflux(k)) * dmup(:,j) + Amu2(i,j,k)*dmu2(:,j)
        endif
        if (mum(j).gt.0) then ! traveling in the +z direction so upwind stencil is dzm
          zterm(-2:0) = Az(i,j) * dzm(:,k)
          rhsd(i) = - sum(f(i,j,k-2:k-1)*zterm(-2:-1))
        else ! traveling in the -z direction so upwind stencil is dzp
          zterm(0:2) = Az(i,j) * dzp(:,k)
          rhsd(i) = - sum(f(i,j,k+1:k+2)*zterm(1:2))
        endif
        rhsd(i) = rhsd(i) -sum(f(i,j-2:j-1,k)*muterm(-2:-1)) - sum(f(i,j+1:j+2,k)*muterm(1:2))  

        diag(i,0) = diag(i,0) + muterm(0) + zterm(0) + Af(i,j,k) +Bf(i,j,k)*nflux(k)
      enddo
      call pentdiag(diag,rhsd, newf(0:nE+1,j,k),nE+2)
    enddo; enddo
    where(newf.lt.0) 
      newf = 0
    endwhere
    f(0:nE+1,1:nmu,1:nz) = (1-tt)*f(0:nE+1,1:nmu,1:nz) + tt*newf
    call updateBC()
  ENDSUBROUTINE

  SUBROUTINE SweepMu()
    use helpers
    use options
    use parmpi
    use beam
    implicit none
    integer i,j,k
    double precision rhsd(0:nmu+1), diag(0:nmu+1,-2:2), pterm(-2:2), zterm(-2:2), forcep, forcet
    double precision tt, newf(nE, 0:nmu+1, nz)

    tt = implicit_theta
    rhsd(0) = 0
    diag(0,:) = 0
    diag(0,0) = 1
    diag(0,1) = -1

    rhsd(nmu+1) = 0
    diag(nmu+1,:) = 0d0
    diag(nmu+1,0) = 1d0
    diag(nmu+1,-1) = -1d0
    do k = 1, nz; do i = nE, 1,-1
      do j = 1, nmu
        forcep = (Ap(i,j,k)+Bp(i,j,k)*nflux(k))*dpm(0,i) ! Current df/dp component to force in the p-hat direction
        forcet = (Amu(i,j,k)+Bmu(i,j,k)*nflux(k))*dmum(0,j) ! Current df/dmu component to force in the theta-hat direction

        if (forcet .gt.0) then !force is increasing theta so upwind derivative stencil is dmum
          diag(j,:) = (Amu(i,j,k) + Bmu(i,j,k)*nflux(k)) * dmum(:,j) + Amu2(i,j,k)*dmu2(:,j)
        else !force is decreasing theta so upwind derivative stencil is dmup
          diag(j,:) = (Amu(i,j,k) + Bmu(i,j,k)*nflux(k)) * dmup(:,j) + Amu2(i,j,k)*dmu2(:,j)
        endif
        if (forcep .gt.0) then ! acceleration so upwind derivative stencil is dpm
          pterm = (Ap(i,j,k) + Bp(i,j,k)*nflux(k)) * dpm(:,i) + Ap2(i,k)*dp2(:,i)
          if (i.eq.1.and.par_leftedge(1)) pterm = 0
        else ! deceleration so upwind derivative stencil is dpp
          pterm = (Ap(i,j,k) + Bp(i,j,k)*nflux(k)) * dpp(:,i) + Ap2(i,k)*dp2(:,i)
          if (i.eq.nE.and.par_rightedge(1)) pterm =0
        endif
        if (mum(j).gt.0) then ! traveling in the +z direction so upwind stencil is dzm
          zterm(-2:0) = Az(i,j)*dzm(:,k)
          rhsd(j) = - sum(f(i,j,k-2:k-1)*zterm(-2:-1))
        else ! traveling in the -z direction so upwind stencil is dzp
          zterm(0:2) = Az(i,j)*dzp(:,k)
          rhsd(j) = - sum(f(i,j,k+1:k+2)*zterm(1:2))
        endif
        diag(j,0) = diag(j,0) + Af(i,j,k) + pterm(0) + zterm(0) + Bf(i,j,k)*nflux(k)
        rhsd(j) = rhsd(j) -sum(f(i-2:i-1,j,k)*pterm(-2:-1)) - sum(f(i+1:i+2,j,k)*pterm(1:2))
      enddo
      call pentdiag(diag, rhsd, newf(i,0:nmu+1,k), nmu+2)
    enddo; enddo
    where(newf.lt.0)
      newf = 0
    endwhere
    f(1:nE,0:nmu+1,1:nz) = tt*newf + (1-tt)*f(1:nE,0:nmu+1,1:nz)
    call updateBC()
  ENDSUBROUTINE

  SUBROUTINE SweepZ()
    use helpers
    use options
    use parmpi
    use beam
    implicit none
    integer i,j,k
    double precision rhsd(0:nz+1), diag(0:nz+1,-2:2), pterm(-2:2), muterm(-2:2), forcep, forcet, newf(nE,nmu, 0:nz+1), tt

    diag(0,:) = 0
    diag(0,0) = 1

    rhsd(nz+1) = 0
    diag(nz+1,:) = 0d0
    diag(nz+1,0) = 1d0
    tt = implicit_theta
    do j = 1, nmu; do i = nE, 1,-1
      rhsd(0) = injectedf(i,j)
      do k = 1, nz
        forcep = (Ap(i,j,k) + Bp(i,j,k)*nflux(k))*dpm(0,i) ! Current df/dp component of the force in the p-hat direction
        forcet = (Amu(i,j,k) + Bmu(i,j,k)*nflux(k))*dmum(0,j) ! Current df/dmu component of the force in the theta-hat direction

        if (mum(j).gt.0) then ! traveling in the +z direction so upwind stencil is dzm
          diag(k,-2:0) = Az(i,j)*dzm(:,k) ; diag(k,1:2) = 0
        else ! traveling in the -z direction so upwind stencil is dzp
          diag(k,0:2) = Az(i,j)*dzp(:,k); diag(k,-2:-1) = 0
        endif          
        if (forcet .gt.0) then !force is increasing theta so upwind derivative stencil is dmum
          muterm = (Amu(i,j,k) + Bmu(i,j,k)*nflux(k)) * dmum(:,j) + Amu2(i,j,k)*dmu2(:,j)
          rhsd(k) = -sum(f(i,j-2:j-1,k)*muterm(-2:-1)) - sum(f(i,j+1:j+2,k)*muterm(1:2))
        else ! force is decreasing theta so upwind derivative stencil is dmup
          muterm = (Amu(i,j,k) + Bmu(i,j,k)*nflux(k)) * dmup(:,j) + Amu2(i,j,k)*dmu2(:,j)
          rhsd(k) = -sum(f(i,j+1:j+2,k)*muterm(1:2)) - sum(f(i,j-2:j-1,k)*muterm(-2:-1))
        endif
        if (forcep .gt.0) then ! acceleration so upwind derivative stencil is dpm
          pterm = (Ap(i,j,k) + Bp(i,j,k)*nflux(k)) * dpm(:,i) + Ap2(i,k)*dp2(:,i)
          if (i.eq.1.and.par_leftedge(1)) pterm = 0
          rhsd(k) = rhsd(k) - sum(f(i-2:i-1,j,k)*pterm(-2:-1)) - sum(f(i+1:i+2,j,k)*pterm(1:2))
        else ! deceleration so upwind derivative stencil is dpp
          pterm = (Ap(i,j,k) + Bp(i,j,k)*nflux(k)) * dpp(:,i) + Ap2(i,k)*dp2(:,i)
          if (i.eq.nE.and.par_rightedge(1)) pterm = 0
          rhsd(k) = rhsd(k) - sum(f(i-2:i-1,j,k)*pterm(-2:-1)) - sum(f(i+1:i+2,j,k)*pterm(1:2))
        endif
        diag(k,0) = diag(k,0) + Af(i,j,k) +  muterm(0) + pterm(0) + Bf(i,j,k)*nflux(k)
      enddo
      call pentdiag(diag,rhsd, newf(i,j,0:nz+1),nz+2)
    enddo; enddo
    where (newf.lt.0)
      newf = 0
    endwhere
    f(1:nE,1:nmu,0:nz+1) = tt*newf + (1-tt)*f(1:nE,1:nmu,0:nz+1)
    call updateBC()
  ENDSUBROUTINE

  SUBROUTINE analytic_soln()
    use parmpi
    use beam
    implicit none
    double precision pz(nE,nmu,nz), pp, injectedf_all(nE_all, nmu), pm_all(nE_all)
    integer i,j,k, ii

    do j = 1, nmu
      call par_Gather1D(injectedf(:,j),injectedf_all(:,j),nE,1)
    enddo
    call par_Gather1D(pm(1:nE), pm_all, nE, 1)

    do k = 1, nz; do j=1, nmu; do i = 1, nE
      pz(i,j,k) = (Ap(i,j,k)/Az(i,j))*(zm(k)-zm(k-1))*(clight/mbeam)
    enddo; enddo; enddo
    call par_CumSum3(pz,nE,nmu,nz)
    do k = 1, nz; do j=1, nmu; do i = 1, nE
      pp = pm(i) - pz(i,j,k)
      ii = minloc(abs(pp-pm_all(1:nE)),1)
      if (pp.ge.pm_all(1).and.pp.le.pm_all(nE_all)) then
        if (pp.lt.pm_all(ii)) ii=ii-1
        f(i,j,k) = ((injectedf_all(ii+1,j)-injectedf_all(ii,j))*pp + injectedf_all(ii,j)*pm_all(ii+1) & 
                   - injectedf_all(ii+1,j)*pm_all(ii))/(pm_all(ii+1)-pm_all(ii))
      else
        f(i,j,k) =0
      endif
    enddo; enddo; enddo
    call updateBC()
  ENDSUBROUTINE

  FUNCTION residual() result(r)
    use helpers
    use beam
    use parmpi
    implicit none
    integer i,j,k,ierr
    double precision a0, ap, am, az, absres, d(2), r

    absres = 0d0
    do k = 1, nz; do j = 1, nmu; do i = 1, nE 
      a0 = (A(0,i,j,k) + B(0,i,j,k)*nflux(k))*f(i,j,k)
      ap = (A(1,i,j,k) + B(1,i,j,k)*nflux(k))*f(i-1,j,k) + (A(2,i,j,k) + B(2,i,j,k)*nflux(k))*f(i+1,j,k)
      am = (A(3,i,j,k) + B(3,i,j,k)*nflux(k))*f(i,j-1,k) + (A(4,i,j,k) + B(4,i,j,k)*nflux(k))*f(i,j+1,k)
      az = A(5,i,j,k)*f(i,j,k-1) + A(6,i,j,k)*f(i,j,k+1) 
      res(i,j,k) = a0 + ap + am + az
      absres = absres + ((abs(a0)+abs(ap)+abs(am)+abs(az)))**2
    enddo; enddo; enddo
    d = (/sum(res**2),absres/)
    call MPI_Allreduce(MPI_IN_PLACE, d, 2, MPI_DOUBLE_PRECISION, MPI_SUM, cart_comm, ierr)
    r= sqrt(d(1)/d(2))
  ENDFUNCTION

  FUNCTION IsConverged(r,d) result(conv)
    use options
    implicit none
    double precision r,d
    logical conv
    if (r.le.tolres.and.d.le.toldiff) then
      conv = .true.
    else
      conv = .false.
    endif
    return
  ENDFUNCTION
  
  SUBROUTINE Calcxi(p, Thetab, xb, xi0,xi1,xip, inc_rel)
    implicit none
    integer,parameter :: nn = 1000
    double precision p,xi0,xi1,xip, Thetab, bsk, xb,u, gma, gmap(nn),up(nn),integ(nn-1), L0,L1
    double precision, external :: bessk
    logical inc_rel
    integer ii

    xb = .5*p**2/Thetab
    if (.not.inc_rel) then ! Classical implementation from Trubnikov 1965
      xip = sqrtpi2*exp(-xb)*sqrt(xb)
      xi1 = erf(sqrt(xb)) - xip
      xi0 = xi1
    else !Relativistic implementation from Pike & Rose 2014 (PHYSICAL REVIEW E 89, 053107 (2014))
      bsk = bessk(2,1d0/Thetab)
      u = p * clight
      gma = sqrt(p**2+1)
      do ii = 1, nn
        up(ii) = (ii - 1d0)/(nn-1d0)*u
      enddo
      gmap = sqrt((up/clight)**2+1)
      integ = exp(-.5*(gmap(1:nn-1)+gmap(2:nn))/Thetab) * (up(2:nn) - up(1:nn-1))
      L1 = sum(integ)
      integ = integ / (.5*(gmap(1:nn-1)+gmap(2:nn)))
      L0 = sum(integ)
      if (bsk.gt.0) then
        xi1 = (gma**2*L1 - Thetab*L0 + (Thetab-gma)*u*exp(-gma/Thetab)) / (bsk*clight)
        xi0 = (gma**2*L0 - Thetab*L1 + (Thetab-gma)*u*exp(-gma/Thetab)) / (bsk*clight)
        xip = (2*gma*L1 + (1/Thetab+2*Thetab)*u*exp(-gma/Thetab)) / (bsk*clight)
      else
        xi1 = gma**2*(erf(sqrt(xb)) - sqrtpi2*exp(-xb)*sqrt(xb))
        xip = 2*gma*(erf(sqrt(xb)) - sqrtpi2*exp(-xb)*sqrt(xb)) + gma**3/Thetab*sqrtpi2*exp(-xb)*sqrt(xb)
        xi0 = xi1
      endif
    endif
  ENDSUBROUTINE

  FUNCTION CoulogEN(b,g,E,EnI)
    implicit none
    double precision b,g, Eni, CoulogEN, E

    CoulogEN = max(log(b*g*sqrt(E) * me / EnI) - .5*b**2,0d0) ! If E < EnI this should be zero
    return
  ENDFUNCTION

  FUNCTION CoulogHN(b,g,EnI)
    implicit none
    double precision b,g, Eni, CoulogHN

    CoulogHN = max(log(2*b**2*g**2* me / EnI) - b**2,0d0) ! If E < EnI this should be zero
    return
  ENDFUNCTION

  FUNCTION CoulogND(b,g,znn)
    implicit none
    double precision b,g, znn, CoulogND

    CoulogND = max(.5*log((g*b)**2/znn**(2./3.)/2/alpha2),0d0)
    return
  ENDFUNCTION

  FUNCTION CoulogCC(mbeam,Zbeam, mi, Zi, ni, xi, b)
    implicit none
    double precision mbeam,Zbeam,mi,Zi,ni,xi,b, u, redm, CoulogCC,rmin,rmax

    u = (1-1/sqrt(xi))*b
    redm = mi*mbeam/(mi + mbeam)
    rmin = max(e2 * abs(Zi * Zbeam)/redm/u/u, hbc2/u/redm)
    rmax = sqrt(mi/ni)*abs(Zi)*b*sqrtpie2
    CoulogCC = max(log(rmax/rmin),0d0)
    return
  ENDFUNCTION

  FUNCTION CalcEtaIon(mi, Zi, ni, ne, T)
    implicit none
    double precision mi, Zi, ni, ne, T, FZ, CalcEtaIon, xi, b, cl
    
    xi = 1.5d0 * mi / me
    b = sqrt(1d0 - 1d0/(xi*kb*T/me +1)**2)
    cl = CoulogCC(me, -1d0, mi, Zi, ni, xi, b)
    FZ = (1d0 + 1.198d0*Zi + 0.222*Zi**2)/ (1d0 + 2.966d0*Zi + 0.753d0*Zi**2) !S.P. Hirshman, Phys. Fluids 20, 589 (1977)
    CalcEtaIon = 4*sqrt(2*pi)/3*Zi*e2*sqrt(me)*cl*FZ/clight/(kb*T)**1.5 * ni / ne
    return
  ENDFUNCTION

  FUNCTION CalcEtaNeut(nn, ne, T)
    implicit none
    double precision nn, ne, T, CalcEtaNeut

    CalcEtaNeut = 5.2d-11*3.759d-6*1d6 * me/e2/clight**2 * nn/ne * sqrt(T)! FROM MARTINEZ-SYKORA 2012 EQ 23-24
    return
  ENDFUNCTION
ENDMODULE

