MODULE loop
! This module encapsulates the loop stratification.
! Module Variables:
!  nIon:	The number of ion species in the ambient plasma.
!  nNeutral:	The number of neutral species in the ambient plasma.
!  zin: 	The loop z-axis measured from the photosphere (cm). 
!  tg:		The gas temperature (K).
!  bfield: 	The magnetic flux density (G).
!  dni: 	The number densities of the ion species composing the ambient plasma (cm 3).
!  dnn: 	The number densities of the neutral species composing the ambient plasma (cm-3).
!  mion: 	The masses of the ion species (eV).
!  Zion:	The charges of the ion species (elementary charge).
!  Zn:		Atomic number of the neutral species.
!  Enion: 	Ionization energy of the neutral species (eV).
!
  implicit none
  double precision,allocatable :: tg(:), bfield(:), dni(:,:), dnn(:,:), mion(:), Zion(:), Zn(:), Enion(:), zin(:)
  double precision, allocatable :: zin_all(:), tg_all(:), bfield_all(:),dni_all(:,:), dnn_all(:,:)
  integer Nion, Nneutral

  contains
   
  SUBROUTINE init_loop()
    use grid
    implicit none
    allocate(zin_all(nz_all), tg_all(nz_all), bfield_all(nz_all), dni_all(Nion,nz_all), dnn_all(Nneutral,nz_all), &
             mion(Nion), Zion(Nion), Zn(Nneutral), Enion(Nneutral)) 
  ENDSUBROUTINE

  SUBROUTINE getmyatm()
    use parmpi
    use grid
    implicit none
    double precision bfield_ref(-1:nz_all+2)

    ! limit the full nz_all arrays to just the parts needed by this process
    allocate(zin(nz), tg(nz), bfield(-1:nz+2), dni(Nion,nz), dnn(Nneutral, nz))
    bfield_ref(1:nz_all) = bfield_all
    bfield_ref(-1) = bfield_all(2)
    bfield_ref(0) = bfield_all(1)
    bfield_ref(nz_all+1:nz_all+2) = bfield_ref(nz_all)

    zin = zin_all(fulloffset(3)+1:fulloffset(3)+nz)
    bfield = bfield_ref(fulloffset(3)-1:fulloffset(3)+nz+2)
    tg = tg_all(fulloffset(3)+1:fulloffset(3)+nz)
    dni = dni_all(:,fulloffset(3)+1:fulloffset(3)+nz)
    dnn = dnn_all(:,fulloffset(3)+1:fulloffset(3)+nz)
  ENDSUBROUTINE
     
  SUBROUTINE deallocate_loop()
    implicit none
    deallocate(zin, tg, bfield, dni, dnn, mion, Zion, Zn, Enion)
    deallocate(zin_all,tg_all,bfield_all,dni_all,dnn_all)
  ENDSUBROUTINE
ENDMODULE
