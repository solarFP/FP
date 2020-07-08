PRO makefpoutstr,fpout,atmfile
@const
nz = fpout.inputparams.nz
nen = fpout.inputparams.nen
nmu = fpout.inputparams.nmu

Em = .5*(fpout.E[1:*]+fpout.E)
mbeam = fpout.inputparams.mbeam

if (fpout.inputparams.inc_relativity) then begin
  gma = fpout.E/(mbeam/1d3) + 1d0
  bta = sqrt(1d0 - 1d0/gma^2)
  gmam = Em/(mbeam/1d3) + 1d0
  btam = sqrt(1d0 - 1d0/gmam^2)
endif else begin
  gma = replicate(1d0,nen+1)
  bta = sqrt(2*fpout.E*1d3/mbeam)
  gmam = replicate(1d0,nen)
  btam = sqrt(2*Em*1d3/mbeam)
endelse

theta = acos(fpout.mu)
thetam = fpout.inputparams.oneD ? theta[0] : .5*(theta[1:*]+theta)
mum = cos(thetam)

zm = .5*(fpout.z[1:*]+fpout.z)

;Calculate useful fluxes
vz = !clight * btam # mum
vzv = rebin(reform(vz,[nen,nmu,1]),[nen,nmu,nz])
vx = !clight * btam # sin(thetam)
vxv = rebin(reform(vx,[nen,nmu,1]),[nen,nmu,nz])
vv = rebin(reform(btam*!clight,[nen,1,1]),[nen,nmu,nz])
esvolv = rebin(reform(fpout.esvol,[nen,nmu,1]),[nen,nmu,nz])
ee = rebin(reform(em*1d3*!ergperev,[nen,1,1]),[nen,nmu,nz])
domega = fpout.esvol / rebin(reform(fpout.E[1:nen] - fpout.E[0:nen-1],[nen,1]),[nen,nmu]) ; solid angle part of esvol
domega = rebin(reform(domega,[nen,nmu,1]),[nen,nmu,nz])

ef = total(total(vzv * fpout.f * esvolv*ee,1),1) ; Energy flux in z direction (erg cm^-2 s^-1)
efx = total(total(vxv * fpout.f * esvolv*ee,1),1); Energy flux in x direction (erg cm^-2 s^-1)
nflux = total(total(fpout.f*vzv*esvolv,1),1); number flux in z direction (particles cm^-2 s^-1)
flux = total(fpout.f*vzv*domega,2) ; number flux distribution in z direction (particles cm^-2 s^-1 keV^-1)
flx = total(fpout.f*vv*domega,2) ; number flux speed distribution. (particles cm^-2 s^-1 keV^-1)

if (abs(mbeam-!me) lt 1d3) then begin ;electron beam
  eph = Em
  fpbrem,em,flx,eph,brem,fpout.atm
  totbrem = eph
  for i = 0, n_elements(eph) -1 do totbrem[i] = total((fpout.z[1:*]-fpout.z)*brem[*,i])
endif else begin
  ; Not electron beam so do not calculate bremsstrahlung
  brem = 0d0
  totbrem = 0d0
  eph = 0d0
endelse

fpout = create_struct(fpout,'nen',nen,'nmu',nmu,'nz',nz,'gma',gma, 'bta',bta,'Em',Em, 'gmam',gmam, 'btam', btam, 'theta',theta, $
                      'thetam', thetam, 'mum',mum, 'zm',zm, 'eflux',ef, 'efluxx', efx, $
                      'nflux',nflux, 'flx',flx, 'flux',flux,'eph',eph, 'brem', brem, 'totbrem',totbrem, 'atmfile',atmfile)
END
