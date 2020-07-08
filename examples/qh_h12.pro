FUNCTION CoulogCC,mbeam,Zbeam, mi, Zi, ni, xi, b
  sqrtpie2 = sqrt(!dpi/!e2)
  hbc2     = 6.1992097d-05       ; hbar c /2 in units of eV cm
  u = (1-1/sqrt(xi))*b
  redm = mi*mbeam/(mi + mbeam)
  rmin = !e2 * abs(Zi * Zbeam)/redm/u/u > hbc2/u/redm
  rmax = sqrt(mi/ni)*abs(Zi)*b*sqrtpie2
  cl = alog(rmax/rmin) >0
  return, cl
END

FUNCTION CalcEtaIon,mi, Zi, ni, nel, T
  xi = 1.5d0 * mi / !me
  b = sqrt(1d0 - 1d0/(xi*!kb*T/!me +1)^2)
  cl = CoulogCC(!me, -1d0, mi, Zi, ni, xi, b)
  FZ = (1d0 + 1.198d0*Zi + 0.222*Zi^2)/ (1d0 + 2.966d0*Zi + 0.753d0*Zi^2) ;S.P. Hirshman, Phys. Fluids 20, 589 (1977)
  EtaIon = 4*sqrt(2*!dpi)/3*Zi*!e2*sqrt(!me)*cl*FZ/!clight/(!kb*T)^1.5 * ni / nel
  return, EtaIon
END

FUNCTION CalcEtaNeut,nn, nel, T
  EtaNeut = 5.2d-11*3.759d-6*1d6 * !me/!e2/!clight^2 * nn/nel * sqrt(T); FROM MARTINEZ-SYKORA 2012 EQ 23-24
  return, EtaNeut
END

FUNCTION CalcEta,mion,Zion,dni,dnn,tg
  nion = n_elements(mion)
  nneut = n_elements(dnn[*,0])

  nel = reform(dni[0,*])
  eta = dblarr(n_elements(tg))
  for i = 1, nion-1 do begin 
    ni = reform(dni[i,*])
    eta += CalcEtaIon(mion[i],Zion[i],ni,nel,tg)
  endfor
  nn = reform(dnn[0,*])
  eta += CalcEtaNeut(nn,nel,tg)
  return, eta
END

FUNCTION qh_h12,eckev, dlt, injflx, Ethkev, atmfile, ee = ee, flux = flux
; implements Holman 2012 return current heating rate (Eq. 24)
@const
Ec = Eckev * 1d3 ; eV
Eth = Ethkev * 1d3; eV

Fe0 = injflx * (dlt-2)/(dlt-1)/(Ec * !ergperev) ; e- /cm^2/s
;print,'Fe0:',fe0
readatm,atmfile
eta = CalcEta(mions,zion,dni,dnn,tg)
save,file = 'eta.sav',eta

rhoc = (Ec - Eth)/(!e2*Fe0) ; Eq. 19
x = zin[0]-zin
dx = [0,zin - zin[1:*]]
rho = total(eta*dx,/cum)
xc = interpol(x,rho,rhoc)
etac = interpol(eta,rho,rhoc)
;insert xc into x
;ii = whereclose_less(x,xc)
;x = [x[0:ii],xc,x[ii+1:*]]
;rho = [rho[0:ii],rhoc,rho[ii+1:*]]
;eta = [eta[0:ii],etac,eta[ii+1:*]]

nx = n_elements(x)
V = dblarr(nx)
ind = where(rho le rhoc, comp = cind)
V[ind] = !e2*Fe0*rho[ind] ;Eq 18
if (cind[0] ne -1) then V[cind] = Ec*(dlt*!e2*Fe0/Ec*(rho[cind]-rhoc) + 1)^(1/dlt) - Eth

;setup energy grid
nen = 500
ee = rebin(reform(exp(dindgen(nen)/(nen-1) * (alog(2d3*Ec) - alog(Eth/2)) + alog(Eth/2)),[nen,1]),[nen,nx])
Flux = dblarr(nen,nx)
vv = rebin(reform(v,[1,nx]),[nen,nx])
ind = where(ee ge Ec - vv , comp = cind)
Flux[ind] = (dlt-1)*Fe0/Eckev * ((ee[ind] + vv[ind])/Ec)^(-dlt); Eq 12
if (cind[0] ne -1) then Flux[cind]=0
ee = reform(ee[*,0])/1d3 ; convert to keV

q = dblarr(nx)
ind = where(rho lt rhoc, comp = cind)
Q[ind] = eta[ind] * !e2 * Fe0^2 ; Eq. 24
if (cind[0] ne-1) then Q[cind] = eta[cind] * !e2 * Fe0^2*(dlt*Eth/Ec + V[cind]/Ec)*(Eth/Ec+V[cind]/Ec)^(1-2*dlt)
Q*= !ergperev ; erg /cm^3/s
return, q
END
