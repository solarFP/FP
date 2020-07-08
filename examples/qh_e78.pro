FUNCTION qh_e78, eckev, dlt, injflx, mu0, atmfile, proton = proton
; Implements Hawley & Fisher (1994) modification to Emslie (1978) thick target heating rate

@const
ec = eckev * 1d3
readatm,atmfile
ea = (dlt-1)/(dlt-2)*ec /5; avg E
if (keyword_set(proton)) then begin
  ma = !mp
endif else begin
  ma = !me
endelse

bta = sqrt(1d0 - 1d0/(ea / ma + 1)^2)
redm = !me*ma / (!me + ma)
hbc2     = 6.1992097d-05
sqrtpie2 = sqrt(!dpi/!e2) 
alph = 0.0072973525693d0 ;fine structure
Kconst = 2*!dpi*!e2^2

tothy = dni[1,*]+dnn[0,*] ; total # hydrogen (neutral + protons)
ndep = n_elements(zin)
cl = dblarr(ndep)
for k = 0, ndep-1 do begin
  xi = !me * ea / !kb/tg[k]
  u = (1-1/sqrt(xi))*bta
  rmin = !e2/redm/u/u > hbc2/u/redm
  rmax = sqrt(!me/dni[0,k])*bta*sqrtpie2
  cl[k] = alog(rmax/rmin) > 0d0
endfor

clp = alog(bta*(ea/ma+1)*sqrt(Ea/ma) * !me / Enion[0]) - .5*bta^2
if (keyword_set(proton)) then clp = alog(2*bta^2*(Ea/ma+1)^2* !me / Enion[0]) - bta^2
clpp = .5*alog(((ea/ma+1)*bta)^2/zn[0]^(2./3.)/2/alph^2)
xion = dni[1,*]/tothy

gam = ma/!me*(xion*cl + (1-xion)*clp)
if (keyword_set(proton)) then begin
  bb = (xion-1)*clp/(clp + xion*(cl-clp))
endif else begin
  bb = (2*xion*cl + (1-xion)*clpp) / (clp + xion*(cl-clp))
endelse

N = total(tothy*[0,(zin-zin[1:*])],/cum) ; column depth
eqN = total(gam / cl * [0,(N[1:*]-N)],/cum) ; equivalent depth (HF94 discussion of eq 2.10) 
eqNc = mu0*Ec^2/(2+.5*bb)/cl/Kconst
b = beta(dlt/2,2d0/(4d0 + bb))
ind = where(eqN/eqNc lt 1)
if (ind[0] ne -1) then begin
  ib = ibeta(dlt/2,2d0/(4d0 + bb[ind]), eqN[ind]/EqNc[ind])
  b[ind] *= ib
endif

qh = .5 * Kconst/mu0 * gam * (dlt-2) * b*injflx/Ec^2*(eqN/eqNc)^(-dlt/2) * tothy ; HF94 eq 2.10
qh = reform(qh)
qh[0] = qh[1]
return,qh
END
