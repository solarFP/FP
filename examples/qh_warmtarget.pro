function nbyw,x
   common nbywcb,Ndepth, T, w0,w, mode, clog,mu0

   if mode eq 0 then w0 = x[0]
   if mode eq 1 then w = x[0]

   if (w lt 0) then w = complex(w)
   K = 2*!dpi*!e2^2
   c = (2d0/3d0)*sqrt(!me/!mp/!dpi) ; defined below Eq 33

   CapW = 1d0/(1d0+2*c*w^1.5)
   CapW0 = 1d0/(1d0+2*c*w0^1.5)
   thrd = 1d0/3d0
   thrd2 = 2d0/3d0
   bta = beta(thrd2,thrd2)
   retval = mu0 * (!kb*T)^2/c/clog/K *(sqrt(w0)*(1d0-w/w0*(CapW/CapW0)^thrd) - thrd2/CapW0^thrd/(2*c)^thrd2/sqrt(w0)*(ibeta(thrd2,thrd2,CapW)*bta- ibeta(thrd2,thrd2,CapW0)*bta)) -Ndepth
   retval = abs(retval)
   return, retval
end

Function Qh_warmtarget, eckev, dlt, eferg, mu0, atmfile, proton = proton
common nbywcb,Ndepth, T, ww0, ww, mode, clog, mu

@const
eflux = eferg/!ergperev
ecut = eckev * 1d3
readatm,atmfile
xi = enion[0]
ea = (dlt-1)/(dlt-2)*ecut /5; avg E
if (keyword_set(proton)) then begin
  ma = !mp
endif else begin
  ma = !me
endelse

nz = n_elements(zin)
qh = dblarr(nz)

mu = mu0
bta = sqrt(1d0 - 1d0/(ea / ma + 1)^2)
redm = !me*ma / (!me + ma)
hbc2     = 6.1992097d-05
sqrtpie2 = sqrt(!dpi/!e2)
alph = 0.0072973525693d0 ;fine structure
Kconst = 2*!dpi*!e2^2

tothy = reform(dni[1,*]+dnn[0,*]) ; total # hydrogen (neutral + protons)
xion = dni[1,*]/tothy
cl = dblarr(nz)
for k = 0, nz-1 do begin
  xi = !me * ea / !kb/tg[k]
  u = (1-1/sqrt(xi))*bta
  rmin = !e2/redm/u/u > hbc2/u/redm
  rmax = sqrt(!me/dni[0,k])*bta*sqrtpie2
  cl[k] = alog(rmax/rmin) > 0d0
endfor

clp = alog(bta*(ea/ma+1)*sqrt(Ea/ma) * !me / Enion[0]) - .5*bta^2
if (keyword_set(proton)) then clp = alog(2*bta^2*(Ea/ma+1)^2* !me / Enion[0]) - bta^2
clpp = .5*alog(((ea/ma+1)*bta)^2/zn[0]^(2./3.)/2/alph^2)
N = total(tothy*[0,(zin-zin[1:*])],/cum) ; column depth

c = (2d0/3d0)*sqrt(!me/!mp/!dpi) ; defined below Eq 33
for k = 0, nz-1 do begin
  ; numerically solve Tamres et al Eq. 38
  Ndepth = N[k]
  T = tg[k]
  wc = ecut / !kb/ T
  clog = cl[k]
  nw = 100
  mode=0 & ww = 1d0
  wm1 = (Ndepth * Kconst * clog * c /mu0/(!kb*T)^2 + 4.8838)^2
  wm = broyden(wm1,'nbyw')
;  print,Ndepth, wm1, wm, sqrt(4*cl[k]*Kconst*Ndepth/2/mu0)/!kb/T,clog,T,k
  if (wm lt wc) then wm = wc
  wh = 1d4*wm
  w0 = exp(dindgen(nw)/(nw-1) * (alog(wh)-alog(wm)) + alog(wm))
  mode =1
  integrand = dblarr(nw)
  if (wm gt wh) then stop
  for i =0, nw-1 do begin
    ; numerically solve Eq. 35 to get omega from N and omega0
    w = w0[i]/2
    ww0 = w0[i]
    w= broyden(w,'nbyw')
    w = w >1d0; set w to minimum of 1.0 
    Phi0 = eflux*(dlt-2)/ecut^2*(w0[i]/wc)^(-dlt)/(1d0-(wh/wc)^(2-dlt))
    integrand[i] = sqrt(w0[i])*(1+2*c*w^1.5d0)^(4d0/3d0)*Phi0/mu0/w^1.5/(1+2*c*w0[i]^1.5)^(1d0/3d0)
  endfor
  Qh[k] = Kconst*clog*total( .5*(integrand[1:*]+integrand[0:nw-2]) * (w0[1:*]-w0[0:nw-2])) *tothy[k]*!ergperev
endfor
print,tsum(zin,qh)
return,qh
END
