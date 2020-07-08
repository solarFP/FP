; Computes X-ray bremsstrahlung from the electron flux distribution
PRO fpbrem, e, flx, eph, xraybrem, atm, Zan = Zan
@const
mc2 = !me/1d3
; returns photons/cm^3/s/kev
;flx is the direction independent flux so flx = integrate(f*beta*c*dmu). units of flx are electrons/cm^2/s/keV
; photon spectrum (# of photons / cm^2 /s/ keV) = n_ambient * l * integral(dsigma/dk(E,k) * electron flux(E) * dE)
nee = n_elements(e)
neph = n_elements(eph)

;Average atomic number
if (n_elements(Zan) eq 0) then Zan = 1.2d0

ndep = n_elements(atm.zin)
xraybrem = dblarr(ndep,n_elements(eph))
dsigmadk = dblarr(n_elements(eph),n_elements(e))

for j = 0, n_elements(e) -1 do begin
  idx = where(eph lt e[j])
  if (idx[0] eq -1) then continue
  brm_bremcross,e[j], eph[idx], Zan,cross
  dsigmadk[idx,j] = cross /mc2 ; units of photons cm^2 / keV / incoming electron
endfor
nh = reform(atm.dnn[0,*] + atm.dni[1,*]) ; hydrogen number density
;calculate spectrum at each loop grid cell
for i=0,ndep -1 do begin
;  now integrate over e
   ;photons /cm^2 /cm/s/kev 
   xraybrem[i,*] = total( ( dsigmadk[*,1:*]*rebin(reform(flx[1:*,i],[1,nee-1]),[neph,nee-1]) $ 
                           +dsigmadk[*,0:nee-2]*rebin(reform(flx[0:nee-2,i],[1,nee-1]),[neph,nee-1]) $ 
                          )*.5 * rebin(reform(e[1:*]-e[0:nee-2],[1,nee-1]),[neph,nee-1]), 2) * nh[i]
endfor

end

