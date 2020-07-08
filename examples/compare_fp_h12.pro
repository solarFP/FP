;RUNS FP in mode to mimic Holman 2012 RC  model:
;  1D mode with only return currect on 
; Compares to 1.5D mode and including all forces. 

; in IDL session run: @compare_fp_h12

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
!path = '../idl:'+!path
@const

doeps = 0
dohighres = 0

; Injected beam parameters
dlt = 4d0 ; delta
ecut = 20d0 ; Cutoff energy in keV
eflux = 1d11 ; injected energy flux (erg cm^-2 s^-1)

; Loop model is 13 Mm half length with a 3.4 MK corona.
atmfile = 'atm.13Mm.3MK.dat'

readatm,atmfile
tt = tg[0]
Eth = (dlt)*!kb*tt/1d3

nen = dohighres ? 1000 : 200
; Call FP
fp,fpout1d,dlt = dlt, ecut = ecut, eflux=eflux, inc_rc = 1, inc_cc=0, reflecttop= 0, inc_magmirror = 0, inc_sync = 0, atmfile = atmfile, oneD = 1, $ 
   patype = 0, pasigma = .14, maxit = 400, impl = 1d0, nen=nen,nmu = 1, toldiff = 1d-4, tolres = 7d-3,emin = .05, emax = 2d3

;better to make the 1.5D cases patype =2. patype = 0 creates very large (actually infinite) fluxes in the theta direction
fp,fpout,dlt = dlt, ecut = ecut, eflux=eflux, inc_rc = 1, inc_cc=0, reflecttop= 0, inc_magmirror = 0, inc_sync = 0, atmfile = atmfile, oneD = 0, $
   patype = 2, pasigma = .14, maxit = 1000, impl = .8d0, nen=nen,nmu = 60, toldiff = 5d-3, tolres = 2d-3,emin = .05, emax = 2d3

fp,fpout_all,dlt = dlt, ecut = ecut, eflux=eflux, inc_rc = 1, inc_cc=1, reflecttop= 0, inc_magmirror = 1, inc_sync = 1, atmfile = atmfile, oneD = 0, $
   patype = 2, pasigma = .14, maxit = 800, impl = .5d0, nen=nen,nmu = 60, toldiff = 1d-2, tolres = 2d-3,emin = .05, emax = 2d3
;save,file = 'fpout.h12.sav',fpout,fpout1d,fpout_all

;restore,'fpout.h12.sav'
qh12 = qh_h12(Ecut,dlt,eflux, eth, atmfile,ee=ee,flux=flux)

device,decom =0
loadct,39
red = 250 & blue = 64 & orange = 210 & white = 255 & black = 0
color0 = doeps ? black : white

olddev = !d.name
oldfont = !p.font
if (doeps) then begin set_plot,'ps' & !p.font =0 & device,file='compare_fp_h12.eps',/enc,/color,/helv,/isolatin1, bits_per_pixel=8,xsize = 10.0,ysize = 8.0,/inches
thick = doeps ? 5 : 3
xm = doeps ? [7.8,1.5] : [10,2]
cs = doeps ? 1.5 : 1.0
!p.multi = [0,1,2]
plot,fpout1d.zm/1d8,fpout1d.heatrate,/ylog, thick = thick, xtit ='Distance along loop (Mm)',ytit = 'Heating rate (erg cm!E-3!N s!E-1!N)', yr= [1d1,1d5],/ysty, chars=cs, /xsty, xmar = xm
oplot,fpout1d.zm/1d8,qh12, co=red, thick = thick
oplot,fpout.zm/1d8,fpout.heatrate,co=orange,thick=thick
oplot,fpout_all.zm/1d8,fpout_all.heatrate,co=blue,thick=thick
legend,['FP_RC_only_1D','FP_RC_only','FP','RCCTTM'],linesty = [0,0,0,0], thick = thick, color = [color0,orange,blue,red], charsize = cs
thick = thick[0]
iz = 28
vline,fpout1d.zm[iz]/1d8, thick = thick, linesty = 1
sdelta = doeps ? '!9d!x' : '!7d!x'
xyouts,7.75,2.5d4,sdelta + ': '+strng(dlt,form = '(F6.1)'), charsize = cs
xyouts,7.75,1d4,'E!Dc!N: '+strng(Ecut,form='(F6.1)')+ ' keV', charsize = cs
xyouts,7.75,4d3,'F!D0!N: '+strng(eflux,form='(E12.1)') + ' erg cm!E-2!n s!E-1!N', charsize = cs

plot,fpout1d.em,fpout1d.flux[*,iz]/1d17,xr = [0,30],thick = thick, xtit = 'Energy (keV)', ytit = 'Flux (10!E17!N e!E-!N cm!E-2!N s!E-1!N keV!E-1!N)', yr = [-1d17,4d17]/1d17, chars= cs, xmar = xm
oplot,ee,flux[*,iz]/1d17,thick = thick, co = red
oplot,ee,flux[*,0]/1d17,thick = thick, co = color0, linesty = 2
oplot,fpout.em,fpout.flux[*,iz]/1d17, thick = thick, co = orange
oplot,fpout_all.em,fpout_all.flux[*,iz]/1d17,thick=thick, co = blue
legend,['FP_RC_only_1D', 'FP_RC_only', 'FP','RCCTTM', 'z = 0'], linesty = [0,0,0,0,2], thick = thick[0], color = [color0,orange,blue,red,color0], chars = cs

if (doeps) then device,/close
set_plot,olddev
!P.font = oldfont
