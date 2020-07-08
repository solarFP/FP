;RUNS FP for proton beam injection
; Compares to warm target of Tamres et al 1986 and cold target of Emslie 1978
; 
;  Return current, magnetic mirroring and synchrotron forces are OFF

; in IDL session run: @compare_fp_t86

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
doeps = 0 ;SET THIS TO 1 TO MAKE AN EPS PLOT

;SET PATH TO INCLUDE FP/idl
!path = '../idl:'+!path

; Ec = 100 keV
eflux = 1d11
dlt = 4d0
Ecut1 = 1d2
atm = 'atm.13Mm.3MK.onlyH.dat'
fp,fpout1,eflux = eflux, dlt = dlt, ecut = ecut1, inc_rc = 0, inc_magmirror = 0, inc_synchro = 0, mbeam = !mp, zbeam = 1d0, atm = atm, maxit = 200, impl=.95, nen = 300, nmu = 100, patype = 0, emin = 1d1, emax = 1d6
fp,fpout1_all,eflux = eflux, dlt = dlt, ecut = ecut1, inc_rc = 1, inc_magmirror = 1, inc_synchro = 1, mbeam = !mp, zbeam = 1d0, atm = atm, maxit = 200, impl=.95, nen = 300, nmu = 100, patype = 0, emin = 1d1, emax = 1d6
qhwt1 = qh_warmtarget(ecut1,dlt,eflux,1d0, atm, /proton)
qhct1 = qh_e78(ecut1,dlt,eflux,1d0, atm, /proton)

; Ec = 1 MeV
eflux = 1d11
dlt = 4d0
Ecut2 = 1d3
atm = 'atm.13Mm.3MK.onlyH.dat'
fp,fpout2,eflux = eflux, dlt = dlt, ecut = ecut2, inc_rc = 0, inc_magmirror = 0, inc_synchro = 0, mbeam = !mp, zbeam = 1d0, atm = atm, maxit = 200, impl=.95, nen = 300, nmu = 100, patype = 0, emin = 1d1, emax = 1d6
fp,fpout2_all,eflux = eflux, dlt = dlt, ecut = ecut2, inc_rc  = 1, inc_magmirror = 1, inc_synchro = 1, mbeam = !mp, zbeam = 1d0, atm = atm, maxit = 200, impl=.95, nen = 300, nmu = 100, patype = 0, emin = 1d1, emax = 1d6
qhwt2 = qh_warmtarget(ecut2,dlt,eflux,1d0, atm, /proton)
qhct2 = qh_e78(ecut2,dlt,eflux,1d0, atm, /proton)

;save,file = 'fpout.proton.sav',fpout1,qhwt1, qhct1, fpout2, qhwt2, qhct2, fpout1_all, fpout2_all
;restore,'fpout.proton.sav'

; Now plot the results
olddev = !d.name
oldfont = !p.font
if (doeps) then begin set_plot,'ps' & !p.font =0 & device,file='compare_fp_t86.eps',/enc,/color,/helv,/isolatin1, bits_per_pixel=8,xsize = 10.0,ysize = 8.0,/inches
loadct,39
if (~doeps) then device,decom = 0
red = 250 & blue = 64 & orange = 210 & white = 255 & black =0
color0 = white
if (doeps) then color0 = black
pm = !p.multi
!p.multi = [0,1,2]
cs1 = 1.5
cs2 = 2.0
xm = [10,2]
if (doeps) then xm = [8,2]
thick = 3
if (doeps) then thick = 5

plot,fpout1.zm/1d8,fpout1.heatrate,/ylog, yr = [1d0,1.4d4], thick = thick,  charsize = cs1,/ysty, xmargin = xm, ymargin = [0,2], xtickname = replicate(' ',30),/xsty, ytickformat = 'tickexp', tit = 'Proton injection in CL'
oplot,fpout1_all.zm/1d8,fpout1_all.heatrate,thick = thick, color = blue
oplot,fpout1.zm/1d8,qhwt1,thick = thick, color = orange
oplot,fpout1.zm/1d8,qhct1,thick = thick, color = red
legend,['FP_CC_only','FP','E78','T86'], thick = replicate(thick,4d0), charsize = cs1, color = [color0, blue, red, orange], linesty = [0,0,0,0]
legend,['!9d!x: '+strng(dlt,format = '(F5.1)'), 'E!Dc!N: '+strng(Ecut1,form = '(F6.1)') +' keV', 'F!D0!N: '+strng(eflux,form='(E10.1)') + ' erg cm!E-2!N s!E-1!N'], box= 0, charsize = cs1, pos = [5.5,7000]

plot,fpout2.zm/1d8,fpout2.heatrate,/ylog, yr = [1d0,1.4d4], thick = thick,  charsize = cs1,/ysty, xmargin = xm, ymargin = [4,0],/xsty, xtit = 'Distance along loop (Mm)', ytickformat= 'tickexp'
oplot,fpout2_all.zm/1d8,fpout2_all.heatrate,thick = thick, color = blue
oplot,fpout2.zm/1d8,qhwt2,thick = thick, color = orange
oplot,fpout2.zm/1d8,qhct2,thick = thick, color = red
legend,['FP_CC_only','FP','E78','T86'], thick = replicate(thick,4d0), charsize = cs1, color = [color0, blue, red, orange], linesty = [0,0,0,0]
legend,['!9d!x: '+strng(dlt,format = '(F5.1)'), 'E!Dc!N: '+strng(Ecut2/1d3,form = '(F6.1)') +' MeV', 'F!D0!N: '+strng(eflux,form='(E10.1)') + ' erg cm!E-2!N s!E-1!N'], box= 0, charsize = cs1, pos = [5.5,7000]

xyouts,.03,.5,'Heating Rate (erg cm!E-3!N s!E-1!N)', charsize = cs2, align = .5, orient = 90,/norm
if (doeps) then device,/close
;RESET VARIABLES CHANGED IN PLOTTING
set_plot,olddev
!P.font = oldfont
!p.multi = pm
