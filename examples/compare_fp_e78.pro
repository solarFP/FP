;RUNS FP in mode to mimic the cold thick target model:
;  Return current, magnetic mirroring and synchrotron forces are OFF

; in IDL session run: @compare_fp_e78

;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;;
doeps = 0 ;SET THIS TO 1 TO MAKE AN EPS PLOT

;SET PATH TO INCLUDE FP/idl
!path = '../idl:'+!path 

; Injected beam parameters
dlt = 4d0 ; delta
ecut = 20d0 ; Cutoff energy in keV
eflux = 1d11 ; injected energy flux (erg cm^-2 s^-1)

; Inject into a loop model that only includes H to most closely match the assumptions made in Emslie 1978
atmfile = 'atm.13Mm.3MK.onlyH.dat'
atmfile_all = 'atm.13Mm.3MK.dat'

; Call FP with all terms switched OFF except Coulomb Collisions
fp,fpout,dlt = dlt, ecut = ecut, eflux=eflux, inc_rc = 0, reflecttop= 1, inc_magmirror = 0, inc_sync = 0, atmfile = atmfile, patype = 2,pasigma = .1,maxit = 400, impl = .95d0, nen=100,nmu = 60, toldiff = 1d-3, tolres = 1d-3
; Call FP with all terms switched ON
fp,fpout_all,dlt = dlt, ecut = ecut, eflux=eflux, inc_rc = 1, reflecttop= 1, inc_magmirror = 1, inc_sync = 1, atmfile = atmfile_all, patype = 2,pasigma = .1,maxit = 500, impl = .9d0, nen=100,nmu = 60, toldiff = 1d-3, tolres = 1d-3
;save,file = 'fpout.sav',fpout, fpout_all
;restore,'fpout.sav'

; Now calculate the heating rate from the thick target model using Emslie 1978
; Construct a Gaussian PA Distribution to match what was used in FP
padist = exp(-.5*(fpout.thetam)^2/fpout.inputparams.pasigma^2)
ind = where(fpout.mum le 0)
padist[ind] =0
nrm = total(padist*(fpout.mu-fpout.mu[1:*]))
padist *= eflux/nrm
qe78 = dblarr(fpout.nz,n_elements(padist))
for i = 0,n_elements(padist)-1 do begin $
  if (padist[i] gt 0) then qe78[*,i] = qh_e78(ecut,dlt,padist[i],fpout.mum[i],atmfile)

; Integrate over the PA distribution
qe = total(qe78 * (replicate(1d0,191)#(fpout.mu-fpout.mu[1:*])),2) 

; Now plot the results
olddev = !d.name
oldfont = !p.font
if (doeps) then begin set_plot,'ps' & !p.font =0 & device,file='compare_fp_e78.eps',/enc,/color,/helv,/isolatin1, bits_per_pixel=8,xsize = 10.0,ysize = 8.0,/inches
loadct,39
if (~doeps) then device,decom = 0
red = 250 & blue = 64 & white = 255 & black =0
color0 = white
if (doeps) then color0 = black
pm = !p.multi
!p.multi = [0,1,2]
cs1 = 1.5
xm = [10,10]
if (doeps) then xm = [8,8]
thick = 3
if (doeps) then thick = 5

plot,fpout.zm/1d8,fpout.atm.tg,/ylog,yr = [1d3,1d7],thick = thick, charsize = cs1,ysty = 9, ytit = 'Temperature (K)', xmargin = xm, ymargin = [0,2], xtickname = replicate(' ',30),/xsty
axis,/yaxis,yr = [1d8,1d17],/ylog,/save,ytit = 'Density (cm!E-3!N)', charsize = cs1
oplot,fpout.zm/1d8,fpout.atm.dni[0,*],thick = thick, co=red
oplot,fpout.zm/1d8,fpout.atm.dnn[0,*],thick = thick, co=blue
legend,['Temperature','e!E-!N density', 'H I density'],color = [color0,red,blue], linesty = [0,0,0], thick = replicate(thick,3), charsize = cs1,pos = [.5d0,1d16]

plot,fpout.zm/1d8,fpout.heatrate,yr = [1d0,8d3], thick =thick, xtit = 'Distance along loop (Mm)', ytit = 'Heating Rate (erg cm!E-3!N s!E-1!N)', charsize = cs1, xmargin= xm, ymargin =[4,0], /ysty,/xsty,/ylog
oplot,fpout.zm/1d8,qe,co=red, thick =thick
oplot,fpout_all.zm/1d8,fpout_all.heatrate,thick =thick, co = 64
legend,['FP_CC_only','E78+HF94','FP'], color= [color0,red,blue], linesty = [0,0,0], thick = replicate(thick,3),charsize = cs1, pos = [.5d0,5000]

x = .50 & y = .46 & s = .04
xyouts,x,y,'!9d!x: '+strng(dlt,form = '(F5.1)'),/norm, charsize = cs1
y-=s
xyouts,x,y,'E!Dc!N: '+strng(ecut,form = '(F5.1)')+' keV',/norm,charsize = cs1
y-=s
xyouts,x,y,'F!D0!N: '+strng(eflux,form = '(E15.1)')+ ' erg cm!E-2!N s!E-1!N',/norm,charsize = cs1
;y-=s
;xyouts,x,y,'atm: '+atmfile, /norm,charsize = cs1
if (doeps) then device,/close
;RESET VARIABLES CHANGED IN PLOTTING
set_plot,olddev
!P.font = oldfont
!p.multi = pm
