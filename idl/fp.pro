PRO writeparam,file,outfile, nen,nmu,atmfile,inc_relativity, inc_cc,inc_synchro,inc_magmirror,inc_rc, oneD, reflecttop, reflectbottom, maxiter, $ 
               tolres, toldiff, implicit_theta, mbeam, zbeam, ecut,dlt,eflux,patype, pasigma,resist_fact, Emin,Emax, restart

  truestr = ".true." & falsestr = ".false."
  openw, lu,file, /get_lun
  printf,lu,"&control"
  printf,lu,"nE = "+strng(nen) + ","
  printf,lu,"nmu = "+strng(nmu) + ","
  printf,lu,"Emin = "+strng(Emin,form = '(E26.17)') + ","
  printf,lu,"Emax = "+strng(Emax,form = '(E26.17)') + ","
  printf,lu,"inc_relativity = " + (inc_relativity ? truestr : falsestr) + ","
  printf,lu,"inc_CC = " + (inc_cc ? truestr : falsestr) + ","
  printf,lu,"inc_synchro = " + (inc_synchro ? truestr : falsestr) + ","
  printf,lu,"inc_magmirror = " + (inc_magmirror ? truestr : falsestr) + ","
  printf,lu,"inc_RC = "+ (inc_rc ? truestr : falsestr) + ","
  printf,lu,"oneD = "+ (oneD ? truestr : falsestr) + ","
  printf,lu,"reflecttop = " + (reflecttop ? truestr : falsestr) + ","
  printf,lu,"reflectbottom = "+ (reflectbottom ? truestr : falsestr) + ","
  printf,lu,"maxiter = "+strng(maxiter)+","
  printf,lu,"writeoutput = .true.,"
  printf,lu,"tolres = "+strng(tolres,format = '(E10.4)')+","
  printf,lu,"toldiff = "+strng(toldiff, format = '(E10.4)')+","
  printf,lu,"implicit_theta = "+strng(implicit_theta,form='(E26.17)')+","
  printf,lu,"atmfile = '"+atmfile+"',"
  printf,lu,"outfile = '"+outfile+"',"
  printf,lu,"mbeam =  "+strng(mbeam, form = '(E26.17)')+","
  printf,lu,"zbeam = "+strng(zbeam, form = '(E26.17)')+","
  printf,lu,"Ecut = "+strng(ecut,form = '(E26.17)')+","
  printf,lu,"dlt = "+strng(dlt,form = '(E26.17)')+","
  printf,lu,"eflux= "+strng(eflux,form = '(E26.17)')+","
  printf,lu,"patype = "+strng(patype)+","
  printf,lu,"pasigma = "+strng(pasigma,form = '(E26.17)')+","
  printf,lu,"resist_fact = "+strng(resist_fact,form = '(E26.17)')+","
  printf,lu,"restart = "+ (restart ? truestr : falsestr)
  printf,lu,"/"
  free_lun,lu
END

PRO fp, nEn = nEn, nmu = nmu, atmfile = atmfile, inc_cc = inc_cc, inc_synchro = inc_synchro, inc_magmirror = inc_magmirror, inc_rc = inc_rc, $
        oneD = oneD, reflecttop = reflecttop, reflectbottom = reflectbottom, maxiter = maxiter, tolres = tolres, toldiff = toldiff, $ 
        implicit_theta = implicit_theta, mbeam = mbeam, zbeam = zbeam, Ecut = Ecut, dlt = dlt, eflux = eflux, patype = patype, $ 
        pasigma = pasigma, resist_fact = resist_fact, restart = restart, writeout = writeout, Emin = Emin, Emax = Emax, $
        nthreads = nthreads, inc_relativity = inc_relativity, outfile = outfile, mpiexec = mpiexec, fpout
; IDL routine to call fortran Fokker-Planck solver

@const
on_error,2
if (n_elements(restart) eq 0) then restart = 0L else restart = long(restart[0])
if (restart) then begin
  if ( (size(fpout,/type) ne 8) ) then begin
    print,'Input structure is not defined. Cannot restart.'
    return
  endif
  dorestart = 1
endif else dorestart = 0

if (dorestart) then resv = fpout.inputparams.nEn & setfpparam,nEn,100L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.nmu & setfpparam,nmu,60L,resv,dorestart
if (dorestart) then resv = fpout.atmfile & setfpparam,atmfile,'atm.dat',resv,dorestart
if (dorestart) then resv = fpout.inputparams.inc_relativity & setfpparam,inc_relativity,1L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.inc_cc & setfpparam,inc_cc,1L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.inc_synchro & setfpparam,inc_synchro,1L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.inc_magmirror & setfpparam,inc_magmirror,0L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.inc_rc & setfpparam,inc_rc,1L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.oneD & setfpparam,oneD,0L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.reflecttop & setfpparam,reflecttop,0L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.reflectbottom & setfpparam,reflectbottom,0L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.maxiter & setfpparam,maxiter,100L,resv,dorestart,/use_cmd
if (dorestart) then resv = fpout.inputparams.tolres & setfpparam,tolres,1d-3,resv,dorestart,/use_cmd
if (dorestart) then resv = fpout.inputparams.toldiff & setfpparam,toldiff,1d-4,resv,dorestart,/use_cmd
if (dorestart) then resv = fpout.inputparams.implicit_theta & setfpparam,implicit_theta,1d0,resv,dorestart,/use_cmd
if (dorestart) then resv = fpout.inputparams.mbeam & setfpparam,mbeam,!me,resv,dorestart
if (dorestart) then resv = fpout.inputparams.zbeam & setfpparam,zbeam,-1d0,resv,dorestart
if (dorestart) then resv = fpout.inputparams.Ecut & setfpparam,Ecut,20d0,resv,dorestart
if (dorestart) then resv = fpout.inputparams.dlt & setfpparam,dlt,5d0,resv,dorestart
if (dorestart) then resv = fpout.inputparams.eflux & setfpparam,eflux,1d11,resv,dorestart
if (dorestart) then resv = fpout.inputparams.patype & setfpparam,patype,2L,resv,dorestart
if (dorestart) then resv = fpout.inputparams.pasigma & setfpparam,pasigma,0.05d0,resv,dorestart
if (dorestart) then resv = fpout.inputparams.resist_fact & setfpparam,resist_fact,1d0,resv,dorestart
if (dorestart) then resv = fpout.inputparams.Emin & setfpparam,Emin,1d0,resv,dorestart
if (dorestart) then resv = fpout.inputparams.Emax & setfpparam,Emax,3d3*Ecut,resv,dorestart
if (n_elements(writeout) eq 0) then writeout = 1L else writeout = long(writeout[0])
if (n_elements(nthreads) eq 0) then nthreads = long(!cpu.tpool_nthreads) else nthreads = long(nthreads)
if (n_elements(mpiexec) eq 0) then mpiexec = 'mpiexec' else mpiexec = string(mpiexec)

id = long(randomu(seed)*1d6) ; pick some random label to keep track of which param.cnt and out.dat go together
paramfile = 'param.'+strng(id)+'.cnt'
if (n_elements(outfile) eq 0) then begin 
  outfile = 'out.'+strng(id)+'.h5' 
  keepout = 0
endif else begin
  outfile = strng(outfile[0])
  keepout = 1
endelse

if (dorestart) then writefpout,fpout,outfile=outfile

writeparam,paramfile, outfile, nen,nmu,atmfile,inc_relativity, inc_cc, inc_synchro, inc_magmirror, inc_rc, oneD, reflecttop, reflectbottom, $
           maxiter, tolres, toldiff, implicit_theta, mbeam, zbeam, ecut,dlt,eflux,patype, pasigma,resist_fact, Emin,Emax, restart 

; Is there a better way to get nz?
readatm,atmfile
nz = n_elements(zin)
nthreads = nthreads  < nz/5L

dir = routine_filepath('fp')
dir = strmid(dir,0,strpos(dir,'fp.pro')-1)
cmd = mpiexec + " -n "+strng(nthreads) + " "+dir+"/fp "+paramfile
spawn,cmd, exit_status=exit_status
if (exit_status ne 0) then begin
  if (file_test(paramfile)) then file_delete,paramfile
  if (~keepout and file_test(outfile)) then file_delete,outfile
  return
endif
readfpout,fpout,atmfile,file = outfile
file_delete,paramfile 
if (~keepout) then file_delete,outfile

END
