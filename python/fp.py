class fpoutClass:
  pass

def readout(outfile,atmfile,fpout):
  import h5py
  import numpy as np
  with h5py.File(outfile, 'r') as h5f:
     fpout.inputparams = dict(h5f['inputparams'].attrs)
     fpout.atm = dict(h5f['atm'].attrs)
     fpout.atm['zin'] = np.array(h5f['atm/zin'])
     fpout.atm['tg'] = np.array(h5f['atm/tg'])
     fpout.atm['bfield'] = np.array(h5f['atm/bfield'])
     fpout.atm['dni'] = np.array(h5f['atm/dni'])
     fpout.atm['dnn'] = np.array(h5f['atm/dnn'])
     fpout.atm['mion'] = np.array(h5f['atm/mion'])
     fpout.atm['Zion'] = np.array(h5f['atm/Zion'])
     fpout.atm['Zn'] = np.array(h5f['atm/Zn'])
     fpout.atm['Enion'] = np.array(h5f['atm/Enion'])
     fpout.E = np.array(h5f['E'])
     fpout.mu = np.array(h5f['mu'])
     fpout.z = np.array(h5f['z'])
     fpout.esvol = np.array(h5f['esvol'])
     fpout.heatrate = np.array(h5f['heatrate'])
     fpout.momrate =  np.array(h5f['momrate'])
     fpout.f = np.array(h5f['f'])
  makefpoutstr(fpout,atmfile)

def makefpoutstr(fpout,atmfile):
  import numpy as np
  import const as c
  fpout.atmfile = atmfile
  fpout.nE = fpout.inputparams['nE']
  fpout.nmu = fpout.inputparams['nmu']
  fpout.nz = fpout.inputparams['nz']
  nE = fpout.nE
  nmu = fpout.nmu
  nz = fpout.nz
  fpout.Em = .5*(fpout.E[1:]+fpout.E[0:nE])
  mbeam = fpout.inputparams['mbeam']
  if (fpout.inputparams['inc_relativity']):
    fpout.gma = fpout.E/(mbeam/1e3) + 1e0
    fpout.bta = np.sqrt(1e0 - pow(fpout.gma,-2))
    fpout.gmam = fpout.Em/(mbeam/1e3) + 1e0
    fpout.btam = np.sqrt(1e0 - pow(fpout.gmam,-2))
  else:
    fpout.gma = np.ones(nE+1)
    fpout.bta = np.sqrt(2*fpout.E*1e3/mbeam)
    fpout.gmam = np.ones(nE)
    fpout.btam = np.sqrt(2*fpout.Em*1e3/mbeam)

  fpout.theta = np.arccos(fpout.mu)
  fpout.thetam = fpout.theta[0] if fpout.inputparams['oneD'] else .5*(fpout.theta[1:]+fpout.theta[0:nmu])
  fpout.mum = np.cos(fpout.thetam)

  fpout.zm = .5*(fpout.z[1:]+fpout.z[0:nz])
  #Calculate useful fluxes
  vz = c.clight *np.outer(fpout.mum, fpout.btam)
  vx = c.clight * np.outer(np.sin(fpout.thetam), fpout.btam)
  domega = fpout.esvol/(fpout.E[1:] - fpout.E[0:nE]) # solid angle part of esvol

  fpout.eflux = np.einsum('ijk,jk,k',fpout.f,vz*fpout.esvol,fpout.Em)*1e3*c.ergperev # Energy flux in z direction (erg cm^-2 s^-1)
  fpout.efluxx = np.einsum('ijk,jk,k',fpout.f,vx*fpout.esvol,fpout.Em)*1e3*c.ergperev # Energy flux in x direction (erg cm^-2 s^-1)
  fpout.nflux = np.einsum('ijk,jk',fpout.f,vz*fpout.esvol)# number flux in z direction (particles cm^-2 s^-1)
  fpout.flux = np.einsum('ijk,jk->ik',fpout.f,vz*domega) # number flux distribution in z direction (particles cm^-2 s^-1 keV^-1)
  fpout.flx = np.einsum('ijk,jk,k->ik',fpout.f,domega,fpout.btam)*c.clight # number flux speed distribution. (particles cm^-2 s^-1 keV^-1)

  if (abs(mbeam-c.me) < 1e3): #electron beam
    import fpbrem
    fpout.Eph = fpout.Em
    fpout.brem = fpbrem.fpbrem(fpout.Em,fpout.flx,fpout.Eph,fpout.atm)
    fpout.totbrem = np.einsum('ij,j',fpout.brem,np.diff(fpout.z))
  else:
  # Not electron beam so do not calculate bremsstrahlung
    fpout.brem = 0e0
    fpout.totbrem = 0e0
    fpout.Eph = 0e0

def writefpout(fpout,outfile = 'out.h5'):
  import h5py
  with h5py.File(outfile,'w') as fle:
    ipg = fle.create_group('inputparams')
    ipg.attrs.update(fpout.inputparams)
    atmattrs = {k:fpout.atm[k] for k in ('nIon','nNeutral') if k in fpout.atm}
    atmg = fle.create_group('atm')
    atmg.attrs.update(atmattrs)
    atmg.create_dataset('zin',data=fpout.atm['zin']); atmg.create_dataset('tg',data=fpout.atm['tg']); atmg.create_dataset('bfield',data=fpout.atm['bfield'])
    atmg.create_dataset('dni',data=fpout.atm['dni']); atmg.create_dataset('dnn',data=fpout.atm['dnn']); atmg.create_dataset('mion',data=fpout.atm['mion'])
    atmg.create_dataset('Zion',data=fpout.atm['Zion']); atmg.create_dataset('Zn',data=fpout.atm['Zn']); atmg.create_dataset('Enion',data=fpout.atm['Enion'])
    fle.create_dataset('E',data = fpout.E); fle.create_dataset('mu',data = fpout.mu); fle.create_dataset('z',data = fpout.z)
    fle.create_dataset('esvol',data = fpout.esvol); fle.create_dataset('heatrate',data = fpout.heatrate); fle.create_dataset('momrate',data = fpout.momrate)
    fle.create_dataset('f',data = fpout.f)

def writeparam(fle,outfile,nE,nmu,atmfile,inc_relativity,inc_cc, inc_synchro, inc_magmirror, inc_rc, oneD, reflecttop, reflectbottom, maxiter, tolres, toldiff,
               implicit_theta, mbeam, Zbeam, Ecut, dlt, Eflux, patype, pasigma, resist_fact, Emin, Emax, restart):

  truestr = ".true."
  falsestr = ".false."

  with open(fle,'w') as pfle:
    pfle.write("&control\n")
    pfle.write("nE = "+str(nE) + ",\n")
    pfle.write("nmu = "+str(nmu) + ",\n")
    pfle.write("Emin = "+  "%.17e" % Emin + ",\n")
    pfle.write("Emax = "+ "%.17e" % Emax + ",\n")
    pfle.write("inc_relativity = " + (truestr if inc_relativity else falsestr)  + ",\n")
    pfle.write("inc_CC = " + (truestr if inc_cc else falsestr) + ",\n")
    pfle.write("inc_synchro = " + (truestr if inc_synchro else falsestr) + ",\n")
    pfle.write("inc_magmirror = " + (truestr if inc_magmirror else falsestr) + ",\n")
    pfle.write("inc_RC = "+ (truestr if inc_rc else falsestr) + ",\n")
    pfle.write("oneD = "+ (truestr if oneD else falsestr) + ",\n")
    pfle.write("reflecttop = " +  (truestr if reflecttop else falsestr) + ",\n")
    pfle.write("reflectbottom = "+ (truestr if reflectbottom else falsestr) + ",\n")
    pfle.write("maxiter = "+str(maxiter)+",\n")
    pfle.write("writeoutput = .true.,\n")
    pfle.write("tolres = "+ "%.4e" % tolres +",\n")
    pfle.write("toldiff = " +"%.4e" % toldiff +",\n")
    pfle.write("implicit_theta = "+"%.17e" % implicit_theta+",\n")
    pfle.write("atmfile = '"+atmfile+"',\n")
    pfle.write("outfile = '"+outfile+"',\n")
    pfle.write("mbeam =  "+ "%.17e" % mbeam+",\n")
    pfle.write("Zbeam = "+"%.17e" % Zbeam+",\n")
    pfle.write("Ecut = "+"%.17e" % Ecut +",\n")
    pfle.write("dlt = "+"%.17e" % dlt+",\n")
    pfle.write("eflux= "+"%.17e" % Eflux+",\n")
    pfle.write("patype = "+str(patype)+",\n")
    pfle.write("pasigma = "+"%.17e" % pasigma+",\n")
    pfle.write("resist_fact = "+"%.17e" % resist_fact +",\n")
    pfle.write("restart = "+ (truestr if restart else falsestr) +"\n")
    pfle.write("/\n")

def solver(nE = 100, nmu=60, atmfile='atm.dat', inc_relativity = True, inc_cc = True, inc_synchro=True, inc_magmirror = False, inc_rc = True,
           oneD = False, reflecttop = False, reflectbottom = False, maxiter=100, tolres=1e-3, toldiff = 1e-4, implicit_theta = 1.0, mbeam = 0e0, Zbeam=-1.0,
           Ecut=20.0, dlt=5.0, Eflux=1e11, patype=2, pasigma=.05, resist_fact=1e0, Emin = 1e0, Emax = 0, restart= False, writeout = True, outfile = '',
           nthreads = 0, mpiexec = 'mpiexec'):
   import inspect
   import os
   import numpy as np
   import const as c
   import subprocess
   from random import random

   if (nthreads == 0):
     nthreads = os.cpu_count()
   if (mbeam == 0e0):
     mbeam = c.me
   if (Emax == 0):
     Emax = 3e3*Ecut

   iid = int(random()*1e6) # pick some random label to keep track of which param.cnt and out.dat go together
   paramfile = 'param.'+str(iid)+'.cnt'
   if (outfile == ''):
     outfile = 'out.'+str(iid)+'.h5'
     keepout = False
   else:
     outfile = str(outfile)
     keepout = True

   if (restart != False):
     if (isinstance(restart,fpoutClass)):
#       writefpout(restart,outfile=outfile)
       dorestart = True
       #Setup input parameters using values from restart
       nE = restart.inputparams['nE']; nmu = restart.inputparams['nmu']; atmfile = restart.atmfile; inc_relativity = restart.inputparams['inc_relativity']
       inc_CC = restart.inputparams['inc_CC']; inc_synchro = restart.inputparams['inc_synchro']; inc_magmirror = restart.inputparams['inc_magmirror']
       oneD = restart.inputparams['oneD']; reflecttop = restart.inputparams['reflecttop']; reflectbottom = restart.inputparams['reflectbottom']
       mbeam = restart.inputparams['mbeam']; Zbeam = restart.inputparams['zbeam']; Ecut = restart.inputparams['Ecut']; dlt = restart.inputparams['dlt']
       Eflux = restart.inputparams['eflux']; patype = restart.inputparams['patype']; pasigma = restart.inputparams['pasigma']
       resist_fact = restart.inputparams['resist_fact']; Emin = restart.inputparams['Emin']; Emax = restart.inputparams['Emax']
       if (maxiter is solver.__defaults__[0]): maxiter = restart.inputparams['maxiter']
       if (tolres == 1e-3): tolres = restart.inputparams['tolres']
       if (toldiff == 1e-4): todiff = restart.inputparams['toldiff']
       if (implicit_theta == 1.0): implicit_theta = restart.inputparams['implicit_theta']
       writefpout(restart,outfile)
     else:
       print('Cannot restart from ' + str(type(restart)))
       return
   else:
     dorestart = False

   import readatm
   atm = readatm.readatm(fle = atmfile)
   if (atm == None): return

   writeparam(paramfile,outfile,nE,nmu,atmfile,inc_relativity,inc_cc, inc_synchro, inc_magmirror, inc_rc, oneD, reflecttop, reflectbottom, maxiter, tolres, toldiff,
              implicit_theta, mbeam, Zbeam, Ecut, dlt, Eflux, patype, pasigma, resist_fact, Emin, Emax, dorestart)

   nz = atm['zin'].size
   nthreads = min(nthreads, int(nz/5))

   thisdir = os.path.dirname(inspect.getfile(writeparam))
   mpiexec = str(mpiexec)
   try:
     out = subprocess.run([mpiexec,'-n',str(nthreads),thisdir+'/fp',paramfile],stderr=subprocess.STDOUT)
   except OSError as err:
     print("OS error: {0}".format(err))
     if (os.path.exists(paramfile)): os.remove(paramfile)
     if (os.path.exists(outfile) and not keepout): os.remove(outfile)
     return None
   fpout = fpoutClass()
   readout(outfile,atmfile,fpout)
   os.remove(paramfile)
   if (not keepout): os.remove(outfile)
   return fpout

