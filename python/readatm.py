def readatm(fle='atm.dat'):
   from scipy.io import FortranFile
   import numpy as np
   try:
     fl = FortranFile(fle, 'r')
     sz = fl.read_ints(np.int32)
     ar = fl.read_reals()
     fl.close()
   except OSError as err:
     print("OS error: {0}".format(err))
     return None
   nz = sz[0]
   nions = sz[1]
   nneutrals = sz[2]
   atm = {'zin':ar[0:nz]}
   atm['tg'] = ar[nz:2*nz]
   atm['bfield'] = ar[2*nz:3*nz]
   atm['dni'] = ar[3*nz:(3+nions)*nz].reshape((nz,nions))
   atm['dnn'] = ar[(3+nions)*nz:(3+nions+nneutrals)*nz].reshape((nz,nneutrals))
   atm['mion'] = ar[(3+nions+nneutrals)*nz:(3+nions+nneutrals)*nz+nions]
   atm['Zion'] = ar[(3+nions+nneutrals)*nz+nions:(3+nions+nneutrals)*nz+2*nions]
   atm['Zn'] = ar[(3+nions+nneutrals)*nz+2*nions:(3+nions+nneutrals)*nz+2*nions+nneutrals]
   atm['Enion'] = ar[(3+nions+nneutrals)*nz+2*nions+nneutrals:(3+nions+nneutrals)*nz+2*nions+2*nneutrals]
   atm['z'] = atm['zin'][0]-atm['zin']
   return(atm)
