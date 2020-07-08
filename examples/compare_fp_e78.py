#RUNS FP in mode to mimic the cold thick target model:
#  Return current, magnetic mirroring and synchrotron forces are OFF

# in python session run: import compare_fp_e78

####################################################
doeps = False #SET THIS TO True TO MAKE AN EPS PLOT

# Injected beam parameters
dlt = 4e0 # delta
Ecut = 20e0 # Cutoff energy in keV
Eflux = 1e11 # injected energy flux (erg cm^-2 s^-1)

# Inject into a loop model that only includes H to most closely match the assumptions made in Emslie 1978
atmfile = 'atm.13Mm.3MK.onlyH.dat'
atmfile_all = 'atm.13Mm.3MK.dat'


###### DO THE RUN ##########
#SET PATH TO INCLUDE FP/python
#TODO: Add check to determine if FP is already in path
import sys
sys.path.append('../python')

import numpy as np
import fp

# Call FP with all terms switched OFF except Coulomb Collisions
fpout = fp.solver(dlt = dlt, Ecut = Ecut, Eflux=Eflux, inc_rc = False, reflecttop= True, inc_magmirror = False, inc_synchro = False, atmfile = atmfile, patype = 2,pasigma = .1, maxiter = 400, implicit_theta = .95e0, nE=100,nmu = 60, toldiff = 1e-3, tolres = 1e-3)
# Call FP with all terms switched ON
fpout_all = fp.solver(dlt = dlt, Ecut = Ecut, Eflux=Eflux, inc_rc = True, reflecttop= True, inc_magmirror = True, inc_synchro = True, atmfile = atmfile_all, patype = 2, pasigma = .1, maxiter = 400, implicit_theta = .9e0, nE=100,nmu = 60, toldiff = 1e-3, tolres = 1e-3)

#The following can be used to store (and restore) the solutions to disk
#import shelve; sf = shelve.open('fpout.e78.db'); sf['fpout'] = fpout; sf['fpout_all'] = fpout_all; sf.close()
#import shelve; sf = shelve.open('fpout.e78.db'); fpout = sf['fpout']; fpout_all = sf['fpout_all']; sf.close()

# Now calculate the heating rate from the thick target model using Emslie 1978
# Construct a Gaussian PA Distribution to match what was used in FP
padist = np.exp(-.5*np.power(fpout.thetam/fpout.inputparams['pasigma'],2))
padist = np.where(fpout.mum > 0, padist,0)
nrm = -np.sum(padist*np.diff(fpout.mu))
padist *= Eflux/nrm
qe78 = np.zeros([fpout.nmu,fpout.nz])
import qh_e78
for i in range(0,fpout.nmu):
  if (padist[i] > 0): qe78[i,:] = qh_e78.qh_e78(Ecut,dlt,padist[i],fpout.mum[i],atmfile)

# Integrate over the PA distribution
qe = -np.einsum('i,ij', np.diff(fpout.mu),qe78)

# Now plot the results
import matplotlib.pyplot as plt
import matplotlib as mpl
import matplotlib.font_manager as font_manager

mpl.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
mpl.rcParams['font.serif']=cmfont.get_name()
mpl.rcParams['mathtext.fontset']='cm'
mpl.rcParams['axes.unicode_minus']=False

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.set_ylabel('Temperature (K)')
l1, = ax1.semilogy(fpout.zm/1e8,fpout.atm['tg'],color = 'black',label = 'Temperature')
axd = ax1.twinx()
axd.set_ylabel('Density (cm$^{-3}$)')
l2, = axd.semilogy(fpout.zm/1e8,fpout.atm['dni'][:,0], color = 'red',label = 'e$^-$ density')
l3, = axd.semilogy(fpout.zm/1e8,fpout.atm['dnn'][:,0], color = 'blue',label = 'H I density')
axd.set_ylim(1e8,1e17)
axd.legend((l1,l2,l3),('Temperature', 'e$^-$ dens', ' H I dens'),loc='center left')

ax2.set_ylabel('Heating rate (erg cm$^{-3}$ s$^{-1}$)')
ax2.set_xlabel('Distance along loop (Mm)')
ax2.semilogy(fpout.zm/1e8,fpout.heatrate,label = 'FP$\_$CC$\_$only', color = 'black')
ax2.semilogy(fpout_all.zm/1e8,fpout_all.heatrate,label = 'FP',color = 'blue')
ax2.semilogy(fpout.zm/1e8,qe,label='E78+HF94',color = 'red')
ax2.set_ylim(1e0,8e3)
ax2.text(6.0,4e3,'$\delta$: '+str(dlt))
ax2.text(6.0,2e3,'E$_c$: '+str(int(Ecut))+ ' keV')
ax2.text(6.0,8e2,'F$_0$: '+'{:7.1E}'.format(Eflux) + ' erg cm$^{-2}$ s$^{-1}$')
ax2.legend()
f.tight_layout()
f.subplots_adjust(wspace = 0, hspace = .05,bottom = .1, top = .98, left =.1,right = .89)
f.show()
if (doeps): f.savefig('compare_fp_e78_py.eps',format = 'eps')
