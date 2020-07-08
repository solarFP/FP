#RUNS FP for proton beam injection
# Compares to warm target of Tamres et al 1986 and cold target of Emslie 1978
#
#  Return current, magnetic mirroring and synchrotron forces are OFF

# in python session run: import compare_fp_t86

##########################################################
doeps = False #SET THIS TO True TO MAKE AN EPS PLOT

###### DO THE RUN ##########
#SET PATH TO INCLUDE FP/python
#TODO: Add check to determine if FP is already in path
import sys
sys.path.append('../python')

import numpy as np
import fp
import const as c

# Ec = 100 keV
eflux = 1e11
dlt = 4e0
Ecut1 = 1e2
atm = 'atm.13Mm.3MK.onlyH.dat'
fpout1 = fp.solver(Eflux = eflux, dlt = dlt, Ecut = Ecut1, inc_rc = False, inc_magmirror = False, inc_synchro = False, mbeam = c.mp, Zbeam = 1e0, atmfile = atm, maxiter = 200, implicit_theta =.95, nE = 300, nmu = 100, patype = 0, Emin = 1e1, Emax = 1e6)
fpout1_all = fp.solver(Eflux = eflux, dlt = dlt, Ecut = Ecut1, inc_rc = True, inc_magmirror = True, inc_synchro = True, mbeam = c.mp, Zbeam = 1e0, atmfile = atm, maxiter = 200, implicit_theta=.95, nE = 300, nmu = 100, patype = 0, Emin = 1e1, Emax = 1e6)
import qh_warmtarget
qhwt1 = qh_warmtarget.qh_warmtarget(Ecut1,dlt,eflux,1e0, atm, proton = True)
import qh_e78
qhct1 = qh_e78.qh_e78(Ecut1,dlt,eflux,1e0, atm, proton = True)

# Ec = 1 MeV
eflux = 1e11
dlt = 4e0
Ecut2 = 1e3
atm = 'atm.13Mm.3MK.onlyH.dat'
fpout2 = fp.solver(Eflux = eflux, dlt = dlt, Ecut = Ecut2, inc_rc = False, inc_magmirror = False, inc_synchro = False, mbeam = c.mp, Zbeam = 1e0, atmfile = atm, maxiter = 200, implicit_theta=.95, nE = 300, nmu = 100, patype = 0, Emin = 1e1, Emax = 1e6)
fpout2_all = fp.solver(Eflux = eflux, dlt = dlt, Ecut = Ecut2, inc_rc  = True, inc_magmirror = True, inc_synchro = True, mbeam = c.mp, Zbeam = 1e0, atmfile = atm, maxiter = 200, implicit_theta=.95, nE = 300, nmu = 100, patype = 0, Emin = 1e1, Emax = 1e6)
qhwt2 = qh_warmtarget.qh_warmtarget(Ecut2,dlt,eflux,1e0, atm, proton = True)
qhct2 = qh_e78.qh_e78(Ecut2,dlt,eflux,1e0, atm, proton = True)

#The following can be used to store (and restore) the solutions to disk
#OLDER VERSIONS OF PYTHON DBM CREATE A 32GB FILE TO STORE THIS SHELF. CRAZY!! IF THIS HAPPENS UPDATE TO dbm.gnu
#import shelve; sf = shelve.open('fpout.t86.db',flag='n'); sf['fpout1'] = fpout1; sf['fpout1_all'] = fpout1_all; sf['fpout2'] = fpout2; sf['fpout2_all'] = fpout2_all; sf['qhct1'] = qhct1; sf['qhct2'] = qhct2; sf['qhwt1'] = qhwt1; sf['qhwt2'] = qhwt2; sf.close()
#import shelve; sf = shelve.open('fpout.t86.db',flag='r'); fpout1 = sf['fpout1']; fpout1_all = sf['fpout1_all']; fpout2 = sf['fpout2']; fpout2_all = sf['fpout2_all']; qhct1 = sf['qhct1']; qhct2 = sf['qhct2']; qhwt1 = sf['qhwt1']; qhwt2 = sf['qhwt2']; sf.close()

#Plot the results
import matplotlib as mpl
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt

mpl.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
mpl.rcParams['font.serif']=cmfont.get_name()
mpl.rcParams['mathtext.fontset']='cm'
mpl.rcParams['axes.unicode_minus']=False

f, (ax1, ax2) = plt.subplots(2, 1, sharex=True)
ax1.set_ylabel('Heating rate (erg cm$^{-3}$ s$^{-1}$)')
l1, = ax1.semilogy(fpout1.zm/1e8,fpout1.heatrate,color = 'black',label = 'FP$\_$CC$\_$only')
l2, = ax1.semilogy(fpout1_all.zm/1e8,fpout1_all.heatrate,color = 'blue',label = 'FP')
l3, = ax1.semilogy(fpout1.zm/1e8,qhct1, color = 'red',label = 'E78')
l4, = ax1.semilogy(fpout1.zm/1e8,qhwt1, color = 'orange',label = 'T86')
ax1.set_ylim(1e0,1e4)
ax1.text(6.0,5e3,'$\delta$: '+str(dlt))
ax1.text(6.0,2e3,'E$_c$: '+str(int(Ecut1))+ ' keV')
ax1.text(6.0,8e2,'F$_0$: '+'{:7.1E}'.format(eflux) + ' erg cm$^{-2}$ s$^{-1}$')
ax1.legend(loc='upper left')

ax2.set_ylabel('Heating rate (erg cm$^{-3}$ s$^{-1}$)')
ax2.set_xlabel('Distance along loop (Mm)')
ax2.semilogy(fpout2.zm/1e8,fpout2.heatrate,label = 'FP$\_$CC$\_$only', color = 'black')
ax2.semilogy(fpout2_all.zm/1e8,fpout2_all.heatrate,label = 'FP',color = 'blue')
ax2.semilogy(fpout2.zm/1e8,qhct2,label='E78',color = 'red')
ax2.semilogy(fpout2.zm/1e8,qhwt2,label='T86',color = 'orange')
ax2.set_ylim(1e0,1e4)
ax2.text(6.0,5e3,'$\delta$: '+str(dlt))
ax2.text(6.0,2e3,'E$_c$: '+str(int(Ecut2))+ ' keV')
ax2.text(6.0,8e2,'F$_0$: '+'{:7.1E}'.format(eflux) + ' erg cm$^{-2}$ s$^{-1}$')
ax2.legend()
f.tight_layout()
f.subplots_adjust(wspace = 0, hspace = .05,bottom = .1, top = .98, left =.1,right = .89)
f.show()
if (doeps): f.savefig('compare_fp_t86_py.eps',format = 'eps')
