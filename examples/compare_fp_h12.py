#RUNS FP in mode to mimic Holman 2012 RC  model:
#  1D mode with only return currect on
# Compares to 1.5D mode and including all forces.

# in IDL session run: @compare_fp_h12

##########################################################
#SET PATH TO INCLUDE FP/python
#TODO: Add check to determine if FP is already in path
import sys
sys.path.append('../python')

import numpy as np
import fp
import const as c

doeps = False
dohighres = False

# Injected beam parameters
dlt = 4e0 # delta
ecut = 20e0 # Cutoff energy in keV
eflux = 1e11 # injected energy flux (erg cm^-2 s^-1)

# Loop model is 13 Mm half length with a 3.4 MK corona.
atmfile = 'atm.13Mm.3MK.dat'

import readatm
atm = readatm.readatm(fle=atmfile)
tt = atm['tg'][0]
Eth = (dlt)*c.kb*tt/1e3

nE = 1000 if dohighres else 200
# Call FP
fpout1d = fp.solver(dlt = dlt, Ecut = ecut, Eflux=eflux, inc_rc = True, inc_cc=False, reflecttop= False, inc_magmirror = False, inc_synchro = False,
                    atmfile = atmfile, oneD = True, patype = 0, pasigma = .14, maxiter = 400, implicit_theta = 1e0, nE=nE,nmu = 1, toldiff = 1e-4,
                    tolres = 7e-3, Emin = .05, Emax = 2e3)

#better to make the 1.5D cases patype =2. patype = 0 creates very large (actually infinite) fluxes in the theta direction
fpout = fp.solver(dlt = dlt, Ecut = ecut, Eflux=eflux, inc_rc = True, inc_cc=False, reflecttop= False, inc_magmirror = False, inc_synchro = False,
                  atmfile = atmfile, oneD = False, patype = 2, pasigma = .14, maxiter = 1000, implicit_theta = .8e0, nE=nE,nmu = 60, toldiff = 5e-3,
                  tolres = 2e-3,Emin = .05, Emax = 2e3)

fpout_all = fp.solver(dlt = dlt, Ecut = ecut, Eflux=eflux, inc_rc = True, inc_cc=True, reflecttop= False, inc_magmirror = True, inc_synchro = True,
                      atmfile = atmfile, oneD = False, patype = 2, pasigma = .14, maxiter = 800, implicit_theta = .5e0, nE=nE,nmu = 60, toldiff = 1e-2,
                      tolres = 2e-3,Emin = .05, Emax = 2e3)

#The following can be used to store (and restore) the solutions to disk
#import shelve; sf = shelve.open('fpout.h12.db',flag='n'); sf['fpout1d'] = fpout1d; sf['fpout_all'] = fpout_all; sf['fpout'] = fpout; sf.close()
#import shelve; sf = shelve.open('fpout.h12.db',flag='r'); fpout1d = sf['fpout1d']; fpout_all = sf['fpout_all']; fpout = sf['fpout']; sf.close()

#load analytic model RCCTTM from Holman 2012
import qh_h12
qh12,ee,flux = qh_h12.qh_h12(ecut, dlt,eflux, Eth, atmfile)

##Make the plot
import matplotlib as mpl
import matplotlib.font_manager as font_manager
import matplotlib.pyplot as plt

#Setup font to mimic latex font
mpl.rcParams['font.family']='serif'
cmfont = font_manager.FontProperties(fname=mpl.get_data_path() + '/fonts/ttf/cmr10.ttf')
mpl.rcParams['font.serif']=cmfont.get_name()
mpl.rcParams['mathtext.fontset']='cm'
mpl.rcParams['axes.unicode_minus']=False

f, (ax1, ax2) = plt.subplots(2, 1)
ax1.set_ylabel('Heating rate (erg cm$^{-3}$ s$^{-1}$)')
ax1.set_xlabel('Distance along loop (Mm)')
l1, = ax1.semilogy(fpout1d.zm/1e8,fpout1d.heatrate,color = 'black',label = 'FP$\_$RC$\_$only$\_$1D')
l2, = ax1.semilogy(fpout.zm/1e8,fpout.heatrate, color = 'orange',label = 'FP$\_$RC$\_$only')
l3, = ax1.semilogy(fpout_all.zm/1e8,fpout_all.heatrate,color = 'blue',label = 'FP')
l4, = ax1.semilogy(fpout1d.zm/1e8,qh12, color = 'red',label = 'RCCTTM')
iz = 28
ax1.set_ylim(1e1,1e5)
ax1.axvline(fpout.zm[iz]/1e8,color ='black', linestyle='dashed',ymax = .4)
ax1.text(6.0,5e3,'$\delta$: '+str(dlt))
ax1.text(6.0,2e3,'E$_c$: '+str(int(ecut))+ ' keV')
ax1.text(6.0,8e2,'F$_0$: '+'{:7.1E}'.format(eflux) + ' erg cm$^{-2}$ s$^{-1}$')
ax1.legend(loc='upper left')

ax2.set_ylabel('Flux (10$^{17}$ e$^{-}$ cm$^{-2}$ s$^{-1}$ keV$^{-1}$)')
ax2.set_xlabel('Energy (keV)')
ax2.plot(fpout1d.Em,fpout1d.flux[iz,:]/1e17,color = 'black',label = 'FP$\_$RC$\_$only$\_$1D')
ax2.plot(fpout.Em,fpout.flux[iz,:]/1e17,color = 'orange',label = 'FP$\_$RC$\_$only')
ax2.plot(fpout_all.Em,fpout_all.flux[iz,:]/1e17,color = 'blue',label = 'FP')
ax2.plot(ee,flux[iz,:]/1e17,color = 'red',label = 'RCCTTM')
ax2.plot(ee,flux[0,:]/1e17,color = 'black',label = 'z = 0 (injected)',linestyle = 'dashed')
ax2.set_ylim(-1,4)
ax2.set_xlim(0,30)
ax2.legend()
f.tight_layout()
f.subplots_adjust(wspace = 0, hspace = .35,bottom = .1, top = .98, left =.1,right = .95)
f.show()
if (doeps): f.savefig('compare_fp_h12_py.eps',format = 'eps')
