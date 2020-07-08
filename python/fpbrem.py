#+
# PROJECT:
#	HESSI
#
# NAME:
#	Brm_BremCross
#
# PURPOSE:
#	Computes the relativistic cross section for electron-ion bremsstrahlung,
#	differential in photon energy.
#
# CATEGORY:
#	HESSI, Spectra, Modeling
#
# CALLING SEQUENCE:
#	Brm_BremCross, eel, eph, z, cross
#
# CALLS:
#	none
#
# INPUTS:
#
#	eel		-	Array of abscissas at which to compute the bremsstrahlung cross
#				section.
#
#	eph		-	Array of photon energies corresponding to the array eel.
#
#	z		-	Mean atomic number of the target plasma.
#
#
# OPTIONAL INPUTS:
#	none
#
# OUTPUTS:
#	cross	-	Array of bremsstrahlung cross sections corresponding to the
#				input array eel.
#
# OPTIONAL OUTPUTS:
#	none
#
# KEYWORDS:
#	none
#
# COMMON BLOCKS:
#	none
#
# SIDE EFFECTS:
#
#
# RESTRICTIONS:
#
#
# PROCEDURE:
#	The cross section is from Equation (4) of E. Haug (Astron. Astrophys. 326,
#	417, 1997).  This closely follows Formula 3BN of H. W. Koch & J. W. Motz
#	(Rev. Mod. Phys. 31, 920, 1959), but requires fewer computational steps.
#	The multiplicative factor introduced by G. Elwert (Ann. Physik 34, 178,
#	1939) is included.
#
# MODIFICATION HISTORY:
#   Version 1, holman@stars.gsfc.nasa.gov, 23 July 2001
#   IDL version:  Sally House, summer intern
#	Documentation for the Fortran version of this code can be found at
#   http://hesperia.gsfc.nasa.gov/hessi/modelware.htm
#   03/11/03    ch2 corrected (C12 replaced by C21) found by G. Emslie
#
#-

def Brm_BremCross( eel, eph, z):
   import numpy as np
   #	Physical coefficients.

   mc2 = 510.98e+00
   alpha = 7.29735308e-03
   twoar02 = 1.15893512e-27

#	Numerical coefficients.

   c11 = 4.e0/3.e0
   c12 = 7.e0/15.e0
   c13 = 11.e0/70.e0
   c21 = 7.e0/20.e0
   c22 = 9.e0/28.e0
   c23 = 263.e0/210.e0

#	Calculate normalized photon and total electron energies.

   k = eph/mc2
   e1 = (eel/mc2)+1.e+00

#	Calculate energies of scatter electrons and normalized momenta.

   e2 = e1-k
   p1 = np.sqrt(e1**2-1.e+00)
   p2 = np.sqrt(e2**2-1.e+00)

#	Define frequently used quantities.

   e1e2 = e1*e2
   p1p2 = p1*p2
   p2sum = p1**2+p2**2
   k2 = k**2
   e1e23 = e1e2**3
   pe = p2sum/e1e23

#	Define terms in cross section.

   ch1 = (c11*e1e2+k2)-(c12*k2/e1e2)-(c13*k2*pe/e1e2)

   ch2 = 1.e0+(1.e0/e1e2)+(c21*pe)+(c22*k2+c23*p1p2**2)/e1e23

#	Collect terms.

   crtmp = ch1*(2.e0*np.log((e1e2+p1p2-1.e0)/k)-(p1p2/e1e2)*ch2)

   crtmp = z**2*crtmp/(k*p1**2)

#	Compute the Elwert factor.

   a1 = alpha*z*e1/p1
   a2 = alpha*z*e2/p2

   fe = (a2/a1)*(1.e0-np.exp(-2.e0*np.pi*a1))/(1.e0-np.exp(-2.e0*np.pi*a2))

#	Compute the differential cross section (units cm**2).

   cross = twoar02*fe*crtmp
   return(cross)

# Computes X-ray bremsstrahlung from the electron flux distribution
def fpbrem(e, flx, eph, atm, Zan = 1.2e0):
  import numpy as np
  import const
  mc2 = const.me/1e3
  # returns photons/cm^3/s/kev
  #flx is the direction independent flux so flx = integrate(f*beta*c*dmu). units of flx are electrons/cm^2/s/keV
  # photon spectrum (# of photons / cm^2 /s/ keV) = n_ambient * l * integral(dsigma/dk(E,k) * electron flux(E) * dE)
  ne = e.size
  neph = eph.size

  #Zan = Average atomic number

  ndep = atm['zin'].size
  xraybrem = np.zeros((neph, ndep))
  #dsigmadk = np.zeros((ne,neph))
  dsigmadk = np.zeros((neph,ne))

  for j in range(0, ne):
    idx = np.where( eph < e[j])
    if (idx[0].size < 1): continue
    cross = Brm_BremCross(e[j], eph[idx], Zan)
    dsigmadk[idx,j] = cross /mc2 # units of photons cm^2 / keV / incoming electron

  nh = atm['dnn'][:,0] + atm['dni'][:,1] # hydrogen number density
  #calculate spectrum at each loop grid cell
#  now integrate over e
   #photons /cm^2 /cm/s/keV
  for i in range(0,ndep):
    xraybrem[:,i] =  np.trapz( dsigmadk * flx[i,:], x = e,axis = 1) * nh[i]

  return(xraybrem)
