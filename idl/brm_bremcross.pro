;+
; PROJECT:
;	HESSI
;
; NAME:
;	Brm_BremCross
;
; PURPOSE:
;	Computes the relativistic cross section for electron-ion bremsstrahlung,
;	differential in photon energy.
;
; CATEGORY:
;	HESSI, Spectra, Modeling
;
; CALLING SEQUENCE:
;	Brm_BremCross, eel, eph, z, cross
;
; CALLS:
;	none
;
; INPUTS:
;
;	eel		-	Array of abscissas at which to compute the bremsstrahlung cross
;				section.
;
;	eph		-	Array of photon energies corresponding to the array eel.
;
;	z		-	Mean atomic number of the target plasma.
;
;
; OPTIONAL INPUTS:
;	none
;
; OUTPUTS:
;	cross	-	Array of bremsstrahlung cross sections corresponding to the
;				input array eel.
;
; OPTIONAL OUTPUTS:
;	none
;
; KEYWORDS:
;	none
;
; COMMON BLOCKS:
;	none
;
; SIDE EFFECTS:
;
;
; RESTRICTIONS:
;
;
; PROCEDURE:
;	The cross section is from Equation (4) of E. Haug (Astron. Astrophys. 326,
;	417, 1997).  This closely follows Formula 3BN of H. W. Koch & J. W. Motz
;	(Rev. Mod. Phys. 31, 920, 1959), but requires fewer computational steps.
;	The multiplicative factor introduced by G. Elwert (Ann. Physik 34, 178,
;	1939) is included.
;
; MODIFICATION HISTORY:
;   Version 1, holman@stars.gsfc.nasa.gov, 23 July 2001
;   IDL version:  Sally House, summer intern
;	Documentation for the Fortran version of this code can be found at
;   http://hesperia.gsfc.nasa.gov/hessi/modelware.htm
;   03/11/03    ch2 corrected (C12 replaced by C21) found by G. Emslie
;
;-

Pro Brm_BremCross, eel, eph, z, cross

;	Physical coefficients.

mc2 = 510.98d+00
alpha = 7.29735308d-03
twoar02 = 1.15893512d-27

;	Numerical coefficients.

c11 = 4.d0/3.d0
c12 = 7.d0/15.d0
c13 = 11.d0/70.d0
c21 = 7.d0/20.d0
c22 = 9.d0/28.d0
c23 = 263.d0/210.d0

;	Calculate normalized photon and total electron energies.

k = eph/mc2
e1 = (eel/mc2)+1.d+00

;	Calculate energies of scatter electrons and normalized momenta.

e2 = e1-k
p1 = SqRt(e1^2-1.d+00)
p2 = SqRt(e2^2-1.d+00)

;	Define frequently used quantities.

e1e2 = e1*e2
p1p2 = p1*p2
p2sum = p1^2+p2^2
k2 = k^2
e1e23 = e1e2^3
pe = p2sum/e1e23

;	Define terms in cross section.

ch1 = (c11*e1e2+k2)-(c12*k2/e1e2)-(c13*k2*pe/e1e2)

ch2 = 1.d0+(1.d0/e1e2)+(c21*pe)+(c22*k2+c23*p1p2^2)/e1e23

;	Collect terms.

crtmp = ch1*(2.d0*ALog((e1e2+p1p2-1.d0)/k)-(p1p2/e1e2)*ch2)

crtmp = z^2*crtmp/(k*p1^2)

;	Compute the Elwert factor.

a1 = alpha*z*e1/p1
a2 = alpha*z*e2/p2

fe = (a2/a1)*(1.d0-Exp(-2.d0*!dpi*a1))/(1.d0-Exp(-2.d0*!dpi*a2))

;	Compute the differential cross section (units cm^2).

cross = twoar02*fe*crtmp

END