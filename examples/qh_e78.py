def qh_e78(eckev, dlt, injflx, mu0, atmfile, proton = False):
# Implements Hawley & Fisher (1994) modification to Emslie (1978) thick target heating rate

  import numpy as np
  import scipy.special as sc
  import const as c
  import readatm
  ec = eckev * 1e3
  atm = readatm.readatm(atmfile)

  ea = (dlt-1)/(dlt-2)*ec /5# avg E
  ma = c.mp if proton else c.me

  bta = np.sqrt(1e0 - 1e0/np.power(ea / ma + 1,2))
  redm = c.me*ma / (c.me + ma)
  hbc2     = 6.1992097e-05
  sqrtpie2 = np.sqrt(np.pi/c.e2)
  alph = 0.0072973525693e0 #fine structure
  Kconst = 2*np.pi*c.e2*c.e2

  tothy = atm['dni'][:,1]+atm['dnn'][:,0] # total # hydrogen (neutral + protons)
  ndep = atm['zin'].size
  cl = np.zeros(ndep)
  for k in range(0, ndep):
    xi = c.me * ea / c.kb/atm['tg'][k]
    u = (1-1/np.sqrt(xi))*bta
    rmin = max([c.e2/redm/u/u, hbc2/u/redm])
    rmax = np.sqrt(c.me/atm['dni'][k,0])*bta*sqrtpie2
    cl[k] = max([np.log(rmax/rmin),0e0])

  clp = np.log(bta*(ea/ma+1)*np.sqrt(ea/ma) * c.me / atm['Enion'][0]) - .5*bta*bta
  if (proton): clp = np.log(2*pow(bta*(ea/ma+1),2)* c.me / atm['Enion'][0]) - bta*bta
  clpp = .5*np.log(pow((ea/ma+1)*bta,2)/pow(atm['Zn'][0],(2./3.))/2/alph/alph)
  xion = atm['dni'][:,1]/tothy

  gam = ma/c.me*(xion*cl + (1-xion)*clp)
  if (proton):
    bb = (xion-1)*clp/(clp + xion*(cl-clp))
  else:
    bb = (2*xion*cl + (1-xion)*clpp) / (clp + xion*(cl-clp))

  dz = np.insert(-np.diff(atm['zin']),0,0)
  N = np.cumsum(tothy*dz) # column depth
  NN = np.insert(np.diff(N),0,0)
  eqN = np.cumsum(gam / cl * NN) # equivalent depth (HF94 discussion of eq 2.10)
  eqNc = mu0*ec*ec/(2+.5*bb)/cl/Kconst
  b = sc.beta(dlt/2,2e0/(4e0 + bb))
  b *= np.where(eqN/eqNc < 1, sc.betainc(dlt/2,2e0/(4e0 + bb), eqN/eqNc),1)

  qh = .5 * Kconst/mu0 * gam * (dlt-2) * b*injflx/ec/ec*np.power(eqN/eqNc,-dlt/2) * tothy # HF94 eq 2.10
  qh[0] = qh[1]
  return(qh)
