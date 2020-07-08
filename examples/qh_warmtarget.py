def nbyw(x,mode,Ndepth,mu0,T,clog,ww):
   import scipy.special as sc
   import const
   import numpy as np

   if (mode == 0):
     w0 = x
     w = ww
   if (mode == 1):
     w = x
     w0 = ww

   if (w < 0): w = 0#complex(w)
   K = 2*np.pi*const.e2**2
   c = (2e0/3e0)*np.sqrt(const.me/const.mp/np.pi) # defined below Eq 33

   CapW = 1e0/(1e0+2*c*w**1.5)
   CapW0 = 1e0/(1e0+2*c*w0**1.5)
   thrd = 1e0/3e0
   thrd2 = 2e0/3e0
   bta = 2.0533902179391776#sc.beta(thrd2,thrd2)

   retval = mu0 * (const.kb*T)**2/c/clog/K *(np.sqrt(w0)*(1e0-w/w0*(CapW/CapW0)**thrd) - thrd2/CapW0**thrd/(2*c)**thrd2/np.sqrt(w0)*(sc.betainc(thrd2,thrd2,CapW)*bta- sc.betainc(thrd2,thrd2,CapW0)*bta)) -Ndepth
   retval = abs(retval)
   return(retval)

def qh_warmtarget(eckev, dlt, eferg, mu0, atmfile, proton = False):

  import scipy.optimize as sc
  import numpy as np
  import const
  import readatm
  eflux = eferg/const.ergperev
  ecut = eckev * 1e3
  atm = readatm.readatm(fle = atmfile)
  xi = atm['Enion'][0]
  ea = (dlt-1)/(dlt-2)*ecut /5# avg E
  ma = const.mp if proton else const.me

  nz = atm['zin'].size
  qh = np.zeros(nz)

  mu = mu0
  bta = np.sqrt(1e0 - 1e0/(ea / ma + 1)**2)
  redm = const.me*ma / (const.me + ma)
  hbc2     = 6.1992097e-05
  sqrtpie2 = np.sqrt(np.pi/const.e2)
  alph = 0.0072973525693e0 #fine structure
  Kconst = 2*np.pi*const.e2**2

  tothy = atm['dni'][:,1]+atm['dnn'][:,0] # total # hydrogen (neutral + protons)
  xion = atm['dni'][:,1]/tothy
  cl = np.zeros(nz)
  for k in range(0, nz):
    xi = const.me * ea / const.kb/atm['tg'][k]
    u = (1-1/np.sqrt(xi))*bta
    rmin = max([const.e2/redm/u/u, hbc2/u/redm])
    rmax = np.sqrt(const.me/atm['dni'][k,0])*bta*sqrtpie2
    cl[k] = max([np.log(rmax/rmin),0e0])

  clp = np.log(bta*(ea/ma+1)*np.sqrt(ea/ma) * const.me / atm['Enion'][0]) - .5*bta**2
  if (proton): clp = np.log(2*bta**2*(ea/ma+1)**2* const.me / atm['Enion'][0]) - bta**2
  clpp = .5*np.log(((ea/ma+1)*bta)**2/atm['Zn'][0]**(2./3.)/2/alph**2)
  dz = np.insert(-np.diff(atm['zin']),0,0)
  N = np.cumsum(tothy*dz) # column depth

  c = (2e0/3e0)*np.sqrt(const.me/const.mp/np.pi) # defined below Eq 33
  for k in range(0,nz):
    # numerically solve Tamres et al Eq. 38
    Ndepth = N[k]
    T = atm['tg'][k]
    wc = ecut / const.kb/ T
    clog = cl[k]
    nw = 100
    mode=0; w = 1e0
    wm = (Ndepth * Kconst * clog * c /mu0/(const.kb*T)**2 + 4.8838)**2
    args = (mode,Ndepth,mu0,T,clog,w)
    wm = sc.fsolve(nbyw,wm,args = args)
    wm = wm[0]
    if (wm < wc): wm = wc
    wh = 1e4*wm
    w0 = np.geomspace(wm,wh,nw)
    mode =1
    integrand = np.zeros(nw)
    for i in range(0,nw):
      # numerically solve Eq. 35 to get omega from N and omega0
      w = w0[i]/2
      args = (mode,Ndepth,mu0,T,clog,w0[i])
      w = sc.fsolve(nbyw,w,args = args)
      w = max(w,1e0) # set w to minimum of 1.0
      Phi0 = eflux*(dlt-2)/ecut**2*(w0[i]/wc)**(-dlt)/(1e0-(wh/wc)**(2-dlt))
      integrand[i] = np.sqrt(w0[i])*(1+2*c*w**1.5)**(4e0/3e0)*Phi0/mu0/w**1.5/(1+2*c*w0[i]**1.5)**(1e0/3e0)

    qh[k] = Kconst*clog*np.trapz(integrand,x=w0)*tothy[k]*const.ergperev

  return(qh)
