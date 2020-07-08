def CoulogCC(mbeam,Zbeam, mi, Zi, ni, xi,b):
  import numpy as np
  import const as c
  sqrtpie2 = np.sqrt(np.pi/c.e2)
  hbc2     = 6.1992097e-05       # hbar c /2 in units of eV cm
  u = (1-1/np.sqrt(xi))*b
  redm = mi*mbeam/(mi + mbeam)
  rmin = np.maximum(c.e2 * abs(Zi * Zbeam)/redm/u/u, hbc2/u/redm)
  rmax = np.sqrt(mi/ni)*abs(Zi)*b*sqrtpie2
  cl = np.maximum(np.log(rmax/rmin),0e0)
  return(cl)

def CalcEtaIon(mi, Zi, ni, nel, T):
  import numpy as np
  import const as c
  xi = 1.5e0 * mi / c.me
  b = np.sqrt(1e0 - 1e0/(xi*c.kb*T/c.me +1)**2)
  cl = CoulogCC(c.me, -1e0, mi, Zi, ni, xi, b)
  FZ = (1e0 + 1.198e0*Zi + 0.222*Zi**2)/ (1e0 + 2.966e0*Zi + 0.753e0*Zi**2) #S.P. Hirshman, Phys. Fluids 20, 589 (1977)
  EtaIon = 4*np.sqrt(2*np.pi)/3*Zi*c.e2*np.sqrt(c.me)*cl*FZ/c.clight/(c.kb*T)**1.5 * ni / nel
  return( EtaIon)

def CalcEtaNeut(nn, nel, T):
  import numpy as np
  import const as c
  EtaNeut = 5.2e-11*3.759e-6*1e6 * c.me/c.e2/c.clight**2 * nn/nel * np.sqrt(T)# FROM MARTINEZ-SYKORA 2012 EQ 23-24
  return(EtaNeut)

def interpol(y,x,x0):
  import scipy.interpolate as interpolate
  f = interpolate.interp1d(x,y)
  return(f(x0))

def CalcEta(mion,Zion,dni,dnn,tg):
  import numpy as np
  import const as c
  nion = mion.size
  nneut =dnn[0,:].size

  nel = dni[:,0]
  eta = np.zeros(tg.size)
  for i in range(1, nion):
    ni = dni[:,i]
    eta += CalcEtaIon(mion[i],Zion[i],ni,nel,tg)
  nn = dnn[:,0]
  eta += CalcEtaNeut(nn,nel,tg)
  return(eta)

def qh_h12(Eckev, dlt, injflx, Ethkev, atmfile):
  # implements Holman 2012 return current heating rate (Eq. 24)
  import numpy as np
  import const as c
  import readatm

  Ec = Eckev * 1e3 # eV
  Eth = Ethkev * 1e3# eV

  Fe0 = injflx * (dlt-2)/(dlt-1)/(Ec * c.ergperev) # e- /cm^2/s
  #print('Fe0:',fe0)
  atm = readatm.readatm(fle=atmfile)
  eta = CalcEta(atm['mion'],atm['Zion'],atm['dni'],atm['dnn'],atm['tg'])

  rhoc = (Ec - Eth)/(c.e2*Fe0) # Eq. 19
  x = atm['zin'][0]-atm['zin']
  dx = np.insert(np.diff(x),0,0)
  rho = np.cumsum(eta*dx)
  xc = interpol(x,rho,rhoc)
  etac = interpol(eta,rho,rhoc)

  nx = x.size
  V= np.where(rho <= rhoc, c.e2*Fe0*rho, Ec*(dlt*c.e2*Fe0/Ec*(rho-rhoc) + 1)**(1/dlt) - Eth) #Eq 18
  #setup energy grid
  nE = 500
  ee = np.outer(np.ones(nx),np.geomspace(Eth/2,2e3*Ec,nE))
  VV = np.outer(V,np.ones(nE))
  Flux = np.where(ee >= Ec - VV, (dlt-1)*Fe0/Eckev * ((ee + VV)/Ec)**(-dlt), 0) #Eq 12
  ee = ee[0,:]/1e3 # convert to keV

  Q = np.where(rho < rhoc, eta * c.e2 * Fe0**2, eta * c.e2 * Fe0**2*(dlt*Eth/Ec + V/Ec)*(Eth/Ec+V/Ec)**(1-2*dlt)) #Eq. 24
  Q *= c.ergperev # erg /cm^3/s
  return(Q,ee,Flux)
