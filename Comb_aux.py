## Noah Wilson
## 12/10/2020
## 4/8/2021 -- Updated a few lines to work with Python 3.
## Auxiliary script meant to be imported by Comb_sim.py

import numpy as np

def get_as(spec):
    # Extracts coefficients for NASA 7-polynomial equations from therm.dat
    # Values in therm.dat were taken from the Extended Third Millennium
    # Thermodynamic Database for Combustion with updates from  Active
    # Thermochemical Tables by Elke Goos, Alexander Burcat and Branko Ruscic.
    # This is a studnet project, database is NOT being used for commercial purposes.

    with open('therm.dat','r') as f:
        Lines = f.readlines()
        lines = [line for line in Lines]
        del(Lines)
        for i in range(len(lines)):
            if '>'+spec in lines[i]: # Finds the given species name in the data
                datalines = lines[i+1:i+4] # Extracts only lines with the values
                break
            else:
                pass
        # The next two lines parses the string data to get floating point values
        onestr = ''.join([data[0:-6] for data in datalines])
        Data = np.array([float(onestr[i:i+15]) for i in range(0,len(onestr),15)])
        # Values are normalized with R in cal/molK,
        # Renorm renormalizes and converts to J/molK
        Renorm = 8.314462582
        as_highT =  Renorm * Data[0:7]
        as_lowT = Renorm * Data[7:-1]
        hf0 = Renorm * Data[-1] # Enthalpy of formation at 298 K
        return as_highT,as_lowT,hf0

def Cp_T(spec,T):
    #Calculates heat capacity in J/mol*K at temperature T for a given species
    #Uses NASA 7-polynomial equations
    lowhigh = {True:spec.As_LowT[0:5],False:spec.As_HighT[0:5]}

    a1,a2,a3,a4,a5 = [a for a in lowhigh[T<=1000]]

    Cpout = a1 + a2*T + a3*T**2 + a4*T**3 + a5*T**4

    return Cpout

def h_T(spec,T):
    # Calculates enthalpy in J/mol at a temperature T for a given species
    #Uses NASA 7-polynomial equations
    lowhigh = {True:spec.As_LowT[0:6],False:spec.As_HighT[0:6]}

    a1,a2,a3,a4,a5,a6 = [a for a in lowhigh[T<=1000]]

    hout = T*(a1 + .5*a2*T + (a3/3)*T**2 + .25*a4*T**3 + .2*a5*T**4 + a6/T)

    return hout

def k_T(spec,speclist,state):
    # Thermal conductivity assuming Le = 1 and using mixture-averaged diffusion
    # coefficients
    i = spec.idx
    N = range(len(speclist))
    D = Dim(state.T,state.P,state.MWmix,speclist)
    return state.rho*D[i]*spec.Cp # Updated thermal conductivity for species i

def totk_T(speclist,MWmix):
    Xs = [(MWmix*sp.Y)/sp.MW for sp in speclist]
    Is = [sp.idx for sp in speclist]
    ks = [sp.k for sp in speclist]

    term1 = sum([Xs[i]*ks[i] for i in Is])
    term2 = 1/sum([Xs[i]/ks[i] for i in Is])
    return .5*(term1+term2)

def Dim(T,P,MWmix,speclist):

    # This calculates the mixture-averaged diffusion coefficients in m^2/s
    # T is Temperature in K
    # P is Pressure in Pa
    # speclist is a list of chemical species class instances from comb_sim.py
    N = len(speclist) # Number of species
    Kdelt = np.identity(N) # Identity matrix used as a Kronecker Delta
    MW = np.array([[sp.MW for sp in speclist]]*N)*1000 # MW of in kg/kmol
    sig = np.array([[sp.sig for sp in speclist]]*N) # sig of in Angstroms
    e = np.array([[sp.e for sp in speclist]]*N) # e is e/kB of in K
    Y = np.array([[sp.Y for sp in speclist]]*N) # Mass fractions
    X = (MWmix*Y)/MW # Mole fractions

    # Correctino factor = 1 if Ya and Yb are both nonzero, = 0 if either Ya or Yb are zero
    Correction = (Y*Y.T)/((Y*Y.T)+1e-12)

    MWAB = 2/((1/MW)+(1/MW.T)) # Molecular weights AB
    sAB = (sig+sig.T)/2 # Collision diameters AB
    Tst = T/np.sqrt(e*e.T) # Dimensionless temperature
    A,B,C,D,E,F,G,H = (1.06036,0.15610,0.19300, # collison integral paramaters
                      0.47635,1.03587,1.52996,
                      1.76474,3.89411)
    OmD = (A/(Tst**B))+(C/np.exp(D*Tst))+(E/np.exp(F*Tst))+(G/np.exp(H*Tst)) #Collision integral
    biD = Correction*0.0266*(T**1.5)/(P*np.sqrt(MWAB)*(sAB**2)*OmD) # Binary coefficients
    np.set_printoptions(suppress=True)

    D = ((1-Y)/np.sum((1-Kdelt)*X/(biD.T+1e-12),axis=1))[0]

    return D

def Vdiff(spec,speclist,state):
    # Mass diffusion velocity with mixture-averaged diffusion
    if state.multi == 1:
        Coeff = -state.D[spec.idx]/(spec.Y*state.MWmix+1e-12)
        Term1 = state.MWmix * spec.dY
        Term2 = spec.Y * state.dMWmix
        v = Coeff*(Term1 + Term2)
    else:
        i = spec.idx
        js = [J for J in range(len(speclist)) if J != i]
        Coeff = state.P*spec.MW/(state.MWmix*state.R*state.rho*state.T*spec.Y+1e-12)
        Sumterm = sum([state.D[i,j]*(state.MWmix*speclist[j].dY + speclist[j].Y*state.dMWmix) for j in js])
        v = Coeff*Sumterm
    return v

def mdots(T, phi, speclist):
    #Uses two-step quasi-global reaction from Simplified Reaction Mechanisms
    #for the Oxidation of Hydrocarbon Fuels in Flames by C. Westbrook et. al.

    F,o2,n2,co2,h2o,co = [sp.con for sp in speclist] # concetrations in mol/cm^3
    conlist = [F,o2,n2,co2,h2o,co]
    # Ea is normalized to the gas constant so has units of K
    A,Ea,a,b,m,n = {'CH4': [2.8e9,24358.,-0.3,1.3,1.,4.],
                'C3H8':[1.0e12,15098.,0.1,1.65,3.,8.]
               }[speclist[0].name] # speclist[0].name is the key of the fuel

    k1 = A*np.exp(-Ea/T)*(F**a)*(o2**b)
    with np.errstate(invalid='ignore'):
        k2f = np.nan_to_num((10**14.6)*np.exp(-20129/T)*co*(h2o**.5)*(o2**0.25))
    k2r = (5e8)*np.exp(-20129/T)*co2

    # Stoichimetric coefficients for the reactions
    # Rows are the reactions (i) and columns are the species (j)
    # j = [0,1,2,3,4,5] --> [Fuel,O2,N2,CO2,H2O,CO]
    # When indexing it is v[reaction,species]
    vjip = np.array([[1.,.5*n+.25*m, 0., 0., 0., 0.],
                     [0.,   0.5, 0., 0., 0., 1.]])

    vjipp = np.array([[0., 0., 0., 0., .5*m, n],
                      [0., 0., 0., 1., 0., 0.]])
    vji = vjipp - vjip

    q1 = k1*np.prod(np.array([conlist[j]**vjip[0,j] for j in range(6)]))
    q2 = k2f*np.prod(np.array([conlist[j]**vjip[1,j] for j in range(6)])) - k2r*np.prod(np.array([conlist[j]**vjipp[1,j] for j in range(6)]))
    # Updates the rate of mass change using production rates w_j in mol/s*m^3 for each species
    # Order : [Fuel, O2, N2, CO2, H2O, CO]
    for j in range(6):
        speclist[j].mdot = 1000*speclist[j].MW*(q1*vji[0,j]+q2*vji[1,j])

def dY(i,x,state,speclist):
    spec = speclist[i]
    N = [J for J in range(len(speclist)) if J != i]

    Coeff = 1/(state.mpp + spec.vdiff*state.rho)

    if x > 0:
        dvdx = (spec.vdiff - spec.vold)/x
    else:
        dvdx = 0

    Term1 = - spec.vdiff*spec.Y*state.d_rho
    Term2 = - state.rho*spec.Y*dvdx
    return Coeff*(spec.mdot - Term1 - Term2)

def Update(speclist,state,x,**ipts):

    Y_key = ['YF','YO2','YN2','YCO2','YH2O','YCO']
    state.T = ipts['T']
    T=state.T
    dT = ipts['dT']
    Y = [ipts[key] for key in Y_key]

    recip_MWmix = 0
    dMWmix = 0
    for i in range(len(speclist)):
        speclist[i].Y = Y[i]
        speclist[i].dY = dY(i,x,state,speclist)
        recip_MWmix += speclist[i].Y/speclist[i].MW
        dMWmix += speclist[i].dY/speclist[i].MW

    state.MWmix = 1/recip_MWmix
    state.dMWmix = -dMWmix/(state.MWmix**2)
    state.rho = state.P*state.MWmix/(state.R*T)
    state.d_rho = state.P*(T*state.dMWmix - state.MWmix*dT)/(state.R*T**2)
    state.Dold = state.D
    state.D = {0:Dij,1:Dim}[state.multi](T,state.P,state.MWmix,speclist)

    for sp in speclist:
        sp.Cp = Cp_T(sp,T) # Update heat capacities
        sp.h = h_T(sp,T) # Update enthalpies
        sp.con = state.rho*sp.Y/(sp.MW*1000) # Update concentrations
        sp.k = k_T(sp,speclist,state)#state.rho*state.D[sp.idx]*sp.Cp #Update conductivity
        sp.vold = sp.vdiff
        sp.vdiff = (sp.Y/(sp.Y+1e-12))*Vdiff(sp,speclist,state)#state.D,state.MWmix,state.dMWmix,state.rho)

    state.totCp = sum([sp.Cp*sp.Y for sp in speclist]) #Update total specific heat
    state.kold = state.totk
    state.totk = totk_T(speclist,state.MWmix) # Update total conductivity
    mdots(T,state.phi,speclist) # Update mass change rates

    state.update = True

def rkstep(ode,ic,state,x=0.001,dx=0.001):
    var = list(ode.keys())
    y0 = list(ic.values())#np.array(zip(*sysf.values())[0])
    eqs = list(ode.values())#list(zip(*sysf.values())[1])

    y1 = dict(zip(var,y0))
    k1 = np.array([f(x,**y1) for f in eqs])
    y2 = dict(zip(var,.5*dx*k1+y0))
    k2 = np.array([f(x+.5*dx,**y2) for f in eqs])

    y3 = dict(zip(var,.5*dx*k2+y0))
    k3 = np.array([f(x+.5*dx,**y3) for f in eqs])

    y4 = dict(zip(var,dx*k3+y0))
    k4 = np.array([f(x+dx,**y4) for f in eqs])

    dy = ((1/6.)*dx*(k1+2*k2+2*k3+k4))
    yout = dict(zip(var,y0+dy))
    return yout

def Dij(T,P,MWmix,speclist):

    # This calculates the multicomponent diffusion coefficients in m^2/s
    # T is Temperature in K
    # P is Pressure in Pa
    # speclist is a list of chemical species class instances from comb_sim.py
    N = len(speclist) # Number of species
    Kdelt = np.identity(N) # Identity matrix used as a Kronecker Delta
    MW = np.array([[sp.MW for sp in speclist]]*N)*1000 # MW of in kg/kmol
    sig = np.array([[sp.sig for sp in speclist]]*N) # sig of in Angstroms
    e = np.array([[sp.e for sp in speclist]]*N) # e is e/kB of in K
    Y = np.array([[sp.Y for sp in speclist]]*N) # Mass fractions
    X = (MWmix*Y)/MW + 1e-12# Mole fractions the 1e-12 is added to avoid singular matrices

    # Correctino factor = 1 if Ya and Yb are both nonzero, = 0 if either Ya or Yb are zero
    Correction = (Y*Y.T)/((Y*Y.T)+1e-12)

    MWAB = 2/((1/MW)+(1/MW.T)) # Molecular weights AB
    sAB = (sig+sig.T)/2 # Collision diameters AB
    Tst = T/np.sqrt(e*e.T) # Dimensionless temperature
    A,B,C,D,E,F,G,H = (1.06036,0.15610,0.19300, # collison integral paramaters
                      0.47635,1.03587,1.52996,
                      1.76474,3.89411)
    OmD = (A/(Tst**B))+(C/np.exp(D*Tst))+(E/np.exp(F*Tst))+(G/np.exp(H*Tst)) #Collision integral
    biD = 0.0266*(T**1.5)/(P*np.sqrt(MWAB)*(sAB**2)*OmD) # Binary coefficients

    coeff = X/(biD*MW.T)[:,None]
    term_1 = (MW*X)[:,:,None]*(1-Kdelt)[:,None]
    term_2 = (MW*X).T[:,:,None]*(Kdelt[:,:,None]-Kdelt)

    L = np.sum(coeff*(term_1-term_2),axis=2) # Compute the L matrix
    F = np.linalg.inv(L) # Take the matrix inverse of L to find F
    Fii=np.array([F.diagonal()]*N)#.T

    D = Correction*X.T*(MWmix/MW)*(F-Fii)

    return D
