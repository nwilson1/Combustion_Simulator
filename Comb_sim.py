## Noah Wilson
## 12/10/2020
## 4/8/2021 -- Updated a few lines to work with Python 3.
## Main script for combustion simulation.

import numpy as np
import matplotlib.pyplot as plt
import Comb_aux as aux

## UI ##
Fuelinpt = input('\nChoose Fuel (CH4:0 C3H8:1) Enter for default (CH4) : ') or 0
Oxinpt = input('\nChoose Oxidizer (Air:0 O2:1) Enter for default (Air) : ') or 0
Phiinpt = input('\nInput equivalence ratio between 0 and 1.  Enter for default (1) : ') or 1
init_Pressure = input('\nInput pressure in Pa. Enter for default (101 kPa) : ') or 101325.
init_Temp = input('\nInput initial temperature in  K. Enter for default (298 K) : ') or 298
Multi = input('\nMulticompoent (0) or mixture averaged (1) diffusion? Enter for default (Averaged) : ') or 1
mflux = input('\nInput inital guess mass flux rate. Enter for default (0.001 kg/(s m^2)) : ') or 0.001
print('\nWarning -- Becomes unstable after multiple steps, only one step is recommended.')
Iterations = input('Enter an integer for number of ODE steps. Enter for default (1) : ') or 1
print('')

## INITIALIZING THE SYSTEM ##

class System:

    _registry = [] # Updates a list of all species available
    List = None # List of species used in current simulation
    def __init__(self,name=None,idx=None,MW=0,sig=None,e=None,type=None,a_stoich=None):

        #SPECIES CONSTANTS#
        self.name = name # Species name as string
        self.idx = idx # species index for the i and j values
        if name != None:
            self._registry.append(self)
            self.As_HighT,self.As_LowT,self.hf0 = aux.get_as(name)
        self.MW = MW*0.001 # Moleucular weight in kg/mol
        self.sig = sig # hard-sphere collision diameter in Angstroms
        self.e = e # Lennard-Jones energy (e/kB) in K
        self.type = type # Species type (fuel, reactant, product)
        self.a_stoich = a_stoich # a-Stoichiometric for a given fuel

        #SPECIES VARIABLES#
        self.mdot = 0 # Rate of mass change m_i''' in kg/s*m^3
        self.vdiff = 0
        self.vold = 0
        self.Y = 0 # Mass fractions
        self.dY = 0 # Mass fraction gradient dY/dx
        self.con = 1 # Species concentrations in mol/m^3
        self.Cp = 1 # Species heat capacity in J/mol*K
        self.k = 1 # Species thermal conductivity in J/m*K
        self.h = 1 # Equal to hf0(298) + sensible enthalpy change in J/mol

        #STATE VARIABLES#
        self.R = 8.31446 # Gas constant in J/mol*K
        self.P = 101325.  # Pressure in Pa, constant due to conserv. of momentum
        self.phi = 1 # Equivalence ratio
        self.T = 298. #Temperature in K
        self.dT = 0
        self.rho = 1 # System density in kg/m^3
        self.d_rho = 0
        self.mpp = .001 # Net mass flow m'' in kg/s*m^2
        self.mtot = 0 # Total mass of the system in kg
        self.totCp = 1 # Total system heat capacity in J/mol*K
        self.totk = 0 #Total system thermal conductivity in
        self.kold = 0 #Prevous value of k to estimate dk/dx
        self.D = 0#np.zeros((6,6)) # Multicomponent diffusion coefficeints in m^2/s
        self.MWmix = 1  # Molecular weight of the mixture in kg/mol
        self.dMWmix = 0
        self.counter = 0
        self.update = False
        self.multi = 1

state = System()
state.T = init_Temp
state.P = init_Pressure
state.phi = Phiinpt
state.mpp = mflux
state.multi = Multi

CH4  = System(name='CH4',idx=0,MW=16.043,sig=3.758,e=148.6,type='Fuel',a_stoich=2.)
C3H8 = System(name='C3H8',idx=0,MW=44.096,sig=5.118,e=237.1,type='Fuel',a_stoich=5.)
Fuels = [sp for sp in System._registry if sp.type=='Fuel']
Fuel = Fuels[Fuelinpt] # Chooses from list of fuels using user input

O2   = System(name='O2',idx=1,MW=31.998,sig=3.467,e=106.7,type='Reac')
N2   = System(name='N2',idx=2,MW=28.013,sig=3.798,e=71.40,type='Reac')

a = state.phi*Fuel.a_stoich # Moles of oxidizer based on fuel and phi
n_N2 = a*[3.76 , 0][Oxinpt] # Amount of N2, 3.76a if Ox is air, 0 if Ox is O2
state.m_tot = Fuel.MW + a*O2.MW+n_N2*N2.MW
Fuel.Y = Fuel.MW/state.m_tot  # Initial Y of fuel
O2.Y = a*O2.MW/state.m_tot    # Initial Y of O2
N2.Y = n_N2*N2.MW/state.m_tot  # Initial Y of N2

CO2  = System(name='CO2',idx=3,MW=44.010,sig=3.941,e=195.2,type='Prod')
H2O  = System(name='H2O',idx=4,MW=18.015,sig=2.641,e=809.1,type='Prod')
CO   = System(name='CO',idx=5,MW=28.010,sig=3.690,e=91.7,type='Prod')

System.List = [Fuel,O2,N2,CO2,H2O,CO] # Updates list of species in simulation
Diff_choice = {0:aux.Dij,1:aux.Dim}
state.D = Diff_choice[state.multi](state.T,state.P,state.MWmix,System.List)

## DIFFERENTIAL EQUATIONS ##

# Species conservation
class Y_i:

    Y_key = ['YF','YO2','YN2','YCO2','YH2O','YCO']

    def __init__(self,i):
        self.I = i

    def dYdx(self,x,m,**ipts):

        if state.update == False:
            aux.Update(System.List,state,x,**ipts)
        T = ipts['T']
        dT = ipts['dT']
        spec = System.List[self.I]
        N = [J for J in range(len(System.List)) if J != self.I]

        Coeff = 1/(state.mpp + spec.vdiff*state.rho)

        if x > 0:
            dvdx = (spec.vdiff - spec.vold)/x
        else:
            dvdx = 0

        Term1 = - spec.vdiff*spec.Y*state.d_rho
        Term2 = - state.rho*spec.Y*dvdx
        dY = Coeff*(spec.mdot - Term1 - Term2)
        return dY

# Mass conservation
def dm(x,**ipts):
    if state.update == False:
        aux.Update(System.List,state,x,**ipts)
    return 0
# Energy conservation
def ddT(x,m,**ipts):
    if state.update == False:
        aux.Update(System.List,state,x,**ipts)
    T = ipts['T']
    dT = ipts['dT']
    subTerm1 = m*state.totCp
    subTerm2 = -(state.totk-state.kold)/x
    subTerm3 = state.rho*sum([sp.Y*sp.vdiff*sp.Cp for sp in System.List])

    Term1 = dT*(subTerm1 + subTerm2 + subTerm3)
    Term2 = sum([sp.h*sp.mdot for sp in System.List])

    #print 'ddT'
    return (Term1 + Term2)/state.totk

def dT(x,**ipts):
    if state.update == False:
        aux.Update(System.List,state,x,**ipts)
    #print 'dT'
    return ipts['dT']

## SOLVING THE ODES ##

ICs = {
          'm':state.mpp,
          'dT':state.dT,
          'T':state.T
          }
ODEs = {
          'm':dm ,
          'dT':ddT ,
          'T':dT
          }

dYfunc = [0]*len(System.List)
for I in range(len(System.List)):

    dYfunc[I] = Y_i(I).dYdx
    ICs[Y_i.Y_key[I]] = System.List[I].Y
    ODEs[Y_i.Y_key[I]] = dYfunc[I]

for i in range(Iterations):
    step = aux.rkstep(ODEs,ICs,state)
    state.update = False
    for key in step.keys():
        ICs[key] = step[key]

print( '\n\n ---- SYSTEM PARAMETERS FOLLWING FINAL ITERATION ---- \n\n')
print( '\nTemperature in K : ', ICs['T'])
print( '\nTotal mass flux rate in kg/(s m^2) : ' ,ICs['m'])
print( '\nMass fracitions (F is fuel): \n')
for key in sorted(ICs.keys()):
    if key in Y_i.Y_key:
        print( '     ',key,' : ',round(ICs[key],4))
print( '\nSpecific heat capacities in J/(mol K):\n\n')
for sp in System.List:
    print( '     ',sp.name,' : ',round(sp.Cp,2))
print( '\nThermal conductivities in J/(m K):\nNote: This is only calculated if the species is present in the system.\n\n')
for sp in System.List:
    print( '     ',sp.name,' : ',round(sp.k,2))
print( '\Mass production/use rates in kg/(s m^3) : \n')
for sp in System.List:
    print( '     ',sp.name,' : ',sp.mdot)
print( '\nDiffusion coefficients in m^2/s :\nNote: The diffusion coefficient is set to zero for any species not present in the system.\nSpecies order is [Fuel O2 N2 CO2 H2O CO]\n\n',state.D,'\n')
