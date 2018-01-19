#!/usr/bin/env python

# Importing everything python needs in order to be smart.
import numpy as np                                      # Numpy tells python how to work with numbers.
import matplotlib.pyplot as plt                         # Matplotlib tells python how to make pretty plots.
import ionospheres_programs as ip                       # Importing my custom library.
import iono                                             # Import the Fortran Subroutines
import shelve                                           # For saving calculated variables/arrays to an external file so that you don't have to keep running your code over again becaue it takes so darn long! 

#import time

#start = time.time()

#***********************************************************************
# Declaring some variables.                                            *
#***********************************************************************

zsize = 1000.0
ztop = zsize - 1.0
z = np.linspace(100, 800, zsize)*1E3                                    # Altitude ranges from 100,000 [m] - 800,000 [m] in increments of 1000.
dz = np.zeros(zsize) ; dz[1:] = z[1:] - z[0:-1] ; dz[0] = dz[1]         # Defining a delta z in [m].
theta = np.linspace(0, np.pi, zsize)                                    # Theta ranges from 0 [rad] - pi [rad] in incriments of zsize. 
T = 800.0*np.tanh(theta*3.0)+200.0                                      # Temperature ranges from 200-1000 [K] in increments of zsize.
Te = 1.0*T
Ti = 1.0*T
dT = np.zeros(zsize) ; dT[1:] = T[1:] - T[0:-1] ; dT[0] = dT[1]         # Temperature increments.
R = 6378.1*1E3                                                          # Radius of Earth in [m].
g0 = 9.81                                                               # Gravity at Earths surface in [m/s^2].
g = g0*(R/(R+z))*(R/(R+z))                                              # Gravity as it varies with altitude in [m/s^2].
k = 1.3806503E-23                                                       # Boltzmann constant in [((m^2)*kg)/((s^2)*K)].
A = 6.022E23                                                            # Avagadro's Number.
planck = 6.626068E-34                                                   # Planck's Constant in [((m^2)*kg)/s].
c = 3.0E8                                                               # Speed of Light in [m/s].
P = 150.0                                                               # F10.7 proxy in [SFU].
eps = 0.25                                                              # Neutral heating efficiency.
eps_e = 0.04                                                            # Electron heating efficiency.
eps_i = 0.0034                                                            # Ion heating efficiency.                                         
dayseconds = 86400.0                                                    # number of seconds in a day.
days = 2.0                                                              # Number of days you want to calculate over.
dt = 60.0                                                               # Definin time incriment as dt [s].
tsize = dayseconds*days/dt                                              # Number of second in a day divided by dt intervals gives us the number of "angular steps" used to calculate chi over the course of a day.
tloop = int(tsize)                                                      # Defining an integer range value for the time dependent loop.
time = np.arange(tsize)*dt                                              # Time
#X = -(np.pi/2.0)                                                        # Setting an initial solar zenith angle value to zero (ie. 12:00PM ~ noon).
X = 0.0
sza = np.cos(X)                                                         # Define the cosine of the solar zenith angle for use in calculating the optical depth.
me = 9.10938188E-31                                                     # Mass of electron in [kg]


# Creating a list of all of the neutral species that live in our
# virtual thermosphere so that we can later assign properties 
# to them.
neutrals = ['O','N','O2','N2','NO']

# Creating a dictionary of masses for each species. 
# Divide by (1000.0*A) to convert from [amu] to [kg].
m = {'O':(16.0/(1000.0*A)),
     'N':(14.0/(1000.0*A)),
     'O2':(32.0/(1000.0*A)),
     'N2':(28.0/(1000.0*A)),
     'NO':(30.0/(1000.0*A))}

m_O = 16.0/(1000.0*A)
m_N = 14.0/(1000.0*A)
m_O2 = 32.0/(1000.0*A)
m_N2 = 28.0/(1000.0*A)
m_NO = 30.0/(1000.0*A)

# Creating a dictionary of initial number density values at the base 
# of the thermospehre. Densities of Each Species are given in units 
# of [m^-3] at 100km.
n_base = {'O':4.9E17,'N':6.4E11,'O2':2.4E18,'N2':1.0E19,'NO':1.0E14}

# Calculating the average mass density at the base of the thermosphere
# in unis of [kg/m^3].
rho_tot = (n_base['O']*m['O'] + 
           n_base['N']*m['N'] + 
           n_base['O2']*m['O2'] + 
           n_base['N2']*m['N2'] + 
           n_base['NO']*m['NO'])

# Calling a definition inside the custom library which fetches the
# solar flux data from Aaron's poorly formatted table.
data = ip.read_solar_flux()

# Using the data fetched to define some more variables
lamb = ((data['lower'] + data['upper'])/2.0)*1E-10      # Use the Average Solar Spectrum angstroms to meters.
e = (planck*c)/lamb                                     # Photon energy.
data['F74113'] *= 1E9*1E4                               # Convert 1E9 photons/cm^2*s to photons/m^2*s.
I_inf = data['F74113']*(1.0 + data['Abi']*(P-80.0))     # Solar irradiance at infinity.

H_ave = np.zeros(29.0)

# Creating an empty dictionary to store our number density arrays.
n = {}        
# Creating an empty array for each species, s.
for s in neutrals:
    n[s] = np.zeros(1000.0)

# Creating an empty dictionary to store our column density arrays.
col_n = {}
# Creating an empty array for each species, s.
for s in neutrals:
    col_n[s] = np.zeros(1000.0)

# Creating an empty dictionary to store our optical depth arrays.
tau = {}
for s in neutrals:
    tau[s] = np.zeros((1000.0,37.0)) 

sumtau = np.zeros((1000.0,37.0))

I = np.zeros((1000.0,37.0))

# Creating an empty dictionary to fill using numpy and calculated Q values.
Q = {}
for s in neutrals:
    Q[s] = np.zeros(1000.0)

# Creating an empty dictionary to fill using numpy and summing up the Q values.
sumQ = np.zeros(1000.0)

#***********************************************************************
# Doing the initial calculations using the initial temperature guess   *
# at solar zenith angle zero.                                          *
#***********************************************************************

# Calculating the number density one last time with our new temperature 
# outside the loop so that we can plot everything.
H = {'O':((k*T)/(m['O']*g)),
     'N':((k*T)/(m['N']*g)),
     'O2':((k*T)/(m['O2']*g)),
     'N2':((k*T)/(m['N2']*g)),
     'NO':((k*T)/(m['NO']*g))}    

# Fetching the average scale heights from the custom libaray.
H_ave = ip.calculate_average_scale_height(H_ave, n_base, m, H, rho_tot)

# Fetching the number densities from the custom libaray.
# n = ip.calculate_number_densities(T, z, n_base, m, k, g)
#n = ip.calculate_number_densities(n, T, z, H, n_base)
n['O'] = iono.nout(n_base['O'], T, z, H['O'])
n['N'] = iono.nout(n_base['N'], T, z, H['N'])
n['O2'] = iono.nout(n_base['O2'], T, z, H['O2'])
n['N2'] = iono.nout(n_base['N2'], T, z, H['N2'])
n['NO'] = iono.nout(n_base['NO'], T, z, H['NO'])

# Creating a dictionary to store the mass densities of each species as
# they vary with altitude in [kg/m^3].
rho = {'O':(n['O']*m['O']),
       'N':(n['N']*m['N']),
       'O2':(n['O2']*m['O2']),
       'N2':(n['N2']*m['N2']),
       'NO':(n['NO']*m['NO'])}

# Calculating and array containing the total number density as it
# changes with altitude in [m^-3]. 
N_tot = n['O'] + n['N'] + n['O2'] + n['N2'] + n['NO']

# Calculating the total weighted mass as it changes with altitude
# in units of [kg].
M = (m['O']*(n['O']/N_tot) + 
     m['N']*(n['N']/N_tot) + 
     m['O2']*(n['O2']/N_tot) + 
     m['N2']*(n['N2']/N_tot) + 
     m['NO']*(n['NO']/N_tot))

# Total weighted molecular mass as it changes with altitude in [amu].
M_molec = M*(1000.0*A)

# Total mass density
RHO = N_tot*M    

# Heat capacity
Cp = 2000.0

# Fetching the column densities from the custom library.
#col_n = ip.calculate_column_densities(col_n, n, H, z)

col_n['O'] = iono.col_nout(n['O'], H['O'], z)
col_n['N'] = iono.col_nout(n['N'], H['N'], z)
col_n['O2'] = iono.col_nout(n['O2'], H['O2'], z)
col_n['N2'] = iono.col_nout(n['N2'], H['N2'], z)
col_n['NO'] = iono.col_nout(n['NO'], H['NO'], z)


# Defining a column mass density just for fun becuase I want to see 
# what the plot looks like and I feel like it is an important quantity 
# in atmospheric science.
col_rho = {'O':(col_n['O']*m['O']),
           'N':(col_n['N']*m['N']),
           'O2':(col_n['O2']*m['O2']),
           'N2':(col_n['N2']*m['N2']),
           'NO':(col_n['NO']*m['NO'])}

# Fetching the optical depths from the custom library.
tau = ip.calculate_optical_depth(tau, data, col_n, zsize)

# Fetching the combined optical depths from the custom library.    
sumtau = ip.combine_optical_depths(sumtau, tau, sza)

# Fetching the intensity from the custom library.
I = ip.calculate_intensity(I, I_inf, sumtau)

# Fetching the energy dissipation from the custom library.
Q = ip.calculate_energy_dissipation(Q, eps, n, data, I, e)

# Fetching the optical depths from the custom library.
Q_nr = np.zeros(zsize)
Q_r = np.zeros(zsize)
Q_en = np.zeros(zsize)
sumQ = ip.combine_energy_dissipation(sumQ, Q, Q_nr, Q_r, Q_en)

# Defineing the thermal conduction in [W/(m*K)].
lam_n = 3.6E-4*np.power(T, .75)

# Defining an array to fill with the maximum temperature value for 
# each zenith angle iteration and setting the zeroth value. This will 
# later be used to make a cool plot of maxumum temperature vs time.
T_max = np.zeros(tloop) 
T_max[0] = T[999]

ni_noon = {}

#***********************************************************************
# Building the ionosphere                                              *
#***********************************************************************

#Creating a list of all of the ion species tht live in our virtual
#thermosphere so that we can later assign properties to them.
ions = ['O+', 'N+', 'O2+', 'N2+', 'NO+', 'O+_new']

ni = {}
for s in ions:
    ni[s] = np.zeros(1000)

# Defining an ion production rate.
Ps = {}
for s in ions:
    Ps[s] = np.zeros(1000.0)

# The rest of the ionosphere is in the loop.
ni_max = {}
for i in ions:
    ni_max[i] = np.zeros(tloop)
    ni_max[i][0] = 0

# Creating an empty dictionary to store our number density arrays.
ni_dawn = {}        
# Creating an empty array for each species, s.
for i in ions:
    ni_dawn[i] = np.zeros(1000.0)

# Creating an empty dictionary to store our number density arrays.
ni_noon = {}        
# Creating an empty array for each species, s.
for i in ions:
    ni_noon[i] = np.zeros(1000.0)

# Creating an empty dictionary to store our number density arrays.
ni_dusk = {}        
# Creating an empty array for each species, s.
for i in ions:
    ni_dusk[i] = np.zeros(1000.0)

# Creating an empty dictionary to store our number density arrays.
ni_midnight = {}        
# Creating an empty array for each species, s.
for i in ions:
    ni_midnight[i] = np.zeros(1000.0)

# Defiing lists/variables for vertical advection specifically as a loss of O+
GOm = {}
GOp = {}
GOm['O+'] = np.zeros(1000)
GOp['O+'] = np.zeros(1000)
v_Op = -2.0
vp = max([v_Op,0])
vm = min([v_Op,0])
dr1 = np.zeros(1000)

# Defining an array to fill with the maximum temperature value for 
# each zenith angle iteration and setting the zeroth value. This will 
# later be used to make a cool plot of maxumum temperature vs time.
Te_max = np.zeros(tloop) 
Te_max[0] = Te[999]

# Defining an array to fill with the maximum temperature value for 
# each zenith angle iteration and setting the zeroth value. This will 
# later be used to make a cool plot of maxumum temperature vs time.
Ti_max = np.zeros(tloop) 
Ti_max[0] = Ti[999]

# Defining and array to fill with maximum ion density values 
# and calling it the F2 peak.
F2_peak = np.zeros(tloop)

#***********************************************************************
# Definining Loss Variables                                            *
#***********************************************************************

Lr_N2 = np.zeros(zsize)
Lr_O2 = np.zeros(zsize)
QTe = np.zeros(zsize)
Lv_O2 = np.zeros(zsize)
Le_O = np.zeros(zsize)
Lf_O = np.zeros(zsize)
d = np.zeros(zsize)
v_OpO = np.zeros(zsize)
v_N2pN2 = np.zeros(zsize)
v_O2pO2 = np.zeros(zsize)
v_ei_Op = np.zeros(zsize)
v_ei_Np = np.zeros(zsize)
v_ei_O2p = np.zeros(zsize)
v_ei_N2p = np.zeros(zsize)

Q_ei = np.zeros(zsize)
Q_en = np.zeros(zsize)
Q_nr = np.zeros(zsize)
Q_r = np.zeros(zsize)


# N2 vibration coeffients from Tables 9.3 and 9.5 on page 279 of Schunk and Nagy       
T93 = np.zeros(10, dtype=[('A', float), ('B', float), ('C', float), ('D',float), ('F', float)])
T93[0] = (2.025, 8.752e-4, 2.954e-7, -9.562e-11, 7.252e-15)
T93[1] = (-7.066, 1.001e-2, -3.066e-6, 4.436e-10, -2.449e-14)
T93[2] = (-8.211, 1.092e-2, -3.369e-6, 4.891e-10, -2.706e-14)
T93[3] = (-9.713, 1.204e-2, -3.732e-6, 5.431e-10, -3.008e-14)
T93[4] = (-10.353, 1.243e-2, -3.850e-6, 5.600e-10, -3.100e-14)
T93[5] = (-10.819, 1.244e-2, -3.771e-6, 5.385e-10, -2.936e-14)
T93[6] = (-10.183, 1.185e-2, -3.570e-6, 5.086e-10, -2.769e-14)
T93[7] = (-12.698, 1.309e-2, -3.952e-6, 5.636e-10, -3.071e-14)
T93[8] = (-14.710, 1.409e-2, -4.249e-6, 6.058e-10, -3.300e-14)
T93[9] = (-17.538, 1.600e-2, -4.916e-6, 7.128e-10, -3.941e-14)

T95 = np.zeros(8, dtype=[('A', float), ('B', float), ('C', float), ('D',float), ('F', float)])
T93[0] = (-3.413, 7.326e-3, -2.200e-6, 3.128e-10, -1.702e-14)
T93[1] = (-4.160, 7.803e-3, -2.352e-6, 3.352e-10, -1.828e-14)
T93[2] = (-5.193, 8.360e-3, -2.526e-6, 3.606e-10, -1.968e-10)
T93[3] = (-5.939, 8.807e-3, -2.669e-6, 3.806e-10, -2.073e-14)
T93[4] = (-8.261, 1.010e-2, -3.039e-6, 4.318e-10, -2.347e-14)
T93[5] = (-8.185, 1.010e-2, -3.039e-6, 4.317e-10, -2.347e-14)
T93[6] = (-10.823, 1.199e-2, -3.620e-6, 5.159e-10, -2.810e-14)
T93[7] = (-11.273, 1.283e-2, -3.879e-6, 5.534e-10, -3.016e-14)

first_sum = np.zeros(zsize, dtype=float)
second_sum = np.zeros(zsize, dtype=float)
Le_N2 = np.zeros(zsize, dtype=float)
E1 = 3353

#***********************************************************************
# Beginning the time dependent loop in which solar zenith angle        *
# advances for each time step.                                         *
#***********************************************************************

for loop in range(tloop):
    # Calling the new temperature profil from the custom library and 
    # recalculating everything again.
    T = ip.T_new(T, sumQ, dz, lam_n, dt, RHO, Cp)
    
    # Filling in the maximum temperature value for each iteration.
    T_max[loop] = T[999]
    
    # Creating a dictionary of scale height's for each species in [m].
    H = {'O':((k*T)/(m['O']*g)),
         'N':((k*T)/(m['N']*g)),
         'O2':((k*T)/(m['O2']*g)),
         'N2':((k*T)/(m['N2']*g)),
         'NO':((k*T)/(m['NO']*g))}    
    
    # Fetching the average scale heights from the custom libaray.
    H_ave = ip.calculate_average_scale_height(H_ave, n_base, m, H, rho_tot)
    
    # Fetching the number densities from the custom libaray.
    # n = ip.calculate_number_densities(T, z, n_base, m, k, g)    
    #n = ip.calculate_number_densities(n, T, z, H, n_base)
    n['O'] = iono.nout(n_base['O'], T, z, H['O'])
    n['N'] = iono.nout(n_base['N'], T, z, H['N'])
    n['O2'] = iono.nout(n_base['O2'], T, z, H['O2'])
    n['N2'] = iono.nout(n_base['N2'], T, z, H['N2'])
    n['NO'] = iono.nout(n_base['NO'], T, z, H['NO'])    
    
    # Creating a dictionary to store the mass densities of each species
    # as they vary with altitude in [kg/m^3].
    rho = {'O':(n['O']*m['O']),
           'N':(n['N']*m['N']),
           'O2':(n['O2']*m['O2']),
           'N2':(n['N2']*m['N2']),
           'NO':(n['NO']*m['NO'])}
    
    # Calculating and array containing the total number density as it
    # changes with altitude in [m^-3].
    N_tot = n['O'] + n['N'] + n['O2'] + n['N2'] + n['NO']
    
    # Calculating the total weighted mass as it changes with altitude 
    # in units of [kg].
    M = (m['O']*(n['O']/N_tot) + 
         m['N']*(n['N']/N_tot) + 
         m['O2']*(n['O2']/N_tot) + 
         m['N2']*(n['N2']/N_tot) + 
         m['NO']*(n['NO']/N_tot))
    
    # Total mass density
    RHO = N_tot*M
    
    # Heat capacity
    Cp = 2000.0
    
    # Total weighted molecular mass as it changes with altitude 
    # in units of [amu].
    M_molec = M*(1000.0*A)
    
    # Fetching the column densities from the custom library.
    #col_n = ip.calculate_column_densities(col_n, n, H, z)
    
    col_n['O'] = iono.col_nout(n['O'], H['O'], z)
    col_n['N'] = iono.col_nout(n['N'], H['N'], z)
    col_n['O2'] = iono.col_nout(n['O2'], H['O2'], z)
    col_n['N2'] = iono.col_nout(n['N2'], H['N2'], z)
    col_n['NO'] = iono.col_nout(n['NO'], H['NO'], z)
    
    
    # Defining a column mass density just for fun becuase I want to  
    # see what the plot looks like and I feel like it is an important 
    # quantity in atmospheric science.
    col_rho = {'O':(col_n['O']*m['O']),
               'N':(col_n['N']*m['N']), 
               'O2':(col_n['O2']*m['O2']),
               'N2':(col_n['N2']*m['N2']), 
               'NO':(col_n['NO']*m['NO'])}
    
    # Setting the new value of solar zenith angle over each iteration.
    X += (dt/dayseconds)*2*np.pi
    
    # Define the new value of the cosine of the solar zenith angle 
    # for each iteration.
    sza = np.cos(X)
    
    # Setting up a minimum value for the solar zenith angle which is 
    # equivalent to recognizing that there is still ionization on the 
    # night side due to a flux of galactic cosmic rays etc...  
    if sza <= 0.05:
        sza = 0.05
        
    # Fetching the optical depths from the custom library.
    tau = ip.calculate_optical_depth(tau, data, col_n, zsize)
    
    # Fetching the combined optical depths from the custom library.    
    sumtau = ip.combine_optical_depths(sumtau, tau, sza)
    
    # Fetching the intensity from the custom library.
    I = ip.calculate_intensity(I, I_inf, sumtau)
    
    # Fetching the energy dissipation from the custom library.
    Q = ip.calculate_energy_dissipation(Q, eps, n, data, I, e)
    
    # Fetching the optical depths from the custom library.
    sumQ = ip.combine_energy_dissipation(sumQ, Q, Q_nr, Q_r, Q_en)
    
    # Defineing the thermal conduction in [W/(m*K)].
    lam_n = 3.6E-4*np.power(T, .75)
    
    #*******************************************************************
    # Building the rest of the ionosphere                             *
    #*******************************************************************
    
    #N+ Ionization Rates (souces) ; Number density of neutral reactant (losses)
    R37 = 3.07E-10*1E-6 ; M37 = n['O2']     # N+ + O2 -> O2+ + N
    R38 = 2.32E-10*1E-6 ; M38 = n['O2']     # N+ + O2 -> NO+ + O 
    R39 = 4.6E-11*1E-6 ; M39 = n['O2']      # N+ + O2 -> O+ + NO
    
    #N2+ Ionization Rates (souces) ; Number density of neutral reactant (losses)
    R43 = 4.1E-10*1E-6 ; M43 = n['NO']      # N2+ + NO -> NO+ + N2
    R44 = 1.3E-10*1E-6 ; M44 = n['O']       # N2+ + O -> NO+ + N
    R45 = 9.8E-12*1E-6 ; M45 = n['O']       # N2+ + O -> O+ + N2
    R46 = 5.0E-11*1E-6 ; M46 = n['O2']      # N2+ + O2 -> O2+ + N2
    
    #O+ Ionization Rates (souces) ; Number density of neutral reactant (losses)
    R47 = 1.2E-12*1E-6 ; M47 = n['N2']      # O+ + N2 -> NO+ + N
    R48 = 2.1E-11*1E-6 ; M48 = n['O2']      # O+ + O2 -> O2+ + O
    R49 = 8.0E-13*1E-6 ; M49 = n['NO']      # O+ + NO -> NO+ + O
    
    #O2+ Ionization Rates (souces) ; Number density of neutral reactant (losses)
    R52 = 4.6E-10*1E-6 ; M52 = n['NO']      # O2+ + NO -> NO+ + O2
    R53 = 1.5E-10*1E-6 ; M53 = n['N']       # O2+ + N -> NO+ + O
    
    # Calling the ion production rate.
    Ps = ip.ion_production_rate(Ps, n, I, data, zsize)
    
    # Ion source terms.
    source = {'O+':(Ps['O+']),
              'N+':(Ps['N+']),
              'O2+':(Ps['O2+']),
              'N2+':(Ps['N2+'])}
    
    ni_sum = ni['O+'] + ni['N+'] + ni['O2+'] + ni['N2+']    
    
    # Ion loss terms.
    RM = {'O+':(R47*M47 + R48*M48 + R49*R49 + ni_sum*3.7E-12*np.power(250.0/T, 0.7)*1E-6),
            'N+':(R37*M37 + R38*M38 + R39*M39 + ni_sum*3.6E-12*np.power(250.0/T, 0.7)*1E-6),
            'O2+':(R52*M52 + R53*M53 + ni_sum*2.4E-7*np.power(300.0/T, 0.70)*1E-6),
            'N2+':(R43*M43 + R44*M44 + R45*M45 + R46*M46 + ni_sum*2.2E-7*np.power(300.0/T, 0.39)*1E-6)}
    
    ni = ip.calculate_ion_number_densities(ni, source, RM, dt)
    
    ni_sum = ni['O+'] + ni['N+'] + ni['O2+'] + ni['N2+']
    
    #plt.plot(ni['O+'],z*1E-3)
    #plt.draw()
    
    for i in ions:
        ni_max[i][loop] = ni[i][zsize-1]
    
    F2_peak[loop] = np.max(ni_sum)
    
    #print loop
    
    #print ni_max['O+'][loop]
    
    #*******************************************************************
    # Adding vertical advection specifically as a loss of O+           *
    #*******************************************************************
    
    GOp['O+'][1:-1] = (ni['O+'][1:-1] - ni['O+'][0:-2])/dz[1:-1]
    GOp['O+'][0] = 0.0
    GOm['O+'][0:-2] = (ni['O+'][1:-1] - ni['O+'][0:-2])/dz[0:-2]
    GOm['O+'][999] = GOm['O+'][998]
    
    dr1 = -dt * (vp * GOp['O+'] + vm * GOm['O+'])
    Oplus_old = ni['O+']
    ni['O+'] = ni['O+'] + dr1
    
    # Calculating electron temperatures separately.
    #print Te
    Te = ip.Te(Te, sumQ, dz, eps, eps_e, ni_sum, N_tot, Q_ei, Q_en, Lr_N2, Lr_O2, Lv_O2, Lf_O)    
    #print Te
    # Filling in the maximum electron temperature value for each iteration.
    Te_max[loop] = Te[999]
    
    # Calculating electron temperatures separately.
    #print Ti
    Ti = ip.Ti(Ti, sumQ, dz, eps, eps_i, ni_sum, N_tot, Q_ei, Q_nr, Q_r)
    #print Ti
    # Filling in the maximum ion temperature value for each iteration.
    Ti_max[loop] = Ti[999]
    
    '''
    if loop == (days-1)*dayseconds/dt + dayseconds/dt * 0.00:
        print "dawn"
        ni_dawn = ni
    if loop == (days-1)*dayseconds/dt + dayseconds/dt * 0.25:
        print "noon"
        ni_noon = ni
    if loop == (days-1)*dayseconds/dt + dayseconds/dt * 0.50:
        print "dusk"
        ni_dusk = ni
    if loop == (days-1)*dayseconds/dt + dayseconds/dt * 0.75:
        print "midnight"
        ni_midnight = ni
    if loop == (days-1)*dayseconds/dt + dayseconds/dt * 1.00:
        print "dawn again"
        ni_dawn = ni
    '''
    
    #*******************************************************************
    # Electron cooling by collisions                                   *
    #*******************************************************************
    
    # Rotation
    Lr_N2 = (3.5E-14*ni_sum*n['N2']*(Te - T)/np.power(Te,.5)) * (1.602E-19)*(1E-12)
    Lr_O2 = (5.2E-15*ni_sum*n['O2']*(Te - T)/np.power(Te,.5)) * (1.602E-19)*(1E-12)
    
    '''
    # N2 vibration
    first_sum[:] = 0
    second_sum[:] = 0
    for n in range(0,10):
        first_sum += (1 - np.exp((n+1)*E1*(1.0/Te - 1.0/T)))*np.power(10,(T93['A'][n] + T93['B'][n]*Te + T93['C'][n]*Te**2 + T93['D'][n]*Te**3 + T93['F'][n]*Te**4 - 16))
    for n in range(0,7):
        second_sum += (1 - np.exp((n+1)*E1*(1.0/Te - 1.0/T)))*np.power(10, (T95['A'][n] + T95['B'][n]*Te + T95['C'][n]*Te**2 + T95['D'][n]*Te**3 + T95['F'][n]*Te**4 - 16))
    
    Le_N2 = ((1e-12)*ni_sum*n['N2']*((1.0 - np.exp(-E1/T))*first_sum + (1.0 - np.exp(-E1/T))*np.exp(-E1/T)*second_sum)) * (1.602E-19)*(1E-12)
    '''
    
    # Vibration
    QTe = np.power(10, (-19.9171 + 0.0267*Te - (3.9960E-5)*np.power(Te,2.0) 
                  + (3.5187E-8)*np.power(Te,3.0) - (1.9228E-11)*np.power(Te,4.0) 
                  + (6.6865E-15)*np.power(Te,5.0) - (1.4791E-18)*np.power(Te,6.0) 
                  + (2.0127E-22)*np.power(Te,7.0) - (1.5346E-26)*np.power(Te,8.0) 
                  + (5.0148E-31)*np.power(Te,8.0)))
    Lv_O2 = (ni_sum*n['O2']*QTe*(1 - np.exp(2239.0*(np.power(Te,-1) - np.power(T,-1))))) * (1.602E-19)*(1E-12)
    
    # O Fine Structure
    D = 5 + np.exp(-326.6*np.power(T,-1)) + 3*np.exp(-227.7*np.power(T,-1))
    S21 = 1.863E-11
    S20 = 1.191E-11
    S10 = 8.249E-16*np.power(Te,0.6)*np.exp(-227.7*np.power(T,-1))
    Lf_O = (ni_sum*n['O']*np.power(D,-1)*(S10*(1 - np.exp(98.9*(np.power(Te,-1) - np.power(T,-1))) 
                                              + S20*(1 - np.exp(326.6*(np.power(Te,-1) - np.power(T,-1)))) 
                                              + S21*(1 - np.exp(227.7*(np.power(Te,-1) - np.power(T,-1))))))) * (1.602E-19)*(1E-12)
    
    #*******************************************************************
    # Resonant Ion-neutral interactions                                *
    #*******************************************************************
    
    v_OpO = 3.67E-11*n['O']*np.power(T,.5)*np.power(1-0.064*np.log10(T),2) * (1E-6) * (1E-3)
    v_O2pO2 = 2.59E-11*n['O2']*np.power(T,.5)*np.power(1-0.073*np.log10(T),2) * (1E-6) * (1E-3)
    v_N2pN2 = 5.14E-11*n['N2']*np.power(T,.5)*np.power(1-0.069*np.log10(T),2) * (1E-6) * (1E-3)
    
    Q_OpO = -(ni['O+']*m['O']*v_OpO)/(m['O'] + m['O'])*(3*k*(Ti - T))
    Q_O2pO2 = -(ni['O2+']*m['O2']*v_O2pO2)/(m['O2'] + m['O2'])*(3*k*(Ti - T))
    Q_N2pN2 = -(ni['N2+']*m['N2']*v_N2pN2)/(m['N2'] + m['N2'])*(3*k*(Ti - T))
    
    Q_r = Q_OpO + Q_O2pO2 + Q_N2pN2
    
    #*******************************************************************
    # Non-Resonant Ion-neutral interactions                            *
    #*******************************************************************
    
    v_OpN = 4.62*n['N'] * (1E-6) * (1E-10) * (1E-3)
    v_OpO2 = 6.64*n['O2'] * (1E-6) * (1E-10) * (1E-3)
    v_OpN2 = 6.82*n['N2'] * (1E-6) * (1E-10) * (1E-3)
    
    v_NpO = 4.42*n['O'] * (1E-6) * (1E-10) * (1E-3)
    v_NpO2 = 7.25*n['O2'] * (1E-6) * (1E-10) * (1E-3)
    v_NpN2 = 7.47*n['N2'] * (1E-6) * (1E-10) * (1E-3)
    
    v_N2pN = 2.95*n['N'] * (1E-6) * (1E-10) * (1E-3)
    v_N2pO = 2.58*n['O'] * (1E-6) * (1E-10) * (1E-3)
    v_N2pO2 = 4.49*n['O2'] * (1E-6) * (1E-10) * (1E-3)
    
    v_O2pO = 2.64*n['O'] * (1E-6) * (1E-10) * (1E-3)
    v_O2pN = 2.31*n['N'] * (1E-6) * (1E-10) * (1E-3)
    v_O2pN2 = 4.13*n['N2'] * (1E-6) * (1E-10) * (1E-4)
    
    Q_OpN = -(ni['O+']*m['O']*v_OpN)/(m['O'] + m['N'])*(3*k*(Ti - T))
    Q_OpO2 = -(ni['O+']*m['O']*v_OpO2)/(m['O'] + m['O2'])*(3*k*(Ti - T))
    Q_OpN2 = -(ni['O+']*m['O']*v_OpN2)/(m['O'] + m['N2'])*(3*k*(Ti - T))
    
    Q_NpO = -(ni['N+']*m['N']*v_NpO)/(m['N'] + m['O'])*(3*k*(Ti - T))
    Q_NpO2 = -(ni['N+']*m['N']*v_NpO2)/(m['N'] + m['O2'])*(3*k*(Ti - T))
    Q_NpN2 = -(ni['N+']*m['N']*v_NpN2)/(m['N'] + m['N2'])*(3*k*(Ti - T))
    
    Q_N2pN = -(ni['N2+']*m['N2']*v_N2pN)/(m['N2'] + m['N'])*(3*k*(Ti - T))
    Q_N2pO = -(ni['N2+']*m['N2']*v_N2pO)/(m['N2'] + m['O'])*(3*k*(Ti - T))
    Q_N2pO2 = -(ni['N2+']*m['N2']*v_N2pO2)/(m['N2'] + m['O2'])*(3*k*(Ti - T))  
    
    Q_O2pO = -(ni['O2+']*m['O2']*v_O2pO)/(m['O2'] + m['O'])*(3*k*(Ti - T))
    Q_O2pN = -(ni['O2+']*m['O2']*v_O2pN)/(m['O2'] + m['N'])*(3*k*(Ti - T))
    Q_O2pN2 = -(ni['O2+']*m['O2']*v_O2pN2)/(m['O2'] + m['N2'])*(3*k*(Ti - T))
    
    Q_nr = Q_OpN + Q_OpO2 + Q_OpN2 + Q_NpO + Q_NpO2 + Q_NpN2 + Q_N2pN + Q_N2pO + Q_N2pO2 + Q_O2pO + Q_O2pN + Q_O2pN2
    
    #*******************************************************************
    # Electron-Ion interactions                                        *
    #*******************************************************************    
    v_ei_Op = 54.5*ni['O+']/np.power(Te, 3.0/2.0) * (1E-6)
    v_ei_Np = 54.5*ni['N+']/np.power(Te, 3.0/2.0) * (1E-6)
    v_ei_O2p = 54.5*ni['O2+']/np.power(Te, 3.0/2.0) * (1E-6)
    v_ei_N2p = 54.5*ni['N2+']/np.power(Te, 3.0/2.0) * (1E-6)
    
    Q_ei_Op = -(ni_sum*me*v_ei_Op)/(m['O'] + me)*(3*k*(Te - Ti))
    Q_ei_Np = -(ni_sum*me*v_ei_Np)/(m['N'] + me)*(3*k*(Te - Ti))
    Q_ei_O2p = -(ni_sum*me*v_ei_O2p)/(m['O2'] + me)*(3*k*(Te - Ti))
    Q_ei_N2p = -(ni_sum*me*v_ei_N2p)/(m['N2'] + me)*(3*k*(Te - Ti))

    Q_ei = Q_ei_Op + Q_ei_Np + Q_ei_O2p + Q_ei_N2p
    
    #*******************************************************************
    # Electron-Neutral interactions                                    *
    #*******************************************************************     
    
    v_en_O = 8.9E-11*n['O']*(1 + 5.7E-4*Te)*np.power(Te,.5) * (1E-6)
    v_en_O2 = 1.82E-10*n['O2']*(1 + 3.6E-2*np.power(Te,.5))*np.power(Te,.5) * (1E-6)
    v_en_N2 = 2.33E-11*n['N2']*(1 - 1.21E-4*Te)*Te * (1E-6)
    
    Q_en_O = -(ni_sum*me*v_en_O)/(m['O'] + me)*(3*k*(Te - T))
    Q_en_O2 = -(ni_sum*me*v_en_O2)/(m['O2'] + me)*(3*k*(Te - T))
    Q_en_N2 = -(ni_sum*me*v_en_O)/(m['N2'] + me)*(3*k*(Te - T))
    
    Q_en = Q_en_O + Q_en_O2 + Q_en_N2
    
    




'''
ni_sum_dawn = ni_dawn['O+'] + ni_dawn['N+'] + ni_dawn['O2+'] + ni_dawn['N2+']
ni_sum_noon = ni_noon['O+'] + ni_noon['N+'] + ni_noon['O2+'] + ni_noon['N2+']
ni_sum_dusk = ni_dusk['O+'] + ni_dusk['N+'] + ni_dusk['O2+'] + ni_dusk['N2+']
ni_sum_midnight = ni_midnight['O+'] + ni_midnight['N+'] + ni_midnight['O2+'] + ni_midnight['N2+']
'''

#***********************************************************************
# Time                                                                 *
#***********************************************************************

#print "We took %f seconds!" % (time.time()-start)


#plt.show()

#***********************************************************************
# Shelve dat shit!!                                                    *
#***********************************************************************
'''
filename='/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/shelve.out'
my_shelf = shelve.open(filename,'n') # 'n' for new

for key in dir():
    try:
        my_shelf[key] = globals()[key]
    except TypeError:
        #
        # __builtins__, my_shelf, and imported modules can not be shelved.
        #
        print('ERROR shelving: {0}'.format(key))
my_shelf.close()

'''
'''
filename='/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/shelve.out'
my_shelf = shelve.open(filename)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

'''


np.savetxt("Te.txt", Te)
np.savetxt("Ti.txt", Ti)
np.savetxt("T.txt", T)
np.savetxt("Q_ei.txt", Q_ei)
np.savetxt("Q_nr.txt", Q_nr)
np.savetxt("Q_r.txt", Q_r)
np.savetxt("Lr_N2.txt", Lr_N2)
np.savetxt("Lr_O2.txt", Lr_O2)
np.savetxt("Lf_O.txt", Lf_O)
np.savetxt("Lv_O2.txt", Lv_O2)
np.savetxt("z.txt",z)





#***********************************************************************
# Plotting the results.                                                *
#***********************************************************************

# Temperature vs Altitude
fig1 = plt.figure(figsize=[12.15,6.06])
ax1 = fig1.add_subplot(111)
ax1.set_title('Temperature vs Altitude')
ax1.set_xlabel('Temperature $[K]$')
ax1.set_ylabel('Altitude $[km]$')
ax1.set_xlim([200,1900])
ax1.plot(T, z*1E-3)
fig1.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Temperature_vs_Altitude.eps')

# Scale Height vs Altitude
fig2 = plt.figure(figsize=[12.15,6.06])
ax2 = fig2.add_subplot(111)
ax2.set_title('Scale Height vs Altitude')
ax2.set_xlabel('Scale Height $[m]$')
ax2.set_ylabel('Altitude $[km]$')
ax2.plot(H['O'], z*1E-3, H['N'], z*1E-3, H['O2'], z*1E-3,\
 H['N2'], z*1E-3, H['NO'], z*1E-3)
ax2.legend(['O','N', 'O2', 'N2', 'NO'])
fig2.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Scale_Height_vs_Altitude.eps')

# Number Density vs Altitude 100km - 800km
fig3 = plt.figure(figsize=[12.15,6.06])
ax3 = fig3.add_subplot(111)
ax3.set_title('Number Density vs Altitude')
ax3.set_xlabel('Number Density $[m^{-3}]$')
ax3.set_ylabel('Altitude $[km]$')
ax3.set_ylim([0,800])
ax3.semilogx(n['O'], z*1E-3, n['N'], z*1E-3, n['O2'], z*1E-3,\
 n['N2'], z*1E-3, n['NO'], z*1E-3)
ax3.legend(['O','N', 'O2', 'N2', 'NO'])
fig3.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Number_Density_vs_Altitude_100_800.eps')
'''
# Number Density vs Altitude 400km - 500km
fig4 = plt.figure(figsize=[12.15,6.06])
ax4 = fig4.add_subplot(111)
ax4.set_title('Number Density vs Altitude')
ax4.set_xlabel('Number Density $[m^{-3}]$')
ax4.set_ylabel('Altitude $[km]$')
ax4.set_xlim([1E14, 1E15])
ax4.set_ylim([400,500])
ax4.semilogx(n['O'], z*1E-3, n['N2'], z*1E-3)
ax4.legend(['O','N', 'O2', 'N2', 'NO'])
fig4.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Number_Density_vs_Altitude_400_500.eps')
'''
# Mass Density vs Altitude
fig5 = plt.figure(figsize=[12.15,6.06])
ax5 = fig5.add_subplot(111)
ax5.set_title('Mass Density vs Altitude')
ax5.set_xlabel('Mass Density $[\\frac{km}{m^{3}}]$')
ax5.set_ylabel('Altitude $[km]$')
ax5.set_ylim([0,800])
ax5.semilogx(rho['O'], z*1E-3, rho['N'], z*1E-3, rho['O2'], z*1E-3,\
 rho['N2'], z*1E-3, rho['NO'], z*1E-3)
ax5.legend(['O','N', 'O2', 'N2', 'NO'])
fig5.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Mass_Density_vs_Altitude.eps')

# Total Weighted Molecular Mass vs Altitude
fig6 = plt.figure(figsize=[12.15,6.06])
ax6 = fig6.add_subplot(111)
ax6.set_title('Total Weighted Molecular Mass vs Altitude')
ax6.set_xlabel('Total Weighted Molecular Mass [amu]')
ax6.set_ylabel('Altitude $[km]$')
ax6.set_ylim([0,800])
ax6.plot(M_molec, z*1E-3)
fig6.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Total_Weighted_Molecular_Mass_vs_Altitude.eps')

# Column Number Density vs Altitude 100km - 800km
fig7 = plt.figure(figsize=[12.15,6.06])
ax7 = fig7.add_subplot(111)
ax7.set_title('Column Number Density vs Altitude')
ax7.set_xlabel('Column Number Density $[m^{-2}]$')
ax7.set_ylabel('Altitude $[km]$')
ax7.set_ylim([0,800])
ax7.semilogx(col_n['O'], z*1E-3, col_n['N'], z*1E-3,\
 col_n['O2'], z*1E-3, col_n['N2'], z*1E-3, col_n['NO'], z*1E-3)
ax7.legend(['O','N', 'O2', 'N2', 'NO'])
fig7.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Column_Number_Density_vs_Altitude.eps')

# Column Mass Density vs Altitude 100km - 800km
fig8 = plt.figure(figsize=[12.15,6.06])
ax8 = fig8.add_subplot(111)
ax8.set_title('Column Mass Density vs Altitude')
ax8.set_xlabel('Column Mass Density $[\\frac{kg}{m^{2}}]$')
ax8.set_ylabel('Altitude $[km]$')
ax8.set_ylim([0,800])
ax8.semilogx(col_rho['O'], z*1E-3, col_rho['N'], z*1E-3,\
 col_rho['O2'], z*1E-3, col_rho['N2'], z*1E-3,\
 col_rho['NO'], z*1E-3)
ax8.legend(['O','N', 'O2', 'N2', 'NO'])
fig8.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Column_Mass_Density_vs_Altitude.eps')

# Tau vs Altitude 100km - 800km
fig9 = plt.figure(figsize=[12.15,6.06])
ax9 = fig9.add_subplot(111)
ax9.set_title('Optical Depth vs Altitude')
ax9.set_xlabel('Optical Depth')
ax9.set_ylabel('Altitude $[km]$')
ax9.set_xlim([1E-6, 1E4])
ax9.set_ylim([0,600])
ax9.semilogx(sumtau, z*1E-3)
fig9.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Optical_Depth_vs_Altitude.eps')

# Flux vs Altitude 100km - 800km
fig10 = plt.figure(figsize=[12.15,6.06])
ax10 = fig10.add_subplot(111)
ax10.set_title('Flux vs Altitude')
ax10.set_xlabel('Flux $[\\frac{photons}{m2*s}]$')
ax10.set_ylabel('Altitude $[km]$')
ax10.set_ylim([0,600])
ax10.plot(I, z*1E-3)
fig10.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Flux_vs_Altitude.eps')

# Energy [W/m3] vs Altitude [km] [100km - 800km]
fig11 = plt.figure(figsize=[12.15,6.06])
ax11 = fig11.add_subplot(111)
ax11.set_title('Energy Dissipation vs Altitude')
ax11.set_xlabel('Energy Dissipation $[\\frac{W}{m3}]$')
ax11.set_ylabel('Altitude $[km]$')
ax11.set_xlim([1E-11, 1E-7])
ax11.set_ylim([0,400])
ax11.semilogx(sumQ, z*1E-3)
fig11.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/Energy_Dissipation_vs_Altitude.eps')

# Maximum Temperature vs Time
fig12 = plt.figure(figsize=[12.15,6.06])
ax12 = fig12.add_subplot(111)
ax12.set_title('Maximum Temperature vs Time')
ax12.set_xlabel('Time $[s]$')
ax12.set_ylabel('Temperature $[K]$')
ax12.plot(time, T_max)
fig12.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/A/T_max_vs_time.eps')

# Ion Production Rate vs Altitude
fig13 = plt.figure(figsize=[12.15,6.06])
ax13 = fig13.add_subplot(111)
ax13.set_title('Ion Production Rate vs Altitude')
ax13.set_xlabel('Ion Production Rate $[\\frac{ions}{s \cdot m^{3}}]$')
ax13.set_ylabel('Altitude $[km]$')
ax13.set_xlim([1E-1, 1E11])
ax13.set_ylim([0,800])
ax13.semilogx(Ps['O+'], z*1E-3, Ps['N+'], z*1E-3, Ps['O2+'], z*1E-3,\
 Ps['N2+'], z*1E-3)
ax13.legend(['O+','N+', 'O2+', 'N2+'])
fig13.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Production_Rate_vs_Altitude_100_800.eps')

# Ion Number Density vs Altitude
fig14 = plt.figure(figsize=[12.15,6.06])
ax14 = fig14.add_subplot(111)
ax14.set_title('Ion Number Density vs Altitude (Noon)')
ax14.set_xlabel('Ion Number Density $[m^{-3}]$')
ax14.set_ylabel('Altitude $[km]$')
ax14.set_xlim([1E-4, 1E14])
ax14.set_ylim([0,800])
ax14.semilogx(ni['O+'], z*1E-3, ni['N+'], z*1E-3, ni['O2+'], z*1E-3,\
 ni['N2+'], z*1E-3)#, Oplus_old, z*1E-3)
ax14.legend(['O+','N+', 'O2+', 'N2+', 'O+$_{old}$'])
fig14.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Number_Density_vs_Altitude.eps')
'''
# Ion Number Density vs Altitude (Dawn)
fig15 = plt.figure(figsize=[12.15,6.06])
ax15 = fig15.add_subplot(111)
ax15.set_title('Ion Number Density vs Altitude (Dawn)')
ax15.set_xlabel('Ion Number Density $[m^{-3}]$')
ax15.set_ylabel('Altitude $[km]$')
ax15.set_xlim([1E-4, 1E14])
ax15.set_ylim([0,800])
ax15.semilogx(ni_dawn['O+'], z*1E-3, ni_dawn['N+'], z*1E-3, ni_dawn['O2+'], z*1E-3,\
 ni_dawn['N2+'], z*1E-3)
ax15.legend(['O+','N+', 'O2+', 'N2+', 'O+$_{old}$'])
fig15.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Number_Density_vs_Altitude_Dawn.eps')

# Ion Number Density vs Altitude (Noon)
fig16 = plt.figure(figsize=[12.15,6.06])
ax16 = fig16.add_subplot(111)
ax16.set_title('Ion Number Density vs Altitude (Noon)')
ax16.set_xlabel('Ion Number Density $[m^{-3}]$')
ax16.set_ylabel('Altitude $[km]$')
ax16.set_xlim([1E-4, 1E14])
ax16.set_ylim([0,800])
ax16.semilogx(ni_noon['O+'], z*1E-3, ni_noon['N+'], z*1E-3, ni_noon['O2+'], z*1E-3,\
 ni_noon['N2+'], z*1E-3)
ax16.legend(['O+','N+', 'O2+', 'N2+', 'O+$_{old}$'])
fig16.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Number_Density_vs_Altitude_Noon.eps')

# Ion Number Density vs Altitude (Dusk)
fig17 = plt.figure(figsize=[12.15,6.06])
ax17 = fig17.add_subplot(111)
ax17.set_title('Ion Number Density vs Altitude (Dusk)')
ax17.set_xlabel('Ion Number Density $[m^{-3}]$')
ax17.set_ylabel('Altitude $[km]$')
ax17.set_xlim([1E-4, 1E14])
ax17.set_ylim([0,800])
ax17.semilogx(ni_dusk['O+'], z*1E-3, ni_dusk['N+'], z*1E-3, ni_dusk['O2+'], z*1E-3,\
 ni_dusk['N2+'], z*1E-3)
ax17.legend(['O+','N+', 'O2+', 'N2+', 'O+$_{old}$'])
fig17.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Number_Density_vs_Altitude_Dusk.eps')

# Ion Number Density vs Altitude (Midnight)
fig18 = plt.figure(figsize=[12.15,6.06])
ax18 = fig18.add_subplot(111)
ax18.set_title('Ion Number Density vs Altitude (Midnight)')
ax18.set_xlabel('Ion Number Density $[m^{-3}]$')
ax18.set_ylabel('Altitude $[km]$')
ax18.set_xlim([1E-4, 1E14])
ax18.set_ylim([0,800])
ax18.semilogx(ni_midnight['O+'], z*1E-3, ni_midnight['N+'], z*1E-3, ni_midnight['O2+'], z*1E-3,\
 ni_midnight['N2+'], z*1E-3)
ax18.legend(['O+','N+', 'O2+', 'N2+', 'O+$_{old}$'])
fig18.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Number_Density_vs_Altitude_Midnight.eps')
'''
# Ion Number Density vs Time
fig19 = plt.figure(figsize=[12.15,6.06])
ax19 = fig19.add_subplot(111)
ax19.set_title('Ion Number Density vs Time')
ax19.set_xlabel('Time $[s]$')
ax19.set_ylabel('Ion Number Density $[m^{-3}]$')
ax19.plot(time, ni_max['O+'], time, ni_max['N+'], time, ni_max['O2+'],\
 time, ni_max['N2+'])
ax19.legend(['O+','N+', 'O2+', 'N2+'])
fig19.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Number_Density_vs_time.eps')

# Total Ion Number Density vs Altitude
fig20 = plt.figure(figsize=[12.15,6.06])
ax20 = fig20.add_subplot(111)
ax20.set_title('Total Ion Number Density vs Altitude (Noon)')
ax20.set_xlabel('Total Ion Number Density $[m^{-3}]$')
ax20.set_ylabel('Altitude $[km]$')
#ax20.set_xlim([1E9, 1E13])
ax20.set_xlim([1E4, 1E13])
ax20.set_ylim([0,800])
ax20.semilogx(ni_sum, z*1E-3)
fig20.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Total_Ion_Number_Density_vs_Altitude.eps')
'''
# Total Ion Number Density vs Altitude (Dawn)
fig21 = plt.figure(figsize=[12.15,6.06])
ax21 = fig21.add_subplot(111)
ax21.set_title('Total Ion Number Density vs Altitude (Dawn)')
ax21.set_xlabel('Total Ion Number Density $[m^{-3}]$')
ax21.set_ylabel('Altitude $[km]$')
#ax21.set_xlim([1E9, 1E13])
ax21.set_xlim([1E4, 1E13])
ax21.set_ylim([0,800])
ax21.semilogx(ni_sum_dawn, z*1E-3)
fig21.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Total_Ion_Number_Density_vs_Altitude_dawn.eps')

# Total Ion Number Density vs Altitude (Noon)
fig22 = plt.figure(figsize=[12.15,6.06])
ax22 = fig22.add_subplot(111)
ax22.set_title('Total Ion Number Density vs Altitude (Noon)')
ax22.set_xlabel('Total Ion Number Density $[m^{-3}]$')
ax22.set_ylabel('Altitude $[km]$')
#ax22.set_xlim([1E9, 1E13])
ax22.set_xlim([1E4, 1E13])
ax22.set_ylim([0,800])
ax22.semilogx(ni_sum_noon, z*1E-3)
fig22.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Total_Ion_Number_Density_vs_Altitude_noon.eps')

# Total Ion Number Density vs Altitude (Dusk)
fig23 = plt.figure(figsize=[12.15,6.06])
ax23 = fig23.add_subplot(111)
ax23.set_title('Total Ion Number Density vs Altitude (Dusk)')
ax23.set_xlabel('Total Ion Number Density $[m^{-3}]$')
ax23.set_ylabel('Altitude $[km]$')
#ax23.set_xlim([1E9, 1E13])
ax23.set_xlim([1E4, 1E13])
ax23.set_ylim([0,800])
ax23.semilogx(ni_sum_dusk, z*1E-3)
fig23.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Total_Ion_Number_Density_vs_Altitude_dusk.eps')

# Total Ion Number Density vs Altitude (Midnight)
fig24 = plt.figure(figsize=[12.15,6.06])
ax24 = fig24.add_subplot(111)
ax24.set_title('Total Ion Number Density vs Altitude (Midnight)')
ax24.set_xlabel('Total Ion Number Density $[m^{-3}]$')
ax24.set_ylabel('Altitude $[km]$')
#ax24.set_xlim([1E9, 1E13])
ax24.set_xlim([1E4, 1E13])
ax24.set_ylim([0,800])
ax24.semilogx(ni_sum_midnight, z*1E-3)
fig24.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Total_Ion_Number_Density_vs_Altitude_midnight.eps')
'''
# Maximum Temperature vs Time
fig25 = plt.figure(figsize=[12.15,6.06])
ax25 = fig25.add_subplot(111)
ax25.set_title('Maximum Temperature vs Time')
ax25.set_xlabel('Time $[days]$')
ax25.set_ylabel('Temperature $[K]$')
ax25.set_ylim([0, 7000])
ax25.plot(time/dayseconds, T_max, time/dayseconds, Te_max, time/dayseconds, Ti_max)
ax25.legend(['T$_n$','T$_e$', 'T$_i$'])
fig25.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/T_max_vs_time.eps')

# F2 Peak vs Time
fig26 = plt.figure(figsize=[12.15,6.06])
ax26 = fig26.add_subplot(111)
ax26.set_title('F2 Peak Density vs Time')
ax26.set_xlabel('Time $[days]$')
ax26.set_ylabel('F2 Peak Density $[m^{-3}]$')
ax26.plot(time/dayseconds, F2_peak)
fig26.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/F2_max_vs_time.eps')

# Electron Cooling Rates vs Altitude
fig27 = plt.figure(figsize=[12.15,6.06])
ax27 = fig27.add_subplot(111)
ax27.set_title('Electron Cooling Rate vs Altitude')
ax27.set_xlabel('Electron Cooling Rate $[m^{-3}]$')
ax27.set_ylabel('Altitude $[km]$')
ax27.set_ylim([0,800])
ax27.set_xlim([1E-29,1E-13])
ax27.semilogx(Lr_N2, z*1E-3, Lr_O2, z*1E-3, Lf_O, z*1E-3, Lv_O2, z*1E-3)
ax27.legend(['L$_{N2}$ Rotational','L$_{O2}$ Rotational','L$_{O}$ Fine Structure', 'L$_{O2}$ Vibrational'])
fig27.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Electron_Cooling_Rates_vs_Altitude.eps')

# Ion Interactions vs Altitude
fig27 = plt.figure(figsize=[12.15,6.06])
ax27 = fig27.add_subplot(111)
ax27.set_title('Ion Interactions vs Altitude')
ax27.set_xlabel('Ion Interactions $[m^{-3}]$')
ax27.set_ylabel('Altitude $[km]$')
ax27.set_ylim([0,800])
ax27.semilogx(np.abs(Q_ei), z*1E-3, np.abs(Q_nr), z*1E-3, np.abs(Q_r), z*1E-3)
ax27.legend(['Q$_{ei}$','Q$_{nr}$','Q$_{r}$'])
fig27.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Neutral_Interactions_vs_Altitude.eps')

'''
# Ion Interactions vs Altitude
fig27 = plt.figure(figsize=[12.15,6.06])
ax27 = fig27.add_subplot(111)
ax27.set_title('Ion Interactions vs Altitude')
ax27.set_xlabel('Ion Interactions $[m^{-3}]$')
ax27.set_ylabel('Altitude $[km]$')
ax27.set_ylim([0,800])
ax27.plot(np.abs(Q_ei), z*1E-3, np.abs(Q_nr), z*1E-3, np.abs(Q_r), z*1E-3)
ax27.legend(['Q$_{ei}$','Q$_{nr}$','Q$_{r}$'])
fig27.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 5/Figures/B/Ion_Neutral_Interactions_vs_Altitude_normalabs.eps')
'''
