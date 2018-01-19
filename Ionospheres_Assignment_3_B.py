# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 20:55:59 2012

@author: jsni
"""

# Importing everything python needs in order to be smart.
import numpy as np                                      # Numpy tells python how to work with numbers.
import matplotlib.pyplot as plt                         # Matplotlib tells python how to make pretty plots.
import ionospheres_programs as ip                       # Importing my custom library.
import shelve                                           # For saving calculated variables/arrays to an external file so that you don't have to keep running your code over again becaue it takes so darn long!

#***********************************************************************
# Unshelving all of the values calculated previously (which are        *
# presumably the equilibrium values)                                   *
#***********************************************************************

filename='/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 3/shelveA.out'
my_shelf = shelve.open(filename)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

#***********************************************************************
# Declaring some variables.                                            *
#***********************************************************************

dayseconds = 86400.0                                                    # number of seconds in a day.
days = 1                                                                # Number of days you want to calculate over.
dt = 300.0                                                              # Definin time incriment as 300 [s].
tsize = dayseconds*days/dt                                              # Number of second in a day divided by 300s intervals gives us the number of "angular steps" used to calculate chi over the course of a day.
tloop = int(tsize)                                                      # Defining an integer range value for the time dependent loop.
time = np.arange(tsize)*dt                                              # Time

#***********************************************************************
# Building the ionosphere                                              *
#***********************************************************************

#Creating a list of all of the ion species tht live in our virtual
#thermosphere so that we can later assign properties to them.
ions = ['O+', 'N+', 'O2+', 'N2+', 'NO+']

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
    n = ip.calculate_number_densities(n, T, z, H, n_base)
    
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
    col_n = ip.calculate_column_densities(col_n, n, H, z)
    
    # Defining a column mass density just for fun becuase I want to  
    # see what the plot looks like and I feel like it is an important 
    # quantity in atmospheric science.
    col_rho = {'O':(col_n['O']*m['O']),
               'N':(col_n['N']*m['N']), 
               'O2':(col_n['O2']*m['O2']),
               'N2':(col_n['N2']*m['N2']), 
               'NO':(col_n['NO']*m['NO'])}
    
    # Fetching the optical depths from the custom library.
    tau = ip.calculate_optical_depth(tau, data, col_n, zsize)
    
    # Setting the new value of solar zenith angle over each iteration.
    X += dt/dayseconds*2*np.pi
    
    # Define the new value of the cosine of the solar zenith angle 
    # for each iteration.
    sza = np.cos(X)
    
    # Setting up a minimum value for the solar zenith angle which is 
    # equivalent to recognizing that there is still ionization on the 
    # night side due to a flux of galactic cosmic rays etc...  
    if sza <= 0.05:
        sza = 0.05
    
    # Fetching the combined optical depths from the custom library.    
    sumtau = ip.combine_optical_depths(sumtau, tau, sza)
    
    # Fetching the intensity from the custom library.
    I = ip.calculate_inensity(I, I_inf, sumtau, zsize)
    
    # Fetching the energy dissipation from the custom library.
    Q = ip.calculate_energy_dissipation(Q, eps, n, data, I, e, zsize)
    
    # Fetching the optical depths from the custom library.
    sumQ = ip.combine_energy_dissipation(sumQ, Q)
    
    # Defineing the thermal conduction in [W/(m*K)].
    lam_n = 3.6E-4*np.power(T, .75)

    #*******************************************************************
    # Building the rest of the  ionosphere                             *
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
    
    # Ion loss terms.
    RM = {'O+':(R47*M47 + R48*M48 + R49*R49),
            'N+':(R37*M37 + R38*M38 + R39*M39),
            'O2+':(R52*M52 + R53*M53),
            'N2+':(R43*M43 + R44*M44 + R45*M45 + R46*M46)}
    
    ni = ip.calculate_ion_number_densities(ni, source, RM, dt)
    
    ni_sum = ni['O+'] + ni['N+'] + ni['O2+'] + ni['N2+']
    
   
    
    for i in ions:
        ni_max[i][loop] = ni[i][zsize-1]
#    print ni_max
#    print( 'Percent completed: %f' % (100.0*loop/np.size(ni)) )




#***********************************************************************
# Shelve dat shit!!                                                    *
#***********************************************************************

filename='/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 3/shelveB.out'
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

#***********************************************************************
# Plotting all the results.                                            *
#***********************************************************************

# Ion Production Rate vs Altitude
fig1 = plt.figure(figsize=[12.15,6.06])
ax1 = fig1.add_subplot(111)
ax1.set_title('Ion Production Rate vs Altitude')
ax1.set_xlabel('Ion Production Rate $[\\frac{ions}{s \cdot m^{3}}]$')
ax1.set_ylabel('Altitude $[km]$')
ax1.set_xlim([1E-1, 1E11])
ax1.set_ylim([0,800])
ax1.semilogx(Ps['O+'], z*1E-3, Ps['N+'], z*1E-3, Ps['O2+'], z*1E-3,\
 Ps['N2+'], z*1E-3)
ax1.legend(['O+','N+', 'O2+', 'N2+'])
fig1.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 3/Figures/B/Ion_Production_Rate_vs_Altitude_100_800.eps')

# Ion Number Density vs Altitude 100km - 800km
fig2 = plt.figure(figsize=[12.15,6.06])
ax2 = fig2.add_subplot(111)
ax2.set_title('Ion Number Density vs Altitude')
ax2.set_xlabel('Ion Number Density $[m^{-3}]$')
ax2.set_ylabel('Altitude $[km]$')
ax2.set_xlim([1E-4, 1E14])
ax2.set_ylim([0,800])
ax2.semilogx(ni['O+'], z*1E-3, ni['N+'], z*1E-3, ni['O2+'], z*1E-3,\
 ni['N2+'], z*1E-3)
ax2.legend(['O+','N+', 'O2+', 'N2+'])
fig2.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 3/Figures/B/Ion_Number_Density_vs_Altitude_100_800.eps')

# Ion Number Density vs Time
fig3 = plt.figure(figsize=[12.15,6.06])
ax3 = fig3.add_subplot(111)
ax3.set_title('Ion Number Density vs Time')
ax3.set_xlabel('Time $[s]$')
ax3.set_ylabel('Ion Number Density $[m^{-3}]$')
ax3.plot(time, ni_max['O+'], time, ni_max['N+'], time, ni_max['O2+'],\
 time, ni_max['N2+'])
ax3.legend(['O+','N+', 'O2+', 'N2+'])
fig3.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 3/Figures/B/Ion_Number_Density_vs_time.eps')

# Total Ion Number Density vs Altitude 100km - 800km
fig4 = plt.figure(figsize=[12.15,6.06])
ax4 = fig4.add_subplot(111)
ax4.set_title('Total Ion Number Density vs Altitude')
ax4.set_xlabel('Total Ion Number Density $[m^{-3}]$')
ax4.set_ylabel('Altitude $[km]$')
#ax4.set_xlim([1E9, 1E13])
ax4.set_xlim([1E4, 1E13])
ax4.set_ylim([0,600])
ax4.semilogx(ni_sum, z*1E-3)
fig4.savefig('/Users/jsni/Desktop/Ionospheres Assignments/\
Assignment 3/Figures/B/Total_Ion_Number_Density_vs_Altitude_100_800.eps')