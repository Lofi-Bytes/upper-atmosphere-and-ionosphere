# -*- coding: utf-8 -*-

"""
Last Edited on 01/26/2012
@author: Jonathan Nickerson
"""

#Importing everything python needs in order to be smart
import numpy as np                                      #Numpy tells python how to work with numbers
import matplotlib.pyplot as plt                         #Matplotlib tells python how to make pretty plots
import ionospheres_programs as ip                       #Custom library
import time

start = time.time()

#Defining a Bunch of Variables
theta = np.arange(1000.0)*np.pi/999.0                   #Theta ranges from 0 [rad] - pi [rad] in incriments of 1000. 
T = 800.0*np.tanh(theta*3.0)+200.0                      #Temperature ranges from 200-1000 [K] in increments of 1000.
z = (np.arange(1000.0)/999.0*700.0+100.0)*1E3           #Altitude ranges from 100,000 [m] - 800,000 [m] in increments of 1000.
delta_z = (z[999] - z[0])/999                           #Defining a delta z in [m].
R = 6378.1*1E3                                          #Radius of Earth in [m].
g0 = 9.81                                               #Gravity at Earths surface in [m/s^2].
g = g0*(R/(R+z))*(R/(R+z))                              #Gravity as it varies with altitude in [m/s^2].
k = 1.3806503E-23                                       #Boltzmann constant in [((m^2)*kg)/((s^2)*K)].
A = 6.022E23                                            #Avagadro's Number.
planck = 6.626068E-34                                   #Planck's Constant in [((m^2)*kg)/s].
c = 3.0E8                                               #Speed of Light in [m/s].
P = 150.0                                               #F10.7 proxy in [SFU].
eps = 0.25                                              #Heating efficiency.

#Creating a list of all of the neutral species that live in our virtual 
#thermosphere so that we can later assign properties to them.
species = ['O','N','O2','N2','NO']
    
#Creating a dictionary of masses for each species. 
#Divide by (1000.0*A) to convert from [amu] to [kg].
m = {'O':(16.0/(1000.0*A)),'N':(14.0/(1000.0*A)),'O2':(32.0/(1000.0*A)),
     'N2':(28.0/(1000.0*A)),'NO':(30.0/(1000.0*A))}

#Creating a dictionary of initial number density values at the base of the 
#thermospehre. Densities of Each Species are given in [m^-3] at 100km.
n_base = {'O':4.9E17,'N':6.4E11,'O2':2.4E18,'N2':1.0E19,'NO':1.0E14}

#Calculating the average mass density at the base of the thermosphere in 
#[kg/m^3].
rho_tot = n_base['O']*m['O'] + n_base['N']*m['N'] + n_base['O2']*m['O2']\
+ n_base['N2']*m['N2'] + n_base['NO']*m['NO']

#Calling a definition inside the custom library which fetches the solar flux 
#data from Aaron's poorly formatted table.
data = ip.read_solar_flux()

#Using the data fetched to define some more variables
lamb = ((data['lower'] + data['upper'])/2.0)*1E-10      #Use the Average Solar Spectrum angstroms to meters
e = (planck*c)/lamb                                     #Photon energy
data['F74113'] *= 1E9*1E4                               #Convert 1E9 photons/cm^2*s to photons/m^2*s
I_inf = data['F74113']*(1.0 + data['Abi']*(P-80.0))     #Solar irradiance at infinity

#~~~~~~~~~~~~~~~~~ Begin temperature dependent loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~

for loop in range(10):
    #Creating a dictionary of scale height's for each species in [m].
    H = {'O':((k*T)/(m['O']*g)),'N':((k*T)/(m['N']*g)),'O2':((k*T)/(m['O2']*g)),
         'N2':((k*T)/(m['N2']*g)),'NO':((k*T)/(m['NO']*g))}    
    #Fetching the average scale heights from the custom libaray.
    H_ave = ip.calculate_average_scale_height(n_base, m, H, rho_tot)
    #Fetching the number densities from the custom libaray.
    #n = ip.calculate_number_densities(T, z, n_base, m, k, g)    
    n = ip.calculate_number_densities(T, z, H, n_base)
    #Creating a dictionary to store the mass densities of each species as they vary
    #with altitude in [kg/m^3].
    rho = {'O':(n['O']*m['O']), 'N':(n['N']*m['N']), 'O2':(n['O2']*m['O2']), 
           'N2':(n['N2']*m['N2']), 'NO':(n['NO']*m['NO'])}
    #Calculating and array containing the total number density as it changes with 
    #altitude in [m^-3]. 
    N_tot = n['O'] + n['N'] + n['O2'] + n['N2'] + n['NO']
    #Calculating the total weighted mass as it changes with altitude in [kg].
    M = m['O']*(n['O']/N_tot) + m['N']*(n['N']/N_tot) + m['O2']*(n['O2']/N_tot)\
    + m['N2']*(n['N2']/N_tot) + m['NO']*(n['NO']/N_tot)
    #Total weighted molecular mass as it changes with altitude in [amu].
    M_molec = M*(1000.0*A)
    #Fetching the column densities from the custom library.
    col_n = ip.calculate_column_densities(n, H, z)
    #Defining a column mass density just for fun becuase I want to see what the 
    #plot looks like and I feel like it is an important quantity in atmospheric 
    #science.
    col_rho = {'O':(col_n['O']*m['O']), 'N':(col_n['N']*m['N']), 
               'O2':(col_n['O2']*m['O2']), 'N2':(col_n['N2']*m['N2']), 
    'NO':(col_n['NO']*m['NO'])}
    #Fetching the optical depths from the custom library.
    tau = ip.calculate_optical_depth(data, col_n)
    #Fetching the combined optical depths from the custom library.
    sumtau = ip.combine_optical_depths(tau)
    #Fetching the intensity from the custom library.
    I = ip.calculate_inensity(I_inf, sumtau)
    #Fetching the energy dissipation from the custom library.
    Q = ip.calculate_energy_dissipation(eps, n, data, I, e)
    #Fetching the optical depths from the custom library.
    sumQ = ip.combine_energy_dissipation(Q)
    #Defineing the thermal conduction in [W/(m*K)].
    lam_n = 3.6E-4*pow(T, .75)    
    #Calling the new temperature profil from the custom library and recalculating 
    #everything again.
    T = ip.T_new(T, sumQ, delta_z, lam_n)

#plt.show()

#~~~~~~~~~~~~~~~~~ End temperature dependent loop ~~~~~~~~~~~~~~~~~~~~~~~~~~~

#Calculating the number density one last time with our new temperature outside the 
#loop so that we can plot everything.
H = {'O':((k*T)/(m['O']*g)),'N':((k*T)/(m['N']*g)),'O2':((k*T)/(m['O2']*g)),
     'N2':((k*T)/(m['N2']*g)),'NO':((k*T)/(m['NO']*g))}    
#Fetching the average scale heights from the custom libaray.
H_ave = ip.calculate_average_scale_height(n_base, m, H, rho_tot)
#Fetching the number densities from the custom libaray.
#n = ip.calculate_number_densities(T, z, n_base, m, k, g)
n = ip.calculate_number_densities(T, z, H, n_base)                

#************************* PLOTS **********************************************

#Temperature vs Altitude
fig1 = plt.figure(figsize=[12.15,6.06])
ax1 = fig1.add_subplot(111)
ax1.set_title('Temperature vs Altitude')
ax1.set_xlabel('Temperature [K]')
ax1.set_ylabel('Altitude [km]')
ax1.set_xlim([200,1900])
ax1.plot(T, z*1E-3)
fig1.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Temperature_vs_Altitude.eps')

#Scale Height vs Altitude
fig2 = plt.figure(figsize=[12.15,6.06])
ax2 = fig2.add_subplot(111)
ax2.set_title('Scale Height vs Altitude')
ax2.set_xlabel('Scale Height [m]')
ax2.set_ylabel('Altitude [km]')
ax2.plot(H['O'], z*1E-3, H['N'], z*1E-3, H['O2'], z*1E-3, H['N2'], z*1E-3,\
 H['NO'], z*1E-3)
ax2.legend(['O','N', 'O2', 'N2', 'NO'])
fig2.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Scale_Height_vs_Altitude.eps')

#Number Density vs Altitude 100km - 800km
fig3 = plt.figure(figsize=[12.15,6.06])
ax3 = fig3.add_subplot(111)
ax3.set_title('Number Density vs Altitude')
ax3.set_xlabel('Number Density [m^-3]')
ax3.set_ylabel('Altitude [km]')
ax3.set_ylim([0,800])
ax3.semilogx(n['O'], z*1E-3, n['N'], z*1E-3, n['O2'], z*1E-3, n['N2'], z*1E-3,\
 n['NO'], z*1E-3)
ax3.legend(['O','N', 'O2', 'N2', 'NO'])
fig3.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Number_Density_vs_Altitude_100_800.eps')

#Number Density vs Altitude 250km - 350km
fig4 = plt.figure(figsize=[12.15,6.06])
ax4 = fig4.add_subplot(111)
ax4.set_title('Number Density vs Altitude')
ax4.set_xlabel('Number Density [m^-3]')
ax4.set_ylabel('Altitude [km]')
ax4.set_xlim([1E14, 1E15])
ax4.set_ylim([400,500])
ax4.semilogx(n['O'], z*1E-3, n['N2'], z*1E-3)
ax4.legend(['O','N', 'O2', 'N2', 'NO'])
fig4.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Number_Density_vs_Altitude_400_500.eps')

#Mass Density vs Altitude
fig5 = plt.figure(figsize=[12.15,6.06])
ax5 = fig5.add_subplot(111)
ax5.set_title('Mass Density vs Altitude')
ax5.set_xlabel('Mass Density [km/m^3]')
ax5.set_ylabel('Altitude [km]')
ax5.set_ylim([0,800])
ax5.semilogx(rho['O'], z*1E-3, rho['N'], z*1E-3, rho['O2'], z*1E-3, rho['N2'],\
 z*1E-3, rho['NO'], z*1E-3)
ax5.legend(['O','N', 'O2', 'N2', 'NO'])
fig5.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Mass_Density_vs_Altitude.eps')

#Total Weighted Molecular Mass vs Altitude
fig6 = plt.figure(figsize=[12.15,6.06])
ax6 = fig6.add_subplot(111)
ax6.set_title('Total Weighted Molecular Mass vs Altitude')
ax6.set_xlabel('Total Weighted Molecular Mass [amu]')
ax6.set_ylabel('Altitude [km]')
ax6.set_ylim([0,800])
ax6.plot(M_molec, z*1E-3)
fig6.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Total_Weighted_Molecular_Mass_vs_Altitude.eps')

#Column Number Density vs Altitude 100km - 800km
fig7 = plt.figure(figsize=[12.15,6.06])
ax7 = fig7.add_subplot(111)
ax7.set_title('Column Number Density vs Altitude')
ax7.set_xlabel('Column Number Density [m^-2]')
ax7.set_ylabel('Altitude [km]')
ax7.set_ylim([0,800])
ax7.semilogx(col_n['O'], z*1E-3, col_n['N'], z*1E-3, col_n['O2'], z*1E-3,\
 col_n['N2'], z*1E-3, col_n['NO'], z*1E-3)
ax7.legend(['O','N', 'O2', 'N2', 'NO'])
fig7.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Column_Number_Density_vs_Altitude.eps')

#Column Mass Density vs Altitude 100km - 800km
fig8 = plt.figure(figsize=[12.15,6.06])
ax8 = fig8.add_subplot(111)
ax8.set_title('Column Mass Density vs Altitude')
ax8.set_xlabel('Column Mass Density [kg/m^2]')
ax8.set_ylabel('Altitude [km]')
ax8.set_ylim([0,800])
ax8.semilogx(col_rho['O'], z*1E-3, col_rho['N'], z*1E-3, col_rho['O2'], z*1E-3,\
 col_rho['N2'], z*1E-3, col_rho['NO'], z*1E-3)
ax8.legend(['O','N', 'O2', 'N2', 'NO'])
fig8.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Column_Mass_Density_vs_Altitude.eps')

#Tau vs Altitude 100km - 800km
fig9 = plt.figure(figsize=[12.15,6.06])
ax9 = fig9.add_subplot(111)
ax9.set_title('Optical Depth vs Altitude')
ax9.set_xlabel('Optical Depth')
ax9.set_ylabel('Altitude [km]')
ax9.set_xlim([1E-6, 1E4])
ax9.set_ylim([0,600])
ax9.semilogx(sumtau, z*1E-3)
fig9.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Optical_Depth_vs_Altitude.eps')

#Flux vs Altitude 100km - 800km
fig10 = plt.figure(figsize=[12.15,6.06])
ax10 = fig10.add_subplot(111)
ax10.set_title('Flux vs Altitude')
ax10.set_xlabel('Flux [photons/m2*s]')
ax10.set_ylabel('Altitude [km]')
ax10.set_ylim([0,600])
ax10.plot(I, z*1E-3)
fig10.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Flux_vs_Altitude.eps')

#Energy [W/m3] vs Altitude [km] [100km - 800km]
fig11 = plt.figure(figsize=[12.15,6.06])
ax11 = fig11.add_subplot(111)
ax11.set_title('Energy Dissipation vs Altitude')
ax11.set_xlabel('Energy Dissipation [W/m3]')
ax11.set_ylabel('Altitude [km]')
ax11.set_xlim([1E-11, 1E-7])
ax11.set_ylim([0,400])
ax11.semilogx(sumQ, z*1E-3)
fig11.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 2/\
Figures/Energy_Dissipation_vs_Altitude.eps')

print "We took %f seconds!" % (time.time()-start)
