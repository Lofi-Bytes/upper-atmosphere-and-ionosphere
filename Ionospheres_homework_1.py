# -*- coding: utf-8 -*-

"""
Created on Wed Jan 11 19:41:05 2012

@author: Jonathan Nickerson
"""

import numpy as np
import matplotlib.pyplot as plt

#Definitions

theta = np.arange(1000.0)*np.pi/999.0      #Theta -> 0 - pi

T = 800.0*np.tanh(theta*3.0)+200.0         #Temperature -> 200-1000 [K]

z = np.arange(1000.0)/999.0*700.0+100.0    #Altitude -> 100-800 km

R = 6378.1                                 #Radius of Earth [km]

g0 = 9.81/1000                             #Gravity at Earths Surface

g = g0*(R/(R+z))*(R/(R+z))                 #Gravity as it Varies With Altitude 
                                           #[km/s^2]
k = 1.3806503E-23/1E6                      #Boltzmann constant 
                                           #[((km^2)*kg)/((s^2)*K)]
A = 6.022E23                               #Avagadro's Number

#Mass of Each Species                       Convert Mass From AMU to kg
m_O = 16.0/(1000.0*A)
m_N = 14.0/(1000.0*A)
m_O2 = 32.0/(1000.0*A)
m_N2 = 28.0/(1000.0*A)
m_NO = 30.0/(1000.0*A)

#Scale Height of Each Species at Initial Altitude
H_O = (k*T)/(m_O*g)           
H_N = (k*T)/(m_N*g)
H_O2 = (k*T)/(m_O2*g)
H_N2 = (k*T)/(m_N2*g)
H_NO = (k*T)/(m_NO*g)

#Densities of Each Species given in [m^3] at 100km
n0_O  = 4.9E17*1E9                       #Multiply to convert from m^3 to km^3
n0_N  = 6.4E11*1E9
n0_O2 = 2.4E18*1E9
n0_N2 = 1.0E19*1E9
n0_NO = 1.0E14*1E9

#Average Mass Density
rho0_ave = n0_O*m_O + n0_N*m_N + n0_O2*m_O2 + n0_N2*m_N2 + n0_NO*m_NO

#Setting Scale Height as an Average (Mass Density Weighted) Up to 120km
H_ave = np.arange(29)
for i in range (1, 29):
    H_ave[i] = (n0_O*m_O*H_O[i] + n0_N*m_N*H_N[i] + n0_O2*m_O2*H_O2[i] + n0_N2*m_N2*H_N2[i] + n0_NO*m_NO*H_NO[i])/rho0_ave    
    H_O[i] = H_ave[i]
    H_N[i] = H_ave[i]
    H_O2[i] = H_ave[i]
    H_N2[i] = H_ave[i]
    H_NO[i] = H_ave[i]
   
#Create an Empty Array
n_O = np.arange(1000.0)
n_N = np.arange(1000.0)
n_O2 = np.arange(1000.0)
n_N2 = np.arange(1000.0)
n_NO = np.arange(1000.0)

#Setting the First Value in the Array to the Initial Density
n_O[0] = n0_O
n_N[0] = n0_N
n_O2[0] = n0_O2
n_N2[0] = n0_N2
n_NO[0] = n0_NO

#Iteratively Caluculating the Value of n_O and Assigning it to Consecutive 
#Positions in the Array. 
for i in range(1, 1000):    
    n_O[i] = (n_O[i-1]*T[i-1])/(T[i])*np.exp(-(z[i]-z[i-1])/H_O[i])
    n_N[i] = (n_N[i-1]*T[i-1])/(T[i])*np.exp(-(z[i]-z[i-1])/H_N[i])    
    n_O2[i] = (n_O2[i-1]*T[i-1])/(T[i])*np.exp(-(z[i]-z[i-1])/H_O2[i])
    n_N2[i] = (n_N2[i-1]*T[i-1])/(T[i])*np.exp(-(z[i]-z[i-1])/H_N2[i])
    n_NO[i] = (n_NO[i-1]*T[i-1])/(T[i])*np.exp(-(z[i]-z[i-1])/H_NO[i])

#Tell us at What Values the Species are Equal in Number Density


#Defining Mass Densities (as they Vary with Altitude)
rho_O = n_O*m_O
rho_N = n_N*m_N
rho_O2 = n_O2*m_O2
rho_N2 = n_N2*m_N2
rho_NO = n_NO*m_NO

#Total number density as it changes with altitude 
N_tot = n_O + n_N + n_O2 + n_N2 + n_NO

#Total weighted mass as it changes with altitude
M = m_O*(n_O/N_tot) + m_N*(n_N/N_tot) + m_O2*(n_O2/N_tot) + m_N2*(n_N2/N_tot) + m_NO*(n_NO/N_tot)

#Total weighted molecular mass as it changes with altitude
M_molec = M*(1000.0*A)

#Plotting Number Density vs Altitude 100km - 800km
fig1 = plt.figure(figsize=[13.5,7.34])
ax1 = fig1.add_subplot(111)
ax1.set_title('Number Density vs Altitude')
ax1.set_xlabel('log(Number Density)')
ax1.set_ylabel('Altitude')
#ax1.set_xlim([0,1E20])
ax1.set_ylim([0,800])
ax1.semilogx(n_O, z, n_N, z, n_O2, z, n_N2, z, n_NO, z)
ax1.legend(['O','N', 'O2', 'N2', 'NO'])
fig1.savefig('./Figures/Density_vs_Altitude_100_800.png')
'''
#Plotting Mass Density vs Altitude
fig2 = plt.figure(figsize=[13.5,7.34])
ax2 = fig2.add_subplot(111)
ax2.set_title('Mass Density vs Altitude')
ax2.set_xlabel('log(Mass Density)')
ax2.set_ylabel('Altitude')
#ax2.set_xlim([1E-21,1E-6])
ax2.set_ylim([0,800])
ax2.semilogx(rho_O, z, rho_N, z, rho_O2, z, rho_N2, z, rho_NO, z)
ax2.legend(['O','N', 'O2', 'N2', 'NO'])
fig2.savefig('/home/jonathan/Desktop/Ionospheres Homeworks/Homework 1/Figures/Mass_Density_vs_Altitude.png')

#Plotting Temperature vs Altitude
fig3 = plt.figure(figsize=[13.5,7.34])
ax3 = fig3.add_subplot(111)
ax3.set_title('Temperature vs Altitude')
ax3.set_xlabel('Temperature')
ax3.set_ylabel('Altitude')
ax3.set_xlim([200,1100])
ax3.plot(T, z)
fig3.savefig('/home/jonathan/Desktop/Ionospheres Homeworks/Homework 1/Figures/Temperature_vs_Altitude.png')

#Plotting Total Weighted Molecular Mass vs Altitude
fig4 = plt.figure(figsize=[13.5,7.34])
ax4 = fig4.add_subplot(111)
ax4.set_title('Total Weighted Molecular Mass vs Altitude')
ax4.set_xlabel('Total Weighted Molecular Mass')
ax4.set_ylabel('Altitude')
#ax4.set_xlim([1E-21,1E-6])
ax4.set_ylim([0,800])
ax4.plot(M_molec, z)
#ax4.legend([])
fig4.savefig('/home/jonathan/Desktop/Ionospheres Homeworks/Homework 1/Figures/Total_Weighted_Molecular_Mass_vs_Altitude.png')

#Plotting Scale Height vs Altitude
fig5 = plt.figure(figsize=[13.5,7.34])
ax5 = fig5.add_subplot(111)
ax5.set_title('Scale Height vs Altitude')
ax5.set_xlabel('Scale Height')
ax5.set_ylabel('Altitude')
ax5.plot(H_O, z, H_N, z, H_O2, z, H_N2, z, H_NO, z)
ax5.legend(['O','N', 'O2', 'N2', 'NO'])
fig5.savefig('/home/jonathan/Desktop/Ionospheres Homeworks/Homework 1/Figures/Scale_Height_vs_Altitude.png')

#Plotting Number Density vs Altitude 100km - 200km
#const = np.arange(1000)*0.0 + 120.0
fig6 = plt.figure(figsize=[13.5,7.34])
ax6 = fig6.add_subplot(111)
ax6.set_title('Number Density vs Altitude')
ax6.set_xlabel('log(Number Density)')
ax6.set_ylabel('Altitude')
ax6.set_xlim([1E19, 1E28])
ax6.set_ylim([100,140])
ax6.semilogx(n_O, z, n_N, z, n_O2, z, n_N2, z, n_NO, z)
#ax6.semilogx(const)
ax6.legend(['O','N', 'O2', 'N2', 'NO'])
fig6.savefig('/home/jonathan/Desktop/Ionospheres Homeworks/Homework 1/Figures/Density_vs_Altitude_100_200.png')
'''