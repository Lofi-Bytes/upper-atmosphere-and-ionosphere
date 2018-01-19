# -*- coding: utf-8 -*-
"""
Created on Sat Feb 18 23:23:29 2012

@author: jsni
"""


    # Assigning an empty array for NO as a place holder.
    data['NO+'] = np.zeros(37) 
    
    
    
    
 + R39*M39 + R45*M45
 
  + R37*M37 + R46*M46 + R48*M48








    nn['NO+'] = n['NO']


def ion_production_rate(n, I, data, ztop):
    # Creating a new number density of the proper key type so python doesn''t shit itself.
    nn = {}
    nn['O+'] = n['O']
    nn['N+'] = n['N']
    nn['O2+'] = n['O2']
    nn['N2+'] = n['N2']
    # Defining an ion production rate.
    Ps = {}
    for i in ions:
        Ps[i] = np.zeros(1000.0)
        for j in range(0, 37):
            Ps[i][0:ztop] += nn[i][0:ztop]*I[0:ztop,j]*data[i][j]
    return Ps
    
    
    
    
    
    def calculate_ion_number_densities(source, RM, time, ztop):
    # Calculating the ion number densities as they vary in time.
    ni = {}
    for s in ions:
        ni[s] = np.zeros(1000)
        for i in range(4320):
            ni[s] = (ni[s-1] + source[s]*(time[i] - time[i-1]))/(1 + RM[s]*(time[i] - time[i-1])) 
    return ni







    for s in neutrals:
        n[s] = np.zeros(1000.0)
        n[s][0] = n_base[s]
        for i in range(1, 1000):
            n[s][i] = (n[s][i-1]*T[i-1])/(T[i])*np.exp(-(z[i]-z[i-1])/H[s][i])  
            
            
            
            
            
            
            
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
Assignment 3/shelve.out'
my_shelf = shelve.open(filename)
for key in my_shelf:
    globals()[key]=my_shelf[key]
my_shelf.close()

#***********************************************************************
# Building the ionosphere                                              *
#***********************************************************************

#Creating a list of all of the ion species tht live in our virtual
#thermosphere so that we can later assign properties to them.
ions = ['O+', 'N+', 'O2+', 'N2+', 'NO+']

# The rest of the ionosphere is in the loop.

ni = {}
for s in ions:
    ni[s] = np.zeros(1000)



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
Ps = ip.ion_production_rate(n, I, data, ztop)

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
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            
            