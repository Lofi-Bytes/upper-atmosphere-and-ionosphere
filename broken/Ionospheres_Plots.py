# -*- coding: utf-8 -*-
"""
Created on Fri Feb 17 00:55:57 2012

@author: jsni
"""

import matplotlib.pyplot as plt                         # Matplotlib tells python how to make pretty plots.

def plot(target, plotname, xs, ys, xlabel=None, ylabel=None, title=None, 
    legendlabels=None, loc=111):
    '''
    X and Y must be lists of x-vectors and corresponding y-vectors.
    Labels must be a list of the same size of corresponding line labels.
    '''
    # Set up plot.
    # from matplotlib.pyplot import Figure, Axes, figure
    # commented out due to redundancy
    
    # DEBUG!  Use these three lines to confirm whether we are inputting tuples
    # print 'X = ', type(xs)
    # print 'Y = ', type(ys)
    # print 'labels = ', legendlabels
    
    
    # Handle the plot target.
    if type(target)==plt.Figure:
         # Target is a figure object.  Place a new axis on that figure.
         fig = target # Our figure is obtained from the target.
         ax  = fig.add_subplot(loc) # Make a new ax at location=loc.
    elif type(target).__base__ == plt.Axes:
         # Target is an Axes, use that axes and its figure.
         ax  = target
         fig = ax.plt.figure
    else:
         # User did not specify target, make plot from scratch.
         fig = plt.figure(figsize=(8,8))
         ax  = fig.add_subplot(loc)
    
    # If I have just one set of x and y values, plot like normal
    if type(xs) != type( () ):
        # print "PLOTTING NORMALLY!"
        ax.plot(xs, ys)
    # I can plot multiple lines with input format (a, b, c)
    # Multiple values will include a legend
    else:
        for x, y, lab in zip(xs, ys, legendlabels):
            print lab
            ax.plot(x, y, label=lab)
        ax.legend(loc=1)
    # After plotting either singular functions or multiple, apply labels
    # and save the plot
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.savefig('/Users/jsni/Desktop/Ionospheres Assignments/Assignment 3/'+plotname+'.eps')
