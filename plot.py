#! /usr/bin/env python3
# Last Modification: 20-06-2017, Eric Andersson <eaherrgarden@gmail.com>
#=======================================================================
# plot.py
#
# Contains all plotting function used for adaptive particles.
#
# Eric Andersson, 20-06-2017
#=======================================================================
def rhopcumdist(plotdir = './plots/', xlim = False, ylim = False, xlog = True, ylog = True, filename = False):
    """ Plots the cumulative particle density distribution.

    Keyword Argument:
        xlim
            If limits is of type tuple it will set the x-axis limit
            to the given values.
        ylim
            If limits is of type tuple it will set the y-axis limit
            to the given values.
        xlog
            True/False whether to use log scale on x-axis or not.
        ylog
            True/False whether to use log scale on y-axis or not.
        filename
            If filename is given as string it will save plot as
            'filename.pdf'.

    """
    # Eric Andersson, 20-06-2017
    import PencilCode as pc
    import numpy as np
    import matplotlib.pyplot as plt

    # Import and set-up data from the pencil-code
    f = pc.read.var()
    rhopcum = np.sort(np.concatenate(f.rhop))
    avgrhop = np.mean(rhopcum)

    # Set axis scales.
    if xlog and ylog:
        plt.xscale('log')
        plt.yscale('log')
    elif xlog:
        plt.xscale('log')
    elif ylog:
        plt.yscale('log')

    # Set axis limits
    if type(xlim) is tuple:
        plt.xlim(xlim)
    if type(ylim) is tuple:
        plt.ylim(ylim)

    # Plot the data
    plt.step(rhopcum[::-1]/avgrhop,
            np.arange(rhopcum.size)/rhopcum.size,
            label = r"$\tau_s = 1.0$, $\epsilon = 0.2$")
    plt.minorticks_on()
    plt.xlabel(r"$\rho_p/<\rho_p>$")
    plt.ylabel(r"P $(<\rho_p)$")
    plt.legend(loc='best')

    # Save the plot
    if type(filename) is str:
        plt.savefig(plotdir + filename + '.pdf', bbox_inches='tight')

    plt.show()
