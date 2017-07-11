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

#======================================================================

def time_averaged_density_histogram(t0, t,
        datadir = './runs/', plotdir = './plots/',
    xlim = (1e-4, 1e3), ylim = (0, 600), xlog = True, ylog = False,
    filename = None, nbins = 50):
    """ Plots the cumulative particle density distribution averaged
    over time to obtain a stable result.

    Positional Arguments:
        t0
            Initial timestep
        t
            Number of timesteps for average

    Keyword Arguments:
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
    # Eric Andersson, 11-07-2017
    import PencilCode as pc
    import numpy as np
    import matplotlib.pyplot as plt
    import sys, os
    # Set up axis
    if xlog and ylog:
        plt.xscale('log')
        plt.yscale('log')
    elif xlog:
        plt.xscale('log')
    elif ylog:
        plt.yscale('log')
    if type(xlim) is tuple:
        if xlog:
            bins = np.logspace(np.log10(xlim[0]),
                                    np.log10(xlim[1]), 50)
        else:
            bins = np.linspace(xlim[0], xlim[1], 50)
    else:
        bins = 'auto'

    if datadir == './runs/':
        runs = [line.rstrip('\n') for line in open(datadir + 'runs.txt')]
    else:
        print('Wrong directory.')
        exit()
    # Allocate array for storing data
    rhop = np.zeros((t, 49))
    # Read data and create plot
    for run in runs:
        print('Plotting ' + run)
        for i in range(t0, t0+t):
            sys.stdout = open(os.devnull, 'w')
            dim = pc.read.dimensions(datadir = datadir + run + '/data/')
            ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
            f = pc.read.var(datadir = datadir + run + '/data/',
                                ivar = i)
            sys.stdout = sys.__stdout__
            if ndim == 1:
                rhop[t0-i][:] = np.histogram(f.rhop/np.mean(f.rhop),
                                            bins = bins)[0]
            elif ndim == 2:
                rhop[t0-i][:] = np.histogram(np.concatenate(
                                                f.rhop/np.mean(f.rhop)),
                                            bins = bins)[0]
            else:
                rhop[t0-i][:] = np.histogram(
                                            np.concatenate(np.concatenate(
                                                f.rhop/np.mean(f.rhop))),
                                            bins = bins)[0]
        avgrhop = np.mean(rhop, axis=0)
        print(rhop.shape)
        print(avgrhop)
        param = pc.read.parameters(datadir = datadir + run + '/data/')
        plt.step(bins[:-1], avgrhop,
            label = r'$\tau_s = {},\ \epsilon = {}$'.format(
                param.taus, param.eps_dtog))

    plt.minorticks_on()
    plt.xlabel(r"$\rho_p/<\rho_p>$")
    plt.ylabel(r"N $(\rho_p)$")
    plt.legend(loc='best')

    # Save the plot
    if type(filename) is str:
        plt.savefig(plotdir + filename + '.pdf', bbox_inches='tight')

    plt.show()


