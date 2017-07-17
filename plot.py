#! /usr/bin/env python3
# Last Modification: 20-06-2017, Eric Andersson <eaherrgarden@gmail.com>
#=======================================================================
# plot.py
#
# Contains all plotting function used for adaptive particles.
#
# Eric Andersson, 20-06-2017
#=======================================================================
def cumulative_density_snap(plotdir = './plots/', xlim = False, ylim = False, xlog = True, ylog = True, filename = False):
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
    plt.ylabel(r"P$\ (<\rho_p)$")
    plt.legend(loc='best')

    # Save the plot
    if type(filename) is str:
        plt.savefig(plotdir + filename + '.pdf', bbox_inches='tight')

    plt.show()

#======================================================================

def cumulative_density(t0=25, t=70,
        datadir = './runs/', plotdir = './plots/',
    xlim = (1e-3, 1e3), ylim = (1e-5, 1e0), xlog = True, ylog = True,
    filename = None):
    """ Plots the cumulative particle density distribution averaged
    over time to obtain a stable result.

    Keyword Arguments:
        t0
            Initial timestep
        t
            Number of timesteps for average
        datadir
            Directory of data files.
        plotdir
            Directory for saving plots
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

    # Set up plot
    fig = plt.figure()
    plt.minorticks_on()
    plt.xlabel(r"$\rho_p/<\rho_p>$")
    plt.ylabel(r"$P\ (<\rho_p)$")

    # Set up axis
    if xlog and ylog:
        plt.xscale('log')
        plt.yscale('log')
    elif xlog:
        plt.xscale('log')
    elif ylog:
        plt.yscale('log')
    plt.xlim(xlim)
    plt.ylim(ylim)

    # Initiate variables
    data = {}
    ndata = -1

    # Read in data
    print('Starting process of reading in data...')
    done = ''
    while done != 'no':
        ndata += 1
        while True:
            try:
                datadir = input('Please give directory for data: ')
                param = pc.read.parameters(datadir=datadir)
                break
            except FileNotFoundError:
                print('The directory does not exist. Try again.')
        Prhop, std = _cumulative_distribution(t0, t, datadir)
        data[ndata] = [Prhop, std, param.taus, param.eps_dtog]
        done = input('Do you wish to add more data? (yes/no): ')

    # Create plot
    for i in data:
        plt.step(data[i][0]/np.mean(data[i][0]),
                np.arange(data[i][0].size)/data[i][0].size,
                label = r'$\tau_s = {},\ \epsilon = {}$'.format(
                    data[i][2], data[i][3]))
    plt.legend(loc='best')

    # Save and show the plot
    if type(filename) is str:
            plt.savefig(plotdir + filename + '.pdf',
                    bbox_inches='tight')
    plt.show()

#======================================================================

def rhop_histogram(t0 = 25, t = 70,
    plotdir = './work/plots/', label = None,
    xlim = (1e-4, 1e3), ylim = (0, 600), xlog = True, ylog = False,
    filename = 'rhop_histogram', add_std = False, normed = False):
    """ Plots an average of number of gridcells with certain density in
    a histogram. The function plots data from Nonlinear Streaming
    instability simulations.

    Keyword Arguments:
        t0
            Initial timestep
        t
            Number of timesteps for average
        datadir
            Directory of data files.
        plotdir
            Directory for saving plots
        label
            Specify whether to use resolution for label, or parameters.
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
        add_std
            Adds standard devtiation to the plotted lines.
        normed
            Set True to create normalized histogram.
    """
    # Eric Andersson, 13-07-2017
    import PencilCode as pc
    import numpy as np
    import matplotlib.pyplot as plt

    # Set up plot
    fig = plt.figure()
    plt.minorticks_on()
    plt.xlabel(r"$\rho_p/<\rho_p>$")
    plt.ylabel(r"$N({\rho_p})$")

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

    # Initiate variables
    data = {}
    ndata = -1

    # Read in data
    print('Starting process of reading in data...')
    done = ''
    while done != 'no':
        ndata += 1
        while True:
            try:
                datadir = input('Please give directory for data: ')
                param = pc.read.parameters(datadir=datadir)
                break
            except FileNotFoundError:
                print('The directory does not exist. Try again.')
        rhop, std = _create_density_histogram(t0, t,
                                        datadir, bins, normed)
        data[ndata] = [rhop, std, param.taus, param.eps_dtog]
        done = input('Do you wish to add more data? (yes/no): ')

    # Write data to plot.
    for i in data:
        if label == 'resolution':
            g = pc.read.grid(datadir=datadir)
            label = r'${}\times{}$'.format(g.x.size-6, g.z.size-6)
        elif label == 'parameters':
            label = r'$\tau_s = {},\ \epsilon = {}$'.format(
                data[i][2], data[i][3])
        else:
            g = pc.read.grid(datadir=datadir)
            label = r'$\tau_s = {},\ \epsilon = {},\ {}\times{}$'.format(
                data[i][2], data[i][3], g.x.size-6, g.z.size-6)
        if add_std:
            col = (i/(ndata+1), i/(ndata+1), i/(ndata+1))
            plt.fill_between(bins[:-1], data[i][0]-data[i][1],
                    data[i][0] + data[i][1], step='pre', alpha = 0.3)
        plt.step(bins[:-1]+0.5*(bins[1]-bins[0]), data[i][0], lw = 1,
                label = label, zorder=9)

        plt.legend(loc='best', prop = {'size':10}, frameon=False)

    # Save and show the plot
    if type(filename) is str:
            plt.savefig(plotdir + filename + '.pdf',
                    bbox_inches='tight')
    plt.show()

#=======================================================================
# LOCAL FUNCTIONS
#=======================================================================

def _create_density_histogram(t0, nt, datadir, bins, normed):
    """ Local function for creating a density histogram of pencil code
    data.

    Positional Arguments:
        t0
            Number of timesteps for warmup
        nt
            Number of timesteps for average
        datadir
            Directory for locating data
        bins
            Shape of the histogram-bins
        normed
            If True, then plot histogram will be normalized
    Return Values
        mean
            Binned averaged density
        sigma
            Standard deviation of the average
    """
    # Eric Andersson, 13-07-2017
    import numpy as np
    import PencilCode as pc
    from scipy.integrate import simps

    # Allocate memory
    Nrhop = np.zeros((nt, bins.size - 1))
    t = np.zeros(nt)

    # Read data
    dim = pc.read.dimensions(datadir = datadir)
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    for i in range(0, nt):
        f = pc.read.var(datadir = datadir, ivar = t0+i)
        if ndim == 1:
            Nrhop[i][:] = np.histogram(f.rhop/np.mean(f.rhop),
                                bins = bins, normed=normed)[0]
        elif ndim == 2:
            Nrhop[i][:] = np.histogram(np.concatenate(
                                            f.rhop/np.mean(f.rhop)),
                                bins = bins, normed=normed)[0]
        else:
            Nrhop[i][:] = np.histogram(np.concatenate(np.concatenate(
                                            f.rhop/np.mean(f.rhop))),
                                bins = bins, normed=normed)[0]
        t[i] = f.t

    mean = simps(y=Nrhop, x=t, axis=0)/(t[-1] - t[0])
    sigma2 = simps(y=Nrhop**2, x=t, axis=0)/(t[-1] - t[0]) - mean**2
    return mean, np.sqrt(sigma2)

#=======================================================================

def _cumulative_distribution(t0, nt, datadir):
    """ Local function that sets up cumulative denisty distribution.

    Positional Arguments:
        t0
            Number of timesteps for warmup
        nt
            Number of timesteps for average
        datadir
            Directory for locating data
    """
    # Eric Andersson, 16-07-2017
    import numpy as np
    import PencilCode as pc
    from scipy.integrate import simps

    # Allocate memory
    dim = pc.read.dimensions(datadir = datadir)
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    t = np.zeros(nt)
    Prhop = np.zeros((nt, dim.nxgrid*dim.nygrid*dim.nzgrid))

    # Read data
    for i in range(0, nt):
        f = pc.read.var(datadir = datadir, ivar = t0+i)
        if ndim == 1:
            Prhop[i][:] = np.sort(f.rhop)[::-1]
        elif ndim == 2:
            Prhop[i][:] = np.sort(
                                np.concatenate(f.rhop))[::-1]
        else:
            Prhop[i][:] = np.sort(np.concatenate(np.concatenate(
                                                        f.rhop)))[::-1]
        t[i] = f.t

    # Calculate mean
    mean = simps(y=Prhop, x=t, axis=0)/(t[-1] - t[0])
    sigma2 = simps(y=Prhop**2, x=t, axis=0)/(t[-1] - t[0]) - mean**2

    return mean, np.sqrt(sigma2)

