#! /usr/bin/env python3
# Last Modification: 20-06-2017, Eric Andersson <eaherrgarden@gmail.com>
#=======================================================================
# plot.py
#
# Contains all plotting function used for adaptive particles.
#
# Eric Andersson, 20-06-2017
#=======================================================================
def cumulative_density(plotdir = './plots/', xlim = False, ylim = False, xlog = True, ylog = True, filename = False):
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

def cumulative_density_mean(t0, t,
        datadir = './runs/', plotdir = './plots/',
    xlim = (1e-3, 1e3), ylim = (1e-5, 1e0), xlog = True, ylog = True,
    filename = None):
    """ Plots the cumulative particle density distribution averaged
    over time to obtain a stable result.

    Positional Arguments:
        t0
            Initial timestep
        t
            Number of timesteps for average

    Keyword Arguments:
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

    # Set up axis
    if xlog and ylog:
        plt.xscale('log')
        plt.yscale('log')
    elif xlog:
        plt.xscale('log')
    elif ylog:
        plt.yscale('log')
    if type(xlim) is tuple:
        plt.xlim(xlim)
    if type(ylim) is tuple:
        plt.ylim(ylim)

    if datadir == './runs/':
        runs = [line.rstrip('\n') for line in open(datadir + 'runs.txt')]
    else:
        print('Wrong directory.')
        exit()

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
            rhop = np.zeros((t, dim.nxgrid*dim.nygrid*dim.nzgrid))
            if ndim == 1:
                rhop[i-t0][:] = np.sort(f.rhop)[::-1]
            elif ndim == 2:
                rhop[i-t0][:] = np.sort(
                        np.concatenate(f.rhop))[::-1]
            else:
                rhop[i-t0][:] = np.sort(np.concatenate(np.concatenate(
                                        f.rhop)))[::-1]
        rhopcum = np.mean(rhop, axis=0)
        param = pc.read.parameters(datadir = datadir + run + '/data/')
        plt.step(rhopcum/np.mean(rhopcum),
            np.arange(rhopcum.size)/rhopcum.size,
            label = r'$\tau_s = {},\ \epsilon = {}$'.format(
                param.taus, param.eps_dtog))

    plt.minorticks_on()
    plt.xlabel(r"$\rho_p/<\rho_p>$")
    plt.ylabel(r"P$\ (<\rho_p)$")
    plt.legend(loc='best')

    # Save the plot
    if type(filename) is str:
        plt.savefig(plotdir + filename + '.pdf', bbox_inches='tight')

    plt.show()



#======================================================================

def density_histogram_mean(t0, t,
        datadir = './runs/', plotdir = './plots/',
    xlim = (1e-4, 1e3), ylim = (0, 600), xlog = True, ylog = False,
    filename = None):
    """ Plots the cumulative particle density distribution averaged
    over time to obtain a stable result.

    Positional Arguments:
        t0
            Initial timestep
        t
            Number of timesteps for average

    Keyword Arguments:
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
    stdrhop = {}
    avgrhop = {}
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
                rhop[i-t0][:] = np.histogram(f.rhop/np.mean(f.rhop),
                                            bins = bins)[0]
            elif ndim == 2:
                rhop[i-t0][:] = np.histogram(np.concatenate(
                                                f.rhop/np.mean(f.rhop)),
                                            bins = bins)[0]
            else:
                rhop[i-t0][:] = np.histogram(
                                            np.concatenate(np.concatenate(
                                                f.rhop/np.mean(f.rhop))),
                                            bins = bins)[0]
        avgrhop[run] = np.mean(rhop, axis=0)
        stdrhop[run] = np.std(rhop, axis=0)

        sys.stdout = open(os.devnull, 'w')
        param = pc.read.parameters(datadir = datadir + run + '/data/')
        sys.stdout = sys.__stdout__
        plt.step(bins[:-1], np.mean(rhop, axis=0),
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

    # Create plot for standard deviations.
    # Set up axis
    if xlog and ylog:
        plt.xscale('log')
        plt.yscale('log')
    elif xlog:
        plt.xscale('log')
    elif ylog:
        plt.yscale('log')
    if type(xlim) is tuple:
        plt.xlim(xlim)

    for run in runs:
        sys.stdout = open(os.devnull, 'w')
        param = pc.read.parameters(datadir = datadir + run + '/data/')
        sys.stdout = sys.__stdout__
        plt.step(bins[:-1], stdrhop[run]/avgrhop[run],
                label = r'$\tau_s = {},\ \epsilon = {}$'.format(
                param.taus, param.eps_dtog))
    plt.minorticks_on()
    plt.xlabel(r"$\rho_p/<\rho_p>$")
    plt.ylabel(r"Relative standard deviation, $\sigma_{\rho_p}/N(\rho_p)$")
    plt.legend(loc='best')

    # Save the plot
    if type(filename) is str:
            plt.savefig(plotdir + filename + '_std.pdf', bbox_inches='tight')
    plt.show()

#=======================================================================

def rhop_histogram(t0 = 25, t = 70,
    plotdir = './lunarc/nobackup/user/ericand/plots/',
    xlim = (1e-4, 1e3), ylim = (0, 600), xlog = True, ylog = False,
    filename = 'rhop_histogram', add_std = False):
    """ Plots the average number of gridcells with certain density in a
    histogram. This function is mean for plotting data from Nonlinear
    Streaming instability simulations and defaults to data stored on
    ccyang's lunarc directory for these kind of simulations.

    Positional Arguments:

    Keyword Arguments:
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

    # Initialte needed variables
    data = {}
    ndata = -1

    # Read in data
    print('Starting process of reading in data...')
    done = ''
    while done != 'no':
        ndata += 1
        datadir = input('Please give directory for data: ')
        param = pc.read.parameters(datadir=datadir)
        rhop, std = _create_density_histogram(t0, t, datadir, bins)
        data[ndata] = [rhop, std, param.taus, param.eps_dtog]
        done = input('Do you wish to add more data? (yes/no): ')

    for i in range(0, ndata+1):
        if add_std:
            col = (i/(ndata+1), i/(ndata+1), i/(ndata+1))
            plt.fill_between(bins[:-1], data[i][0]-data[i][1],
                    data[i][0] + data[i][1], step='pre', color=col)
        plt.step(bins[:-1]+0.5*(bins[1]-bins[0]), data[i][0], lw = 2,
            label = r'$\tau_s = {},\ \epsilon = {}$'.format(
                data[i][2], data[i][3]), zorder=9)

    plt.legend(loc='best')

    # Save the plot
    if type(filename) is str:
            plt.savefig(plotdir + filename + '.pdf',
                    bbox_inches='tight')
    plt.show()

#=======================================================================
# LOCAL FUNCTIONS
#=======================================================================

def _create_density_histogram(t0, nt, datadir, bins):
    """ Local function for creating a density histogram of pencil code
    data.

    Keyword Arguments:
        t0
            Number of timesteps for warmup
        nt
            Number of timesteps for average
        datadir
            Directory for locating data
        bins
            Shape of the histogram-bins

    Return Values
        mean
            Binned averaged density
        sigma
            Standard deviation of the average
    """
    # Eric Andersson, 13-07-2017
    import numpy as np
    import PencilCode as pc
    import sys, os
    from scipy.integrate import simps

    # Allocate memory
    Nrhop = np.zeros((nt, bins.size - 1))
    t = np.zeros(nt)

    for i in range(t0, t0+nt):
        f = pc.read.var(datadir = datadir, ivar = i)
        dim = pc.read.dimensions(datadir = datadir)
        ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
        if ndim == 1:
            Nrhop[i-t0][:] = np.histogram(f.rhop/np.mean(f.rhop),
                                bins = bins)[0]
        elif ndim == 2:
            Nrhop[i-t0][:] = np.histogram(np.concatenate(
                                            f.rhop/np.mean(f.rhop)),
                                bins = bins)[0]
        else:
            Nrhop[i-t0][:] = np.histogram(np.concatenate(np.concatenate(
                                            f.rhop/np.mean(f.rhop))),
                                bins = bins)[0]
        t[i-t0] = f.t

    mean = simps(y=Nrhop, x=t, axis=0)/(t[-1] - t[0])
    sigma2 = simps(y=Nrhop**2, x=t, axis=0)/(t[-1] - t[0]) - mean**2
    return mean, np.sqrt(sigma2)
