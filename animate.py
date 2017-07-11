#! /usr/bin/env python3
# Last Modification: 07-07-2017, Eric Andersson <eaherrgarden@gmail.com>
#=======================================================================
# animate.py
#
# Contains all animation functions used for adaptive particles.
#
# Eric Andersson, 30-06-2017
#=======================================================================

def grid_density_1D(filename, dim, animdir = './animations/',
        xlim = None, ylim = None, xlog = None, ylog = None, fps = 24,
        show_off = False):
    """ Creates an animation of how the grid-density changes in time.

    Positional Arguments:

        filename
            Name of the output file
        dim
            Dimension along which rho will be plotted

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
        fps
            Frames per second of video.
        show_off
            If true, turns off playing the animation after finilizing

    """
    # Eric Andersson, 30-06-2017
    import PencilCode as pc
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    # Parameters used in run
    print("Setting up parameters...")
    param = pc.read.parameters()
    taus = param.taus               # Dimensionless stopping time
    eps = param.eps_dtog            # Local mass density ratio
    print("Done.")

    # Read in data
    print("Reading data...")
    f = pc.read.var()
    t, rho = pc.read.slices("rhop")
    print("Done.")

    # Check dimension
    if dim == 'x':
        plane = 'xz'
        dimindex = 0
        x = f.x
    elif dim == 'y':
        plane = 'xy'
        dimindex = 1
        x = f.y
    else:
        plane = 'xz'
        dimindex = 1
        x = f.x

    # Set up figure
    print("Setting up figure...")
    fig  = plt.figure()
    ax = fig.add_subplot(111)
    line, = ax.step([], [],
            label = r'$\tau_s = {},\ \epsilon = {}$'.format(
                taus, eps))
    plt.minorticks_on()
    plt.xlabel(dim)
    plt.ylabel(r'$\rho$')
    plt.legend(loc = 'upper right')

    # Set up time text
    time_text = ax.text(0.02, 0.95, '', transform=ax.transAxes)

    # Set axis scales
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

    print("Done.")

    # Function that sets up background of each frame
    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text


    # Function that is called everytime a new frame is produced
    def animate(i):
        print("\rAnimating ({:6.1%})......".format(i/t.size),
                end='', flush=True)
        time_text.set_text(r'Time = {0:.2f} $\Omega t$'.format(t[i]))
        line.set_data(x, rho[i][plane][dimindex])
        return line, time_text

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=t.size, interval=200, blit=False)
    anim.save(filename + '.mp4', fps = fps,
            extra_args=['-vcodec', 'libx264'])
    print("Done.")
    if show_off:
        plt.clear()
    else:
        plt.show()
#=======================================================================

def cumulative_rhop_2D(animdir = './animations/', xlim = None, ylim = None, xlog = True, ylog = True, filename = 'animation', fps = 24, show_off = False):
    """ Plots the cumulative particle density distribution.

    Keyword Argument:
        animdir
            Directory for saving animation.
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
            Name of output file.
        fps
            Frames per second in animation
        show_off
            If true, turns off playing the animation after finilizing

    """
    # Eric Andersson, 06-07-2017
    import PencilCode as pc
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    ### WARNING ###
    # ONLY WORKS IN 2D SIMULATIONS!
    dim = pc.read.dimensions()
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim == 3:
        raise NotImplementedError(
                "More than 2 dimensions is not implemented yet")

    # Import and set-up data from the pencil-code
    print("Reading and setting up data")
    t, rhop = pc.read.slices('rhop')
    plane = rhop.dtype.names
    for p in plane:
        if (rhop[0][p].shape[0] != 1) and (rhop[0][p].shape[1] != 1):
            plane = p
            break
    print("Done.")
    print("Working along the " + plane + "-plane")


    # Parameters used in run
    print("Setting up parameters...")
    param = pc.read.parameters()
    taus = param.taus               # Dimensionless stopping time
    eps = param.eps_dtog            # Local mass density ratio
    print("Done.")

    # Set up figure
    print("Setting up figure...")
    fig  = plt.figure()
    ax = fig.add_subplot(111)
    line, = ax.step([], [],
            label = r'$\tau_s = {},\ \epsilon = {}$'.format(
                taus, eps))
    plt.minorticks_on()
    plt.xlabel(r'$\rho_p/<\rho_p>$')
    plt.ylabel(r'$P(<\rho_p)$')
    plt.legend(loc = 'lower left')

    # Set up time text
    time_text = ax.text(0.7, 0.95, '', transform=ax.transAxes)

    # Set axis scales
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

    print("Done.")

    # Function that sets up background of each frame
    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text


    # Function that is called everytime a new frame is produced
    def animate(i):
        print("\rAnimating ({:6.1%})......".format(i/t.size),
                end='', flush=True)
        time_text.set_text(r'Time = {0:.2f} $\Omega t$'.format(t[i]))
        rhopcum = np.sort(np.concatenate(rhop[i][plane]))
        avgrhop = np.mean(rhopcum)
        line.set_data(rhopcum[::-1]/avgrhop,
            np.arange(rhopcum.size)/rhopcum.size)
        return line, time_text

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=t.size, interval=200,
                                blit=False, repeat = False)
    anim.save(animdir + filename + '.mp4', fps = fps,
            extra_args=['-vcodec', 'libx264'])
    print("Done.")
    if show_off:
        plt.clear()
    else:
        plt.show()

#======================================================================
def rhop_histogram_2D(animdir = './animations/', xlim = (0.01, 10),
        ylim = (0, 600), xlog = False, ylog = False,
        filename = 'animation', fps = 24, show_off = False):
    """ Plots the cumulative particle density distribution.

    Keyword Argument:
        animdir
            Directory for saving animation.
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
            Name of output file.
        fps
            Frames per second in animation
        show_off
            If true, turns off playing the animation after finilizing
    """
    # Eric Andersson, 06-07-2017
    import PencilCode as pc
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    # ONLY WORKS IN 2D SIMULATIONS!
    dim = pc.read.dimensions()
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim == 3:
        raise NotImplementedError(
                "Detected 3D. For 3D use animate.rhop_histogram_3D")

    # Import and set-up data from the pencil-code
    print("Reading and setting up data")
    t, rhop = pc.read.slices('rhop')
    plane = rhop.dtype.names
    for p in plane:
        if (rhop[0][p].shape[0] != 1) and (rhop[0][p].shape[1] != 1):
            plane = p
            break
    print("Done.")
    print("Working along the " + plane + "-plane")


    # Parameters used in run
    print("Setting up parameters...")
    param = pc.read.parameters()
    taus = param.taus               # Dimensionless stopping time
    eps = param.eps_dtog            # Local mass density ratio
    print("Done.")

    # Set up figure
    print("Setting up figure...")
    fig  = plt.figure()
    ax = fig.add_subplot(111)
    print("Done.")

    # Set up axis
    if ylog:
        ax.set_yscale('log')
    if type(xlim) is tuple:
        if xlog:
            bins = np.logspace(np.log10(xlim[0]),
                                    np.log10(xlim[1]), 50)
        else:
            bins = np.linspace(xlim[0], xlim[1], 50)
    else:
        bins = 'auto'

    # Function that sets up background of each frame
    def init():
        ax.text(0.75, 0.85, r'$\tau_s = {},\ \epsilon = {}$'.format(
                taus, eps), transform=ax.transAxes)
        plt.minorticks_on()
        plt.xlabel(r'$\rho_p/<\rho_p>$')
        plt.ylabel(r'$N_p$')
        if type(ylim) is tuple:
            plt.ylim(ylim)
        # Set axis scales
        if xlog and ylog:
            plt.xscale('log')
            plt.yscale('log')
        elif xlog:
            plt.xscale('log')
        elif ylog:
            plt.yscale('log')
        # Set up time text
        time_text = ax.text(0.75, 0.95,
            r'Time = {0:.2f} $\Omega t$'.format(t[0]),
            transform=ax.transAxes)
        ax.hist(np.concatenate(
                    rhop[0]['xz'])/np.mean(rhop[0][plane]),
                bins = bins)

    # Function that is called everytime a new frame is produced
    def animate(i):
        print("\rAnimating ({:6.1%})......".format(i/t.size),
                end='', flush=True)
        # Reset figure.
        ax.clear()
        plt.minorticks_on()
        plt.xlabel(r'$\rho_p/<\rho_p>$')
        plt.ylabel(r'$N_{\rm cells}$')
        if type(ylim) is tuple:
            plt.ylim(ylim)
        # Set axis scales
        if xlog and ylog:
            plt.xscale('log')
            plt.yscale('log')
        elif xlog:
            plt.xscale('log')
        elif ylog:
            plt.yscale('log')
        # Plot data
        ax.text(0.75, 0.85, r'$\tau_s = {},\ \epsilon = {}$'.format(
                taus, eps), transform=ax.transAxes)
        ax.text(0.75, 0.95,
            r'Time = {0:.2f} $\Omega t$'.format(t[i]),
            transform=ax.transAxes)
        ax.hist(np.concatenate(
                        rhop[i]['xz'])/np.mean(rhop[i][plane]),
                bins = bins)

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=t.size, interval=200,
                                blit=False, repeat = False)
    anim.save(animdir + filename + '.mp4', fps = fps,
            extra_args=['-vcodec', 'libx264'])
    print("Done.")
    if show_off:
        plt.clear()
    else:
        plt.show()

#======================================================================

def rhop_histogram_3D(tmax, animdir = './animations/',
        xlim = (0.01, 10), ylim = (0, 600), xlog = False, ylog = False,
        filename = 'animation', fps = 24, show_off = False, dsnap = 1):
    """ Plots the cumulative particle density distribution.

    Keyword Argument:
        animdir
            Directory for saving animation.
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
            Name of output file.
        fps
            Frames per second in animation
        show_off
            If true, turns off playing the animation after finilizing
        dsnap
            Time between each time stamp.
    """
    # Eric Andersson, 06-07-2017
    import PencilCode as pc
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation
    import sys
    import os

    print('Iniziating ...')
    # Check dimensions
    dim = pc.read.dimensions()
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim != 3:
        print('Could not detect 3 dimensions. Will run animate.rhop_histogram_2D instead for better cadence.')
        rhop_histogram_2D(animdir = animdir, xlim = xlim,
        ylim = ylim, xlog = xlog, ylog = ylog,
        filename = filename, fps = fps, show_off = show_off)
        exit()

    # Parameters used in run
    print("Setting up parameters...")
    param = pc.read.parameters()
    taus = param.taus               # Dimensionless stopping time
    eps = param.eps_dtog            # Local mass density ratio
    print("Done.")

    # Set up figure
    print("Setting up figure...")
    fig  = plt.figure()
    ax = fig.add_subplot(111)
    print("Done.")

    # Set up axis
    if ylog:
        ax.set_yscale('log')
    if type(xlim) is tuple:
        if xlog:
            bins = np.logspace(np.log10(xlim[0]),
                                    np.log10(xlim[1]), 50)
        else:
            bins = np.linspace(xlim[0], xlim[1], 50)
    else:
        bins = 'auto'

    # Function that sets up background of each frame
    def init():
        ax.text(0.75, 0.85, r'$\tau_s = {},\ \epsilon = {}$'.format(
                taus, eps), transform=ax.transAxes)
        plt.minorticks_on()
        plt.xlabel(r'$\rho_p/<\rho_p>$')
        plt.ylabel(r'$N_p$')
        if type(ylim) is tuple:
            plt.ylim(ylim)
        # Set axis scales
        if xlog and ylog:
            plt.xscale('log')
            plt.yscale('log')
        elif xlog:
            plt.xscale('log')
        elif ylog:
            plt.yscale('log')
        sys.stdout = open(os.devnull, 'w')
        f = pc.read.var(ivar = 0)
        sys.stdout = sys.__stdout__
        ax.hist(np.concatenate(np.concatenate(
                    f.rhop))/np.mean(f.rhop),
                bins = bins)
        # Set up time text
        time_text = ax.text(0.75, 0.95,
            r'Time = {0:.2f} $\Omega t$'.format(f.t),
            transform=ax.transAxes)

    # Function that is called everytime a new frame is produced
    def animate(i):
        print("\rAnimating ({:6.1%})......".format(i/(tmax/dsnap)),
                end='', flush=True)
        #Read data
        sys.stdout = open(os.devnull, 'w')
        f = pc.read.var(ivar = i)
        sys.stdout = sys.__stdout__
        # Reset figure.
        ax.clear()
        plt.minorticks_on()
        plt.xlabel(r'$\rho_p/<\rho_p>$')
        plt.ylabel(r'$N_{\rm cells}$')
        if type(ylim) is tuple:
            plt.ylim(ylim)
        # Set axis scales
        if xlog and ylog:
            plt.xscale('log')
            plt.yscale('log')
        elif xlog:
            plt.xscale('log')
        elif ylog:
            plt.yscale('log')
        # Plot data
        ax.text(0.75, 0.85, r'$\tau_s = {},\ \epsilon = {}$'.format(
                taus, eps), transform=ax.transAxes)
        ax.text(0.75, 0.95,
            r'Time = {0:.2f} $\Omega t$'.format(f.t),
            transform=ax.transAxes)
        ax.hist(np.concatenate(
            np.concatenate(f.rhop))/np.mean(f.rhop),
                bins = bins)

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=int(tmax/dsnap), interval=200,
                                blit=False, repeat = False)
    anim.save(animdir + filename + '.mp4', fps = fps,
            writer='imagemagick', extra_args=['-vcodec', 'libx264'])
    print("Done.")
    if show_off:
        plt.clear()
    else:
        plt.show()

#======================================================================
