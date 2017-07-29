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

def cumulative_density(animdir = './animations/', xlim = None, ylim = None, xlog = True, ylog = True, filename = 'animation', fps = 24, show_off = False):
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

def density_histogram(plane, animdir = './animations/',
        xlim = (0.001, 100), ylim = (0, 600), xlog = True, ylog = False,
        filename = 'histogram', fps = 24, normed = False,
        show_off = False, add_mean = False):
    """ Plots the cumulative particle density distribution.

    Positional Argument:
        plane
            The plane in which density is measured

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
        add_mean
            If true, adds mean and standard deviation in background
    """
    # Eric Andersson, 22-07-2017
    import PencilCode as pc
    import AdaptiveParticles as ap
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation as animation

    # Read in data
    print('Reading in data...')
    if xlog:
        bins = np.logspace(np.log10(xlim[0]),
                                    np.log10(xlim[1]), 50)
    else:
        bins = np.linspace(xlim[0], xlim[1], 50)
    while True:
        try:
            datadir = input('Please give directory for data: ')
            p = pc.read.parameters(datadir=datadir)
            pd = pc.read.pardim(datadir=datadir)
            g = pc.read.grid(datadir=datadir)
            t, rhop = pc.read.slices('rhop')
            break
        except FileNotFoundError:
            print('The directory does not exist. Try again.')
    res = [g.x.size-6, g.z.size-6]
    npar= int(pd.npar / (res[0]*res[1]))
    data = [p.taus, p.eps_dtog, res, npar]
    print('Done.')

    # Set up figure
    print("Setting up figure...")
    fig  = plt.figure()
    ax = fig.add_subplot(111)
    if add_mean == True:
        mean, std = ap.compute.rhop_histogram(t0 = 300,
                datadir = datadir, normed = normed, plane = plane,
                bins=bins)
        plt.fill_between(bins[:-1] + 0.5*(bins[1]-bins[0]),
                            mean - std, mean + std,
                            step='pre', alpha = 0.3)
        plt.step(bins[:-1]+0.5*(bins[1]-bins[0]), mean, '--',
                lw = 1, zorder=3)

    ax.text(0.75, 0.90, r'$\tau_s = {},\ \epsilon = {}$'.format(
                data[0], data[1]), transform=ax.transAxes)
    ax.text(0.75, 0.85, r'${}\times{}$'.format(
                data[2][0], data[2][1]), transform=ax.transAxes)
    ax.text(0.75, 0.80, r'$np = {}$'.format(
                data[3]), transform=ax.transAxes)
    # Set axis scales
    if xlog and ylog:
        plt.xscale('log')
        plt.yscale('log')
    elif xlog:
        plt.xscale('log')
    elif ylog:
        plt.yscale('log')
    if type(ylim) is tuple:
        plt.ylim(ylim)
    if type(xlim) is tuple:
        plt.xlim(xlim)
    plt.minorticks_on()
    plt.xlabel(r'$\rho_p/<\rho_p>$')
    if normed == True:
        plt.ylabel(r'$N_p/N_{tot}$')
    else:
        plt.ylabel(r'$N_p$')

    # Set up plots
    line, = ax.step([],[], '-k', zorder = 9)
    time_text = ax.text(0.75, 0.95, '', transform=ax.transAxes)
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
        line.set_data(bins[:-1]+0.5*(bins[1]-bins[0]),
                np.histogram(np.concatenate(rhop[i][plane]),
                    normed = normed, bins = bins)[0])
        return line, time_text

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=t.size, interval=200,
                                blit=False, repeat = False,
                                save_count=0)
    anim.save(animdir + filename + '.mp4', fps = fps,
                    extra_args=['-vcodec', 'libx264'])
    print("Done.")
    if show_off:
        plt.clear()
    else:
        plt.show()

#======================================================================
