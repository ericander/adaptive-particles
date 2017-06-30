#! /usr/bin/env python3
# Last Modification: 30-06-2017, Eric Andersson <eaherrgarden@gmail.com>
#=======================================================================
# animate.py
#
# Contains all animation functions used for adaptive particles.
#
# Eric Andersson, 30-06-2017
#=======================================================================

def grid_density(filename, animdir = './animations/',
        xlim = False, ylim = False, xlog = False, ylog = False):
    """ Creates an animation of how the grid-density changes in time.

    Positional Arguments:
        filename
            Name of the output file

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
    # Eric Andersson, 30-06-2017
    import PencilCode as pc
    import numpy as np
    import matplotlib.pyplot as plt
    import matplotlib.animation

    # Set up data
    """ To do!
    taus =          # Dimensionless stopping time
    eps =           # Local mass density ratio
    x =             # Position of grid cells
    rho =           # Density in grid cells
    dt =            # Time step of simulation
    """

    # Set up figure
    fig  = plt.figure()
    line, = plt.step([], [],
            label = r'$\tau_s = {taus},\ \epsilon = {eps}$'.format(
                taus, eps))
    plt.minorticks_on()
    plt.xlabel('x')
    plt.ylabel(r'$\rho$')
    leg = plt.legend(loc = 'best')

    # Set up time text
    plt.draw()
    p = leg.get_window_extent()
    time_text = plt.annotate('',
                (p.p0[0], p.p0[1] - 20), (p.p0[0], p.p0[1]-20),
                xycoords='figure pixels', zorder = 9)

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

    # Function that sets up background of each frame
    def init():
        line.set_data([], [])
        time_text.set_text('')
        return line, time_text

    # Function that is called everytime a new frame is produced
    def animate(i):
        print("\rAnimating ({:6.1%})......".format(i/nsteps),
                end='', flush=True)
        time_text.set_text(r'Time = {0:.2f} $\Omega t$'.format(i*dt))
        line.set_data(x, rho[i])
        return line, time_text

    anim = animation.FuncAnimation(fig, animate, init_func=init,
                                frames=nsteps, interval=1, blit=True)
    anim.save(filename + '.mp4', fps = 24,
            extra_args=['-vcodec', 'libx264'])
    print("Done.")
#=======================================================================
