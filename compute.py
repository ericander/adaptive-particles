#! /usr/bin/env python3
# Last Modification: 21-07-2017, Eric Andersson <eaherrgarden@gmail.com>
#=======================================================================
# compute.py
#
# Fazilitates functions for computing values, parameters etc for
# Adaptive particles.
#
# Eric Andersson, 21-07-2017
#=======================================================================

def rhop_histogram(t0, plane, datadir, bins, normed):
    """ Local function for creating a density histogram of pencil code
    data. Works only for 2D simulations.

    Positional Arguments:
        t0
            Number of timesteps for warmup
        plane
            Choose 'xy', 'xz' or 'yz' plane to work along.
        datadir
            Directory for locating data
        bins
            Shape of the histogram-bins
        normed
            If True, then plot histogram will be normalized
    Return Values
        mean
            Binned mean density
        sigma
            Standard deviation of the mean
    """
    # Eric Andersson, 21-07-2017
    import numpy as np
    import PencilCode as pc
    from scipy.integrate import simps

    # Read data
    dim = pc.read.dimensions(datadir = datadir)
    ndim = (dim.nxgrid > 1) + (dim.nygrid > 1) + (dim.nzgrid > 1)
    if ndim == 3:
        raise ExceedDimensionsError(
                "Works only for 2D simulations.")

    t, rhop = pc.read.slices('rhop', datadir = datadir)

    # Allocate memory
    Nrhop = np.zeros((t.size, bins.size - 1))
    print(
    'Computes the average from t = ' + str(t[t0]) + ' to end of simulation.')
    
    for i in range(t0, t.size):
        Nrhop[i - t0][:] = np.histogram(np.concatenate(
                        rhop[i][plane]),
                                bins = bins, normed=normed)[0]

    mean = simps(y=Nrhop, x=t, axis=0)/(t[-1] - t[t0])
    sigma2 = simps(y=Nrhop**2, x=t, axis=0)/(t[-1] - t[t0]) - mean**2

    return mean, np.sqrt(sigma2)
