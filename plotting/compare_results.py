#!/usr/bin/env python
# -*- coding: utf-8 -*-

#from __future__ import print_function, absolute_import, unicode_literals

__all__ = ["corner","plot_two"]
__version__ = "0.0.1"
__author__ = "Stephanie Douglas (sdouglas@astro.columbia.edu)"
__contributors__ = [
    # Alphabetical by first name.
    "Adrian Price-Whelan @adrn",
    "Brendon Brewer @eggplantbren",
    "Stephanie Douglas @stephtdouglas",
    "Dan Foreman-Mackey @exoplaneteer"
    "Ekta Patel @ekta1224",
    "Emily Rice @emilurice",
    "Geoff Ryan @geoffryan",
    "Phil Marshall @drphilmarshall",
    "Pierre Gratier @pirg",
]

import logging

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as pl
from matplotlib.ticker import MaxNLocator
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.patches import Ellipse
import matplotlib.gridspec as gridspec
import matplotlib.cm as cm
from matplotlib.cbook import flatten


def corner(medians, errors, spt, param_labels, run_labels,
           extents=None, truths=None, truth_color="#4682b4",
           scale_hist=False, quantiles=[], verbose=True, 
           spec_grid=None, **kwargs):
    """
    Make a *sick* corner plot showing the projections of a data set in a
    multi-dimensional space. kwargs are passed to hist2d() or used for
    `matplotlib` styling.

    Parameters
    ----------
    medians: array_like (ndim, num_runs, num_samples)
        Results of a set of mcmc runs.  There should be one row for each parameter,
        and then moving across each row gives the median value from a given run.

    errors: array_like (ndim, num_runs, 2, num_samples)
        errors associated with the median values given above. Assymetric error
        bars are expected.

    spts: array_like (num_samples)
        spectral type of input object

    param_labels : iterable (ndim,) 
        A list of names for the dimensions.

    extents : iterable (ndim,) (optional)
        A list of length 2 tuples containing lower and upper bounds (extents)
        for each dimension, e.g., [(0.,10.), (1.,5), etc.]

    truths : iterable (ndim,) (optional)
        A list of reference values to indicate on the plots.

    truth_color : str (optional)
        A ``matplotlib`` style color for the ``truths`` makers.

    scale_hist : bool (optional)
        Should the 1-D histograms be scaled in such a way that the zero line
        is visible?

    quantiles : iterable (optional)
        A list of fractional quantiles to show on the 1-D histograms as
        vertical dashed lines.

    verbose : bool (optional)
        If true, print the values of the computed quantiles.

    spec_grids : matplotlib.gridspec.GridSpec (optional)
        array of 

    """

    K = len(medians)
    logging.debug(str(medians))
    logging.debug(str(errors))
    logging.debug('K {} spt {}'.format(K,spt))
    factor = 2.0           # size of one side of one panel
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.25 * factor  # size of top/right margin
    whspace = 0.1         # w/hspace size
    plotdim = factor * K + factor * (K - 1.) * whspace
    dim = lbdim + plotdim + trdim

    if spec_grid is None:
        fig = pl.figure(figsize=(8,10))
        spec_grid = gridspec.GridSpec(K,K)
    else:
        fig=pl.gcf()

    #set up a full grid for ease; unused spots will be set to invisible later
    setup_axes = [[pl.subplot(spec_grid[i,j]) for j in np.arange(K)] 
                  for i in np.arange(K)]
    axes = np.array(setup_axes).reshape((K,K))

    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    pl.subplots_adjust(left=lb, bottom=lb, right=tr, top=tr,
                        wspace=whspace, hspace=whspace)

    if extents is None:
        extents = [[min(flatten(x)), max(flatten(x))] for x in medians]
        for i in range(K):
            if param_labels[i]=='teff':
                extents[i]= [1400,2400]
            elif param_labels[i]=='logg':
                extents[i] = [3.0,5.8]

        # Check for parameters that never change.
        m = np.array([e[0] == e[1] for e in extents], dtype=bool)
        if np.any(m):
            raise ValueError(("It looks like the parameter(s) in column(s) "
                              "{0} have no dynamic range. Please provide an "
                              "`extent` argument.")
                             .format(", ".join(map("{0}".format,
                                                   np.arange(len(m))[m]))))

    for i in range(K):
        ax = axes[i, i]

        # Plot ith quantity vs. SpT, with comparison as needed
        logging.debug('i {} {}'.format(i,type(i)))
        logging.debug(str(medians[i]))
        single_plot_setup(param_labels[i],ax)
        logging.debug('len {}'.format(len(medians[i])))
        mshape = np.shape(medians[i])
        logging.debug('shape {} {}'.format(mshape, sum(mshape)))
        spt_one = np.ones(mshape)*spt
        #spt_two = np.arange(sum(mshape)).reshape(mshape)*0.1
        spt_array = spt_one #+ spt_two
        spt_errors = np.zeros(np.shape(errors[i]))
        plot_two(medians[i],spt_array,spt_errors,errors[i],run_labels,ax=ax)
        if i==0:
            ax.legend(loc=3,numpoints=1)


        #if truths is not None:
        #    ax.axvline(truths[i], color=truth_color)

        # Set up the axes.
        ax.set_xlim(extents[i])

        # Not so DRY.
        if i < K - 1:
            ax.set_xticklabels([]) 
        else:
            [l.set_rotation(45) for l in ax.get_xticklabels()]
            if param_labels is not None:
                ax.set_xlabel(param_labels[i])
#                ax.xaxis.set_label_coords(0.5, -0.3)

        for j in range(K):
            print 'j',j
            ax = axes[i, j]
            if j > i:
                ax.set_visible(False)
                ax.set_frame_on(False)
                continue
            elif j == i:
                continue

            # call plot_two
            print i, param_labels[i], j, param_labels[j]
            plot_two(medians[j],medians[i],errors[i],errors[j],run_labels,ax=ax)
            ax.set_xlim(extents[j])
            ax.set_ylim(extents[i])


            #if truths is not None:
            #    ax.plot(truths[j], truths[i], "s", color=truth_color)
            #    ax.axvline(truths[j], color=truth_color)
            #    ax.axhline(truths[i], color=truth_color)

#            ax.xaxis.set_major_locator(MaxNLocator(5))
#            ax.yaxis.set_major_locator(MaxNLocator(5))

            if i < K - 1:
                ax.set_xticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_xticklabels()]
                if param_labels is not None:
                    ax.set_xlabel(param_labels[j])
#                    ax.xaxis.set_label_coords(0.5, -0.3)

            if j > 0:
                ax.set_yticklabels([])
            else:
                [l.set_rotation(45) for l in ax.get_yticklabels()]
                if param_labels is not None:
                    ax.set_ylabel(param_labels[i])
#                    ax.yaxis.set_label_coords(-0.3, 0.5)

    return fig, axes

def plot_two(x,y,yerr,xerr,run_labels,*args,**kwargs):
    """
    Plot a set of discrete values

    x : arraylike (1D)
        length num_samples

    y : arraylike (1D)
        length num_samples

    yerr : arraylike (1D or 2D)
        if 2D, shape must be (num_samples,2)

    xerr : arraylike (1D or 2D)
        if 2D, shape must be (num_samples,2)

    run_labels : array of strings
        names of individual runs to be compared

    """
    ax = kwargs.pop("ax", pl.gca())
    extent = kwargs.pop("extent", None)


    num_runs = len(x)

    markers = ['o','v','s','^','*','<','D','>']
    cmap = cm.get_cmap("Paired")
    color_norm = Normalize(vmin=0,vmax=7)
    scalar_map = cm.ScalarMappable(norm=color_norm,cmap=cmap)
    

    logging.debug('{} runs'.format(num_runs))
    for i in range(num_runs):
        plot_color = scalar_map.to_rgba(i)
        logging.debug('x {}'.format(x[i]))
        logging.debug('y {}'.format(y[i]))
        logging.debug('xerr {}'.format(xerr[i].reshape((-1,2)).T))
        logging.debug('yerr {}'.format(yerr[i].reshape((-1,2)).T))
        ax.errorbar(x[i],y[i],
            yerr[i].reshape((-1,2)).T,xerr[i].reshape((-1,2)).T,
            marker=markers[i],color=plot_color,ms=9,mec='None',elinewidth=2,
            capsize=5,ecolor=plot_color,barsabove=True,mew=2,
            label=run_labels[i],linewidth=0)


def single_plot_setup(param_name,ax):
    """
    """
    logging.debug(param_name)

    spt = np.arange(6,22) # M6=6, L0=10, T0=21

    if param_name.lower()=='teff':
        teff = (2265.9 + 347.82*spt - 60.558*spt**2 + 3.151*spt**3
            -0.060481*spt**4 + 0.00024506*spt**5) #eqn 3 in Stephens 2009
        ax.plot(teff,spt,color='DarkGrey',linewidth=6)
        ax.text(1450,19,"Stephens+2009 (eqn 3)",color='#505050',fontsize='small')
    elif param_name.lower()=='logg':
        base_ones = np.ones(len(spt))
        for i in np.array([5,4,3.5]):
            ax.plot(base_ones*i,spt,color='DarkGrey',linewidth=3)
#    else:
#        print "unknown parameter name; not plotting"

    new_yticks = np.arange(6,24,2)
    ax.set_yticks(new_yticks)
    new_ylabels = np.empty(len(new_yticks),'S2')
    m_loc = (new_yticks<10)
    new_ylabels[m_loc] = ['M{}'.format(s) for s in new_yticks[m_loc]]
    l_loc = ((new_yticks<20) & (new_yticks>=10))
    logging.debug('spt ticks {}'.format(new_yticks[l_loc]))
    new_ylabels[l_loc] = ['L{}'.format(s) for s in (new_yticks[l_loc]-10)]
    t_loc = (new_yticks>=20)
    new_ylabels[t_loc] = ['T{}'.format(s) for s in (new_yticks[t_loc]-20)]
    ax.set_yticklabels(new_ylabels)
    ax.tick_params(labelleft=False,labelright=True)
    ax.set_ylim(20,9)
