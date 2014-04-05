# coding: utf-8

""" Plotting functions for visualizing emcee output. """

from __future__ import division, print_function

__author__ = "adrn <adrn@astro.columbia.edu>"

# Standard library
import os, sys

# Third-party
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import astropy.units as u
from astropy.io.misc import fnunpickle

def emcee_plot(xs, labels=None, truths=None, extents=None, fig=None):
    """ Plot posterior probability distributions and chain traces from the 
        chains associated with the given sampler. 
        
        Parameters
        ----------
        xs : array_like (nwalkers, nsamples, ndim)
            The samples. This should be a 3-dimensional array. 
        labels : iterable (ndim,) (optional)
            A list of names for the dimensions.        
        truths : iterable (ndim,) (optional)
            A list of reference values to indicate on the plots.
        extents : iterable (ndim,) (optional)
            A list of length 2 tuples containing lower and upper bounds (extents)
            for each dimension, e.g., [(0.,10.), (1.,5), etc.]
    """
    
    nwalkers, nsamples, ndim = xs.shape
    
    if fig is None:
        fig = plt.figure(figsize=(8,10))
    
    # I want the trace plots to span two columns, the histograms one column
    gs = gridspec.GridSpec(ndim, 3)
            
    # For each parameter, I want to plot each walker on one panel, and a histogram
    #   of all links from all walkers
    for ii in range(ndim):
        walkers = xs[:,:,ii]
        flatchain = np.hstack(walkers)
        ax1 = plt.subplot(gs[ii, :2])
        
        steps = np.arange(nsamples)
        for walker in walkers:
            ax1.plot(steps, walker,
                     drawstyle="steps", color="#555555", alpha=0.5)
        
        if labels:
            ax1.set_ylabel(labels[ii], fontsize=36, labelpad=18,
                           rotation="horizontal", color="k")
        
        # Don't show ticks on the y-axis
        ax1.yaxis.set_ticks([])
        
        # For the plot on the bottom, add an x-axis label. Hide all others
        if ii == ndim-1:
            ax1.set_xlabel("step number", fontsize=24, labelpad=18, color="k")
        else:
            ax1.xaxis.set_visible(False)
        
        ax2 = plt.subplot(gs[ii, 2])
        
        # Same y-bounds as the walkers plot, so they line up
        #ax1.set_ylim(np.min(these_chains[:,ii]), np.max(these_chains[:,ii]))
        if extents:
            ax1.set_ylim(extents[ii])
        else:
            mu,sigma = np.median(flatchain), np.std(flatchain)
            ax1.set_ylim(mu-10*sigma, mu+10*sigma)
            
        ax2.set_ylim(ax1.get_ylim())
        ax2.xaxis.set_visible(False)
        ax2.yaxis.tick_right()
        
        # Create a histogram of all samples. Make 100 bins between the y-axis 
        #   bounds defined by the 'walkers' plot.
        ax2.hist(flatchain, 
                 orientation='horizontal',
                 bins=np.linspace(ax1.get_ylim()[0],ax1.get_ylim()[1],100),
                 facecolor="#67A9CF",
                 edgecolor="none")
        
        if truths:
            ax2.axhline(truths[ii], color="#016C59", linestyle="--")
        
        # For the first plot, add titles and shift them up a bit
        if ii == 0:
            t = ax1.set_title("Walkers", fontsize=30, color="k")
            t.set_y(1.01) 
            t = ax2.set_title("Posterior", fontsize=30, color="k")
            t.set_y(1.01) 
        
        # Adjust axis ticks, e.g. make them appear on the outside of the plots and
        #   change the padding / color.
        ax1.tick_params(axis='x', pad=2, direction='out', colors="#444444", labelsize=14)
        ax2.tick_params(axis='y', pad=2, direction='out', colors="#444444", labelsize=14)
        
        # Removes the top tick marks
        ax1.get_xaxis().tick_bottom()
    
    fig.subplots_adjust(hspace=0.02, wspace=0.0, bottom=0.075, top=0.9, left=0.12, right=0.88)
    return fig
