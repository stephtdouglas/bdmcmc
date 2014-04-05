# Create a full-page plot with 
#    -spectrum with 200 random draws
#    -corner plot with locations of random draws included


import datetime
import logging

import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np

from bdmcmc.plotting.plot_random import plot_random
from bdmcmc.plotting import triangle

def page_plot(chains,model):
    """
    Plot data plus random sample of models, plus corner plot of run

    Parameters
    ----------
    chains : 

    model : ModelGrid instance

    """

    K=3 # ndim of model, obviously a placeholder

    burn_in = 0
    cropchain = chains[:,burn_in:,:].reshape(
            (-1,K))

    # set up figure and axes
    plt.figure(figsize=(8,10))
    # split the figure into 2 parts, top and bottom
    split_grid = gridspec.GridSpec(2,1,height_ratios=[1,2]) 

    # the top grid will contain the model and random draws
    rand_grid = gridspec.GridSpecFromSubplotSpec(1,1,
        subplot_spec=split_grid[0])
    rand_ax = plt.subplot(rand_grid[:,:])

    # the bottom grid will contain the corner plot
    corner_grid = gridspec.GridSpecFromSubplotSpec(K,K,
        subplot_spec=split_grid[1])

    # plot the random samples plot in the top axis
    random_samp = plot_random(cropchain,model,ax=rand_ax)
    #rand_ax.


    # plot the corner plot
    fig,axes_array = triangle.corner(cropchain,spec_grid=corner_grid,
        labels=model.params, quantiles=[0.16,0.5,0.84])


    # overplot the random samples
    for i in range(K):
        for j in np.arange(i):
            ax = axes_array[i,j]

            if (i==(K-1)) or (j==(K-1)):
                plot_color='DarkOrange'
            else:
                plot_color='r'
            ax.plot(random_samp[:,j],random_samp[:,i],'.',
                color=plot_color,alpha=0.5)



    factor = 2.0           # size of one side of one panel
    lbdim = 0.5 * factor   # size of left/bottom margin
    trdim = 0.1  * factor  # size of top/right margin
    whspace = 0.075         # w/hspace size
    plotdim = factor * K + factor * (K - 3.) * whspace
    dim = lbdim + plotdim + trdim

    lb = lbdim / dim
    tr = (lbdim + plotdim) / dim
    plt.subplots_adjust(left=lb, bottom=lb, right=tr, top=0.95,
                        wspace=whspace, hspace=whspace)
