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

def quantile(x,quantiles):
    xsorted = sorted(x)
    qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
    return zip(quantiles,qvalues)

def unc_string(quant_name,quant_array):
    error_string = '{0}={1:.2f}\n         +{2:.2f}/-{3:.2f}'.format(quant_name,
        quant_array[1][1],quant_array[2][1]-quant_array[1][1],
        quant_array[1][1]-quant_array[0][1])
    return error_string


def page_plot(chains,model,plot_title,extents=None):
    """
    Plot data plus random sample of models, plus corner plot of run

    Parameters
    ----------
    chains : array_like (nwalkers, nsteps, ndim)
         output from emcee run

    model : ModelGrid instance

    plot_title : string

    extents : iterable (ndim,) (optional)
        A list of length 2 tuples containing lower and upper bounds
        for each dimension, e.g., [(0.,10.), (1.,5.),etc]

    """

    logging.debug(str(np.shape(chains)))
    K = np.shape(chains)[-1]

    burn_in = 0
    cropchain = chains[:,burn_in:,:].reshape((-1,K))
#    cropchain = chains.reshape((-1,K))

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
    rand_ax.set_title(plot_title)


    # plot the corner plot
    logging.debug(len(model.params))
    if len(model.params)!=K: # need to make a better check on this
        labels = list(np.append(model.params,'ln(s)'))
    else:
        labels=model.params
    logging.info(labels)
    fig,axes_array = triangle.corner(cropchain,spec_grid=corner_grid,
        labels=labels, quantiles=[0.16,0.5,0.84],extents=extents)


    # overplot the random samples
    for i in range(K):
        for j in np.arange(i):
            ax = axes_array[i,j]

            if (i==(K-1)) or (j==(K-1)):
                plot_color='DarkOrange'
            else:
                plot_color='r'
            ax.plot(random_samp[:,j],random_samp[:,i],'.',
                color=plot_color,alpha=0.5,mec='None')

    ax = axes_array[0,2]
    ax.set_xlim(0,1)
    ax.set_ylim(0,1)
    text_y_pos = 0.6
    for i in range(K):
        quantiles = quantile(cropchain[:,i],[.16,.5,.84])
        string_unc = unc_string(labels[i],quantiles)
        ax.text(0.0,text_y_pos,string_unc)
        text_y_pos -= 0.4
        ax.set_visible(True)
        ax.tick_params(which='both',length=0,width=0,labelleft=False,
            labelbottom=False,labelright=False,labeltop=False)

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

    return rand_ax,axes_array
