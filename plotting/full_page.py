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

def page_plot(chains):

    K=3 # ndim of model, obviously a placeholder

    burn_in = 0
    cropchain = chains[:,burn_in:,:].reshape(
            (-1,K))

    # set up figure and axes
    plt.figure(figsize=(8,10))
    # split the figure into 2 parts, top and bottom
    split_grid = gridspec.GridSpec(2,1) 

    # the top grid will contain the model and random draws
    rand_grid = gridspec.GridSpecFromSubplotSpec(1,1,
        subplot_spec=split_grid[0])
    rand_ax = plt.subplot(top_grid[:,:])

    # the bottom grid will contain the corner plot
    corner_grid = gridspec.GridSpecFromSubplotSpec(K,K,
        subplot_spec=split_grid[1])

    # plot the random samples plot in the top axis
    random_sample = plot_random(cropchain)


    # plot the corner plot



    # overplot the random samples
    # OR do I want to let triangle calculate random samples
    # and overplot them, then return those and plot them?
    # that seems easier and requires less axis manipulation
