# Module containing wrappers for plotting particular spectra
# Stephanie Douglas, 31 December 2013
################################################################################
"""
All functions take:
  ax - matplotlib axis instance
  spectrum - a dictionary containing 'wavelength', 'flux', and 'unc'
  plot_noise (bool) - whether to plot the uncertainty spectrum or not

functions:
    plot_low
"""


import matplotlib.pyplot as plt
import numpy as np

def plot_low(ax,spectrum,plot_noise=False):
    """
    """

    ax.step(spectrum['wavelength'],spectrum['flux'],'k-')
    if plot_noise:
        ax.step(spectrum['wavelength'],spectrum['unc'],'-',color='Grey')
