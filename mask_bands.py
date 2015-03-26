# Module to mask out regions of the fit
# Stephanie Douglas, 28 February 2014
################################################################################

import logging

import numpy as np
from astropy import units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

class BandMask(object):
    """
    Class containing wavelength and pixel values for masking regions
    of a spectrum from the fit/showing in a plot.
    """

    def __init__(self,wavelength_grid=None):
        """
        Parameters
        ----------
        wavelength_grid: astropy.units Quanitity array (optional)
            if provided, then pixel values for the mask will be created
            otherwise, only wavelength values will be stored
    
        """

        self.wave_mask = np.array([])

        if wavelength_grid is not None:
            self.wave = wavelength_grid

    def mask_Hband(self):
        self.wave_mask = np.append(self.wave_mask,[1.58,1.75]).reshape((-1,2))

    def mask_FeH(self):
        self.mask_Hband()
        self.wave_mask = np.append(self.wave_mask,[0.99,1.01]).reshape((-1,2))
        self.wave_mask = np.append(self.wave_mask,[1.2,1.33]).reshape((-1,2))
        
    def plot_mask(self,ax=None):
        if ax==None:
            ax = plt.gca()

        ylims = ax.get_ylim()

        for wave_range in self.wave_mask:
            ax.add_patch(Rectangle((wave_range[0],ylims[0]),(wave_range[1]-
                wave_range[0]),(ylims[1]-ylims[0]),
                fc='#DCDCDC',ec='none',fill=True))

    def make_pixel_mask(self):
        self.pixel_mask = np.array([],int)
        for region in self.wave_mask:
            new_pixels = np.where((self.wave>=region[0]*u.um) & (self.wave<=region[1]*u.um))[0]
            self.pixel_mask = np.append(self.pixel_mask,new_pixels)
