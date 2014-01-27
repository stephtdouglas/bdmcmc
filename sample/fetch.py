# Module to fetch my objects from the database
# Stephanie Douglas, 31 December 2013
################################################################################

# Third Party
import numpy as np
import asciitable as at
import astropy.units as u

# Mine
from bdmcmc.config import *
from bdmcmc.spectra import *


class BDList(object):
    """
    """

    def __init__(self,infile):

        sample_list = at.read(infile)
        self.unums = sample_list['unum']
        self.spts = sample_list['Type']
        self.obs_dates = sample_list['ObsDate']

        self.brown_dwarfs = {}
        for u in self.unums:
            self.brown_dwarfs[u] = BrownDwarf(u)
            self.brown_dwarfs[u].spt = self.spts[self.unums==u][0]
            self.brown_dwarfs[u].obs_date = self.obs_dates[self.unums==u][0]
            #print u, self.brown_dwarfs[u].sid

        self.fetch_specs()


    def fetch_specs(self):
        for u in self.unums:
           self.brown_dwarfs[u].get_low()
           #print u, self.brown_dwarfs[u].sid, self.brown_dwarfs[u].specs.keys()
           for i in range(58,66):
               self.brown_dwarfs[u].get_order(i,
                    obs_date=self.brown_dwarfs[u].obs_date)





def fetch_12():
    my_bds = BDList('/home/stephanie/ldwarfs/Ldwarf_sample2014.csv')

    # Use Burgasser+08's spectrum from the SpeX Prism Library
    oldfile = at.read('/home/stephanie/ldwarfs/summerAMNH/LdwarfSpectra/spex_prism_gd165b_090629.txt')
    old_spectrum = {'wavelength':oldfile['col1']*u.um,
        'flux':oldfile['col2']*u.dimensionless_unscaled,
        'unc':oldfile['col3']*u.dimensionless_unscaled}
    my_bds.brown_dwarfs['U13079'].specs['low'] = old_spectrum

    # select appropriate spectrum for GD 165 B
    source_id = my_bds.brown_dwarfs['U40039'].sid
    my_bds.brown_dwarfs['U40039'].specs['low'] = spectrum_query(source_id,
         '','',filename='spex_prism_G196-3B_U40039.fits')


    return my_bds
