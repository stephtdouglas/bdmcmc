# Module for making a BrownDwarf object including all its spectra
# 2 December 2013, Stephanie Douglas
################################################################################


import numpy as np
from config import *

class BrownDwarf(object):
    """
    """


    def __init__(self,name):
        if name[0]=='U':
            query_result = db.query.execute(
                "SELECT id FROM sources WHERE unum='{}'".format(
                name)).fetchall()
        elif ((len(name)==7) & ((name[4]=='+') | (name[4]=='-'))):
            query_result = db.query.execute(
                "SELECT id FROM sources WHERE shortname='{}'".format(
                name)).fetchall()
        else:
            query_result = db.query.execute(
                "SELECT id FROM sources WHERE names='{}'".format(
                name)).fetchall()
        
        if len(query_result)==1:
            self.sid = query_result[0][0]
        else:
            print 'Object {} not found!'.format(name)
            self.sid = np.inf


    def get_low(self,obs_date=None,filename=None):
        """
        Get low-resolution (SpeX) observation from database

        Inputs:
            obs_date (optional) (required if there are more than 2)
                format as 'YYYY-MM-DD'
            filename (optional if there are more than 2 w/ same obs_date)
                IF THERE ARE TWO OBSERVATIONS WITH THE SAME DATE IT WILL 
                PICK THE FIRST BY DEFAULT.

        Stores wavelength, flux, and noise as arrays in the variables
        self.low_wave, self.low_flux, self.low_noise
        """

        spex_id = 6
        irtf_id = 7

        query_add=''
        if obs_date!=None:
            query_add = " obs_date='{}' AND ".format(obs_date)
        if filename!=None:
            query_add = query_add + " filename='{}' AND ".format(filename)
        print query_add

        base_query = "SELECT wavelength, flux, unc FROM spectra WHERE "
        print base_query
        end_query = " telescope_id={} AND instrument_id={}".format(
            irtf_id,spex_id)
        print end_query
        full_query = base_query+query_add+end_query
        print full_query

        query_result = db.dict.execute(full_query).fetchall()

        self.low_wave = query_result[0]['wavelength']
        self.low_flux = query_result[0]['flux']
        self.low_noise = query_result[0]['unc']

