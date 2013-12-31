# Module for getting spectra from the database
# Includes the BrownDwarf class, for tying an object to its spectra
# 2 December 2013, Stephanie Douglas
################################################################################


import numpy as np

# config loads database and makes it available as db
from config import * 

def spectrum_query(source_id,telescope_id,instrument_id,
    obs_date=None,filename=None,order=None):
    """
    Inputs:
        telescope_id (int) 
          [(7, 'NASA IRTF'),(8, 'Keck I'),(9, 'Keck II'),(10, 'KP 4m'),
          (11, 'KP 2.1m'),(12, 'KP Bok'),(13, 'MMT'),(14, 'CTIO 1.5m'),
          (15, 'CTIO 4m'),(16, 'Gemini-North'),(17, 'Gemini-South'),
          (18, 'ESO VLT U2'),(19, 'ARC 3.5m'),(20, 'Subaru')]

        instrument_id (int)
          [(1, 'R-C Spec'),(2, 'GMOS-N'),(3, 'GMOS-S'),(4, 'FORS'),
          (5, 'LRIS'),(6, 'SPeX, IRTF Spectrograph'),(7, 'LDSS3-Two'),
          (8, 'FOCAS'),(9, 'NIRSPEC')]

        obs_date (optional) (required if there are more than 2)
            format as 'YYYY-MM-DD'
        filename (optional if there are more than 2 w/ same obs_date)
            IF THERE ARE TWO OBSERVATIONS WITH THE SAME DATE IT WILL 
            PICK THE FIRST BY DEFAULT.
        order (optional) order of high-resolution spectrum

    Outputs:
        result - a dictionary containing 'wavelength','flux', & 'noise'
    """
    query_add=''
    if obs_date!=None:
        query_add = " obs_date='{}' AND ".format(obs_date)
    if filename!=None:
        query_add = query_add + " filename='{}' AND ".format(filename)
    if order!=None:
        query_add = query_add + " wavelength_order='{}' AND ".format(
            order)
    #print query_add

    base_query = "SELECT wavelength, flux, unc FROM spectra WHERE "
    #print base_query
    end_query = " source_id={} AND telescope_id={} AND instrument_id={}".format(
        source_id,telescope_id,instrument_id)
    #print end_query
    full_query = base_query+query_add+end_query
    #print full_query

    q_result = db.dict.execute(full_query).fetchall()

    try:
        return q_result[0]
    except:
        return {'wavelength':[],'flux':[],'unc':[]}


class BrownDwarf(object):
    """
    Creates a brown dwarf object, containing its db id number
    plus spectra if desired.


    Call as x = spectra.BrownDwarf(<name>) where <name> is
        a U-number, a short name, or another name
        that can be used to find the object in the database

    Functions:
        __init__
        get_low
        get_order
        get_6165

    Variables:
        self.sid (int)
        self.specs (dict)

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
        self.name = name

        self.specs = {}



    def get_low(self,obs_date=None,filename=None):
        """
        Get low-resolution (SpeX) observation from database;
        calls spectra.spectrum_query

        Inputs:
            obs_date (optional) (required if there are more than 2)
                format as 'YYYY-MM-DD'
            filename (optional if there are more than 2 w/ same obs_date)
                IF THERE ARE TWO OBSERVATIONS WITH THE SAME DATE IT WILL 
                PICK THE FIRST BY DEFAULT.

        Stores wavelength, flux, and unc as arrays in the dictionary
        self.specs, e.g. self.specs['low']['wavelength']
        """

        spex_id = 6
        irtf_id = 7

        self.specs['low'] = spectrum_query(self.sid,
            irtf_id, spex_id, obs_date, filename)
        


    def get_order(self,order,obs_date=None):
        """
        Get one order of a high_resolution spectrum (NIRSPEC) 
        from the database; calls spectra.spectrum_query

        Inputs: 
            order (int) number in range 58-65 for NIRSPEC
            obs_date (optional) in format 'DDmmmYY', e.g. '04dec08'

        Stores wavelength, flux, and unc as arrays in the dictionary
        self.specs, keyed by the order number
        e.g. self.specs[61]['wavelength']
        """

        inst_id = 9 # NIRSPEC
        scope_id = 9 # Keck II

        self.specs[order] = spectrum_query(self.sid,scope_id,inst_id,
            obs_date=obs_date,order=order)


    def get_6165(self,obs_date=None):
        """
        Get orders 61 and 65 from the database, which adds them to
        self.specs, then join them into one spectrum

        Inputs:
            obs_date (optional) in format 'DDmmmYY', e.g. '04dec08'

        Stores wavelength, flux, and unc as arrays in the dictionary
        self.specs, e.g. self.specs[6165]['wavelength']
        """

        self.get_order(61,obs_date)
        self.get_order(65,obs_date)

        self.specs[6165] = {'wavelength':np.append(
            self.specs[65]['wavelength'],self.specs[61]['wavelength']),
            'flux':np.append(self.specs[65]['flux'],
            self.specs[61]['flux']),
            'unc':np.append(self.specs[65]['unc'],self.specs[61]['unc'])}
