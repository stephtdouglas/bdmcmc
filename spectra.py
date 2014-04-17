# Module for getting spectra from the database
# Includes the BrownDwarf class, for tying an object to its spectra
# 2 December 2013, Stephanie Douglas
################################################################################

import logging

import numpy as np
import astropy.units as u

# config loads database and makes it available as db
from config import * 

def spectrum_query(source_id,telescope_id,instrument_id,return_header=False,
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

        header (bool; default=False) whether to return the header

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
    #logging.debug(query_add)

    base_query = ("SELECT wavelength_units, flux_units,"
        +" wavelength, flux, unc, header FROM spectra WHERE ")
    #logging.debug(base_query)
    if telescope_id=='' and instrument_id=='':
        end_query = " source_id={}".format(source_id)
    else:
        end_query = (" source_id={} AND telescope_id={} AND instrument_id={}"
            ).format(source_id,telescope_id,instrument_id)
    #logging.debug(end_query)
    full_query = base_query+query_add+end_query
    logging.debug(full_query)

    q_result = db.dict.execute(full_query).fetchall()
    logging.debug(q_result)
    #print '0',q_result[0]

    try:
        q_result = q_result[0]
    except:
        return {'wavelength':[],'flux':[],'unc':[]}
        logging.info('spectrum not found for %s', full_query)

    wave_unit = q_result['wavelength_units']
    flux_unit = q_result['flux_units']
    flux_unit = flux_unit.replace('ergs','erg s').replace('Wm','W m')
    flux_unit = u.Unit(flux_unit.replace('normalized',''))
    wave_unit = u.Unit(wave_unit)
    logging.info('%s %s', wave_unit.to_string('fits'), 
        flux_unit.to_string('fits'))

    q_result_out = {'wavelength':q_result['wavelength']*wave_unit,
         'flux':q_result['flux']*flux_unit}
    if q_result['unc']==None:
         q_result_out['unc'] = np.ones(len(q_result['flux']))*flux_unit
         logging.info('no unc array found for %s', full_query)
    else:
         q_result_out['unc'] = q_result['unc']*flux_unit

    if return_header:
        q_result_out['header'] = q_result['header']

    logging.debug(str(q_result_out.keys()))
    return q_result_out

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
        """
        Parameters
        ----------
        name: string
            identifier from the database
            -if the first character is 'U', assumes it's a u-number
            -if it is 7 characters long with a + or - at the 5th
                spot, assumes it's a shortname
            -otherwise, it will search under names

        Creates
        -------
        self.name: string
            equals the input name

        self.sid: integer
            source_id for further database queries

        self.specs: dictionary
            initializes an empty dictionary to hold spectra
        """
        if name[0]=='U':
            query_result = db.query.execute(
                "SELECT id, unum, shortname FROM sources WHERE unum='{}'"
                .format(name)).fetchall()
            logging.info('using unum %s', name)
        elif ((len(name)==9) & ((name[4]=='+') | (name[4]=='-'))):
            query_result = db.query.execute(
                "SELECT id, unum, shortname FROM sources WHERE shortname='{}'"
                .format(
                name)).fetchall()
            logging.info('using shortname %s', name)
        else:
            query_result = db.query.execute(
                "SELECT id, unum, shortname FROM sources WHERE names='{}'"
                .format(name)).fetchall()
            logging.info('using names %s',name)
        if len(query_result)==1:
            self.sid = query_result[0][0]
            self.unum = query_result[0][1]
            self.shortname = query_result[0][2]
            logging.info('source_id %s',str(self.sid))
        else:
            logging.info('Object %s not found!',name)
            self.sid = np.inf
            self.unum = ''
            self.shortname = ''
        self.name = name

        self.get_spt()

        self.specs = {}



    def get_low(self,obs_date=None,filename=None):
        """
        Get low-resolution (SpeX) observation from database;
        calls spectra.spectrum_query

        Parameters
        ----------
        obs_date: string (default=None) 
            (required if there are 2 or more possibilities)
            format as 'YYYY-MM-DD'
        filename: string (default=None) 
            (required if there are more than 2 w/ same obs_date)
            IF THERE ARE TWO OBSERVATIONS WITH THE SAME DATE IT WILL 
            PICK THE FIRST BY DEFAULT.

        Creates
        -------
        self.specs['low']: dictionary
            'wavelength', 'flux', and 'unc' as arrays,
            and 'slit_width' as a float, all in a nested dictionary
            access as, e.g., self.specs['low']['wavelength']
        """

        spex_id = 6
        irtf_id = 7

        self.specs['low'] = spectrum_query(self.sid,
            irtf_id, spex_id, True, obs_date, filename)

        if 'header' in self.specs['low'].keys():
            header = self.specs['low'].pop('header')
            self.specs['low']['slit_width'] = header['SLTW_ARC']*u.arcsec
        else:
            self.specs['low']['slit_width'] = -9.
            logging.info('no header available; no slit_width information')



    def get_order(self,order,obs_date=None):
        """
        Get one order of a high_resolution spectrum (NIRSPEC) 
        from the database; calls spectra.spectrum_query

        Parameters
        ----------
        order: integer 
            number in range 58-65 for desired NIRSPEC order
        obs_date: string (default=None) 
            (required if there are 2 or more possibilities)
            format as 'DDmmmYY', e.g. '04dec08'

        Creates
        -------
        self.specs[##]: dictionary
            'wavelength', 'flux', and 'unc' as arrays in a nested dictionary
            within self.specs, keyed by the order number
            access as, e.g., self.specs[61]['wavelength']

        """

        inst_id = 9 # NIRSPEC
        scope_id = 9 # Keck II

        self.specs[order] = spectrum_query(self.sid,scope_id,inst_id,
            obs_date=obs_date,order=order)


    def get_6165(self,obs_date=None):
        """
        Get orders 61 and 65 from the database, which adds them to
        self.specs, then join them into one spectrum

        Parameters
        ----------
        obs_date: string (default=None) 
            (required if there are 2 or more possibilities)
            format as  'DDmmmYY', e.g. '04dec08'

        Creates
        -------
        self.specs[6165]: dictionary
            'wavelength', 'flux', and 'unc' as arrays in a nested dictionary
            access as, e.g., self.specs[6165]['wavelength']
        """

        self.get_order(61,obs_date)
        self.get_order(65,obs_date)

        self.specs[6165] = {'wavelength':np.append(
            self.specs[65]['wavelength'],self.specs[61]['wavelength']),
            'flux':np.append(self.specs[65]['flux'],
            self.specs[61]['flux']),
            'unc':np.append(self.specs[65]['unc'],self.specs[61]['unc'])}


    def get_spt(self):
        spts = db.query.execute(("SELECT spectral_type, regime FROM "+
            "spectral_types WHERE source_id={}".format(self.sid))
            ).fetchall()
        if len(spts)>0:
            self.spt = spts[0][0]
            self.spt_regime = spts[0][1]
            for x in spts:
                if x[1]=='OPT':
                    self.spt,self.spt_regime = x[0],x[1]
        else:
            self.spt, self.spt_regime = -99, None

        if type(self.spt)==str:
            self.spt=float(self.spt[0:4])
