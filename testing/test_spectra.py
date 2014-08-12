import logging, os

import numpy as np
import astropy.units as u

import bdmcmc.spectra, bdmcmc.sample.select_sxd

reload(bdmcmc.spectra)
reload(bdmcmc.sample.select_sxd)

logging.basicConfig(level=logging.DEBUG)

bd = bdmcmc.spectra.BrownDwarf('U20012')
bd.get_low()

# check the SXD sample
spex_sxd = bdmcmc.sample.select_sxd.fetch_sxd()

for name,date,sid in spex_sxd:
     if (name!=None) and (name!='0559-1404') and (name!='1254-0122'):
         bd = bdmcmc.spectra.BrownDwarf(name)
         bd.get_low(obs_date=date)
         print name,date,bd.specs['low']['slit_width']
