import logging, os

import numpy as np
import astropy.units as u

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model

reload(bdmcmc.spectra)
reload(bdmcmc.make_model)

logging.basicConfig(level=logging.DEBUG)

mbase_path = '/vega/astro/users/sd2706/'
if os.path.exists(mbase_path)==False:
    mbase_path = '/home/stephanie/ldwarfs/'

bd = bdmcmc.spectra.BrownDwarf('U20012')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel(mbase_path+'modelSpectra/SpeX_marley_nolowg.pkl')

test_model_p = [4.794,2.34,2139.]
test_norm_p = [1e-16,2e-16,1.4e-16]
test_s = 5e-16
test_p = np.append(np.append(test_model_p,test_norm_p),test_s)

logging.debug("testing without snap")
mg1 = bdmcmc.make_model.ModelGrid(bd.specs['low'],am.model,am.params)
print mg1(test_p)

am2 = bdmcmc.get_mod.AtmoModel(mbase_path+'modelSpectra/SpeX_marley_nolowg.pkl')
mg3 = bdmcmc.make_model.ModelGrid(bd.specs['low'],am2.model,am2.params)

# Test find_nearest2
num_models = len(mg3.plims[mg3.params[0]]["vals"])
param_arrays = [[mg3.plims[mg3.params[i]]["vals"][j] for i in range(mg3.ndim)]
                for j in range(num_models)]
mg3.find_nearest2(param_arrays,test_model_p)
mg3.check_grid_coverage()

logging.debug("testing WITH snap")
mg2 = bdmcmc.make_model.ModelGrid(bd.specs['low'],am.model,am.params,snap=True)
print mg2(test_p)
