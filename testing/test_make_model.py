import logging

import astropy.units as u

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model

reload(bdmcmc.make_model)

logging.basicConfig(level=logging.DEBUG)

bd = bdmcmc.spectra.BrownDwarf('U20012')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/modelSpectra/SpeX_marley.pkl')

mg = bdmcmc.make_model.ModelGrid(bd.specs['low'],am.model,am.params)

test_model_p = [4.794,2.34,2139.]
test_norm_p = [1e-16,2e-16,1.4e-16]
test_s = 5e-16
test_p = np.append(np.append(test_model_p,test_norm_p),test_s)

print mg(test_p)
