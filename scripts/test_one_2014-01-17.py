
# this one is using the ignore_bands branch to test that branch
# It's also using ndim*100 walkers and ndim*50 steps
import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
import astropy.units as u
import matplotlib.pyplot as plt

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/summerAMNH/modelSpectra/dustylow.pkl',wave_unit=u.AA)

bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params)

bdsamp.mcmc_go()

bdsamp.plot_chains()
plt.savefig('test_one_2014-01-17.png')

bdsamp.plot_triangle()
plt.savefig('test_one_tri_2014-01-17.png')
