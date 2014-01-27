# Just testing that the 'agg' backend outputs things properly
# It's using the master branch with ndim*2 walkers and ndim*10 steps
import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model, bdmcmc.sample
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/dustylow.pkl',wave_unit=u.AA)

bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,am.params)

bdsamp.mcmc_go()

bdsamp.plot_chains()
plt.savefig('test_agg_ch_2014-01-21.png')

bdsamp.plot_triangle()
plt.savefig('test_agg_tri_2014-01-21.png')
