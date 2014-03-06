# running a fit with the resolution-dependent-smoothed models

import logging
from datetime import date

import numpy as np
import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model
import astropy.units as u
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from matplotlib.backends.backend_pdf import PdfPages

logging.basicConfig(level=logging.INFO)

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()

#am = bdmcmc.get_mod.AtmoModel('/vega/astro/users/sd2706/modelSpectra/SpeX_dusty.pkl')
am = bdmcmc.get_mod.AtmoModel('/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl')

am.model['wsyn'] = bd.specs['low']['wavelength']
high_grav = np.where(am.model['logg']>5.55)[0]
while len(high_grav)>0:
    i = high_grav[0]
    am.model['logg'] = np.delete(am.model['logg'],i)
    am.model['teff'] = np.delete(am.model['teff'],i)
    am.model['fsyn'] = np.delete(am.model['fsyn'],i,0)
    high_grav = np.where(am.model['logg']>5.55)[0]
am.model['fsyn'] = am.model['fsyn']*u.dimensionless_unscaled

logging.debug('script lengths dw {} mw {} mf {}'.format(
    len(bd.specs['low']['wavelength']),len(am.model['wsyn']),
    len(am.model['fsyn'][2])))

mg = bdmcmc.make_model.ModelGrid(bd.specs['low'],am.model,am.params,smooth=False)

pp = PdfPages('grid_interp_{}.pdf'.format(date.isoformat(date.today())))

teffs = np.arange(1600,2400,100)
logg = np.arange(3.5,5.5,1)


plt.figure(figsize=(10,7.5))
ax = plt.subplot(111)

my_cmap = plt.get_cmap('gist_rainbow')
c_norm = colors.Normalize(vmin=0,vmax=100)
scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=my_cmap)


for t in teffs:
    for l in logg:
        ax.cla()
        ax.set_xlim(0.5,3.5)
        ax.set_xlabel('Wavelength (micron)')
        ax.set_ylabel('Flux (arbitrary units)')
        ax.set_title('logg = {}'.format(l))
        ax.step(bd.specs['low']['wavelength'],bd.specs['low']['flux'],color='k')

        for add in range(0,100,10):
            p = [l,t+add]
            test_flux = mg.interp_models(p)
            ax.step(bd.specs['low']['wavelength'], test_flux,
                 color=scalar_map.to_rgba(add),label=('T={}'.format(t+add)))
        ax.legend()
        pp.savefig()
pp.close()
