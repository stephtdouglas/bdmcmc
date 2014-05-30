import datetime
import logging

## Third-party
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import astropy.units as u
import cPickle
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.colors import Normalize
import matplotlib.cm as cm
from matplotlib.font_manager import FontProperties

import bdmcmc.bdfit, bdmcmc.spectra, bdmcmc.get_mod
import bdmcmc.make_model
import bdmcmc.mask_bands as mb
from bdmcmc.plotting.plot_random import plot_random
import bdmcmc.plotting.compare_results as cr

logging.basicConfig(level=logging.INFO)

def one_reg(ax, cropchain, model,rand_color):
    rs = plot_random(cropchain, model, ax,rand_color,False)


def plot_four(bd_name,chain_filename,plot_title,
    model_filename='/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl',
    mask_H=True):

    # Set up brown dwarf
    bd = bdmcmc.spectra.BrownDwarf(bd_name)
    bd.get_low()
    if mask_H:
        mask = mb.BandMask(bd.specs['low']['wavelength'])
        mask.mask_Hband()
        mask.make_pixel_mask()

        reverse_mask = np.delete(np.arange(len(bd.specs['low']['wavelength'])),
            mask.pixel_mask)
        mask.pixel_mask = reverse_mask

        bd.specs['low']['wavelength'] = bd.specs['low']['wavelength'][
            mask.pixel_mask]
        bd.specs['low']['flux'] = bd.specs['low']['flux'][mask.pixel_mask]
        bd.specs['low']['unc'] = bd.specs['low']['unc'][mask.pixel_mask]


    # Set up bands
    wav = bd.specs['low']['wavelength']
    full = np.where((wav>=0.9*u.um) & (wav<2.5*u.um))[0]
    Jband = np.where((wav>=0.9*u.um) & (wav<1.4*u.um))[0]
    Hband = np.where((wav>=1.4*u.um) & (wav<1.9*u.um))[0]
    Kband = np.where((wav>=1.9*u.um) & (wav<2.5*u.um))[0]

    bands = {'J':Jband,'H':Hband,'K':Kband,'full':full}
    band_names = ['full','J','H','K']

    # get model
    am = bdmcmc.get_mod.AtmoModel(model_filename)


    # Set up figure
    fig = plt.figure(figsize=(10,8))
    gs = gridspec.GridSpec(2,3)
    ax1 = plt.subplot(gs[0,:])
    #axb = [plt.subplot(gs[1,i]) for i in range(3)]
    #axes = np.append(ax1,axb)
    ax2 = plt.subplot(gs[1,:])
    axes = np.append(ax1,ax2)
    logging.debug(axes)

    # Match colors to compare_results
    cmap = cm.get_cmap("paired")
    color_norm = Normalize(vmin=0,vmax=4)
    scalar_map = cm.ScalarMappable(norm=color_norm,cmap=cmap)
    plot_colors = [scalar_map.to_rgba(i) for i in range(4)]
    rand_color = [plot_colors[i] for i in [3,2,0,1]] # reset from ['H','K','J','full'][i]


    for i in range(4):
        b = band_names[i]
        logging.debug(b)
        band_spectrum = {
            'wavelength':bd.specs['low']['wavelength'][bands[b]],
            'flux':bd.specs['low']['flux'][bands[b]],
            'unc':bd.specs['low']['unc'][bands[b]]}

        mg = bdmcmc.make_model.ModelGrid(band_spectrum,am.model,
            am.params,smooth=False)

        cpfile = open(chain_filename.replace('full',b),'rb')
        chains = cPickle.load(cpfile)
        cpfile.close()

        cropchain = chains.reshape((-1,np.shape(chains)[-1]))
        
        #one_reg(axes[i],cropchain,mg)

        #if i>1:
        #    axes[i].set_ylabel('')
        #    axes[i].tick_params(labelleft=False,labelright=False)

        if i>0:
            ax = axes[1]
        else:
            ax = axes[0]
        one_reg(ax,cropchain,mg,rand_color[i])
        
        if i>0:
              ax.set_xlim(axes[0].get_xlim())
              yl = ax.get_ylim()
              ax.plot((1.4,1.4),yl,'k-',lw=2)
              ax.plot((1.9,1.9),yl,'k-',lw=2)


chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/1228-1547 Marley H 2014-05-13_chains.pkl'
results_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/1228-1547_2014-05-13_all.pkl'

chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/0036+1821 Marley H 2014-05-13_chains.pkl'
results_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/0036+1821_2014-05-13_all.pkl'


model_file = '/home/stephanie/ldwarfs/modelSpectra/SpeX_marley_nolowg.pkl'

obj_name = chain_file.split('/')[-1].split()[0]
date = chain_file.split('/')[-2].split('_')[-1]

#plot_four(obj_name,chain_file,'',model_file)
#plt.show()
#plt.savefig('region_ex_{}_{}.png'.format(obj_name,date),dpi=600,bbox_inches='tight')

infile = open(results_file,'rb')
results = cPickle.load(infile)
infile.close()
param_labels = ['log(g)','Fsed','Teff','ln(s)']
fig, axes = cr.corner(results[0],results[1],0.0,param_labels,['H','K','J','full'])

for i in range(3):
    axes[i,i].set_visible(False)
    axes[i,i].set_frame_on(False)
for i in range(4):
    if i<3:
        axes[2,i].set_xticklabels(axes[3,i].get_xticks())
        axes[2,i].set_xlabel(param_labels[i])
        axes[2,i].set_ylim(1500,2050)
        [l.set_rotation(45) for l in axes[2,i].get_xticklabels()]
    axes[3,i].set_visible(False)
    axes[3,i].set_frame_on(False)

axes[1,0].set_yticks(axes[2,1].get_xticks())
fontP = FontProperties()
fontP.set_size('x-small')
axes[1,0].legend(loc='best',numpoints=1,prop = fontP)

plt.show()
plt.savefig('region_res_{}_{}.png'.format(obj_name,date),dpi=600,bbox_inches='tight')
