import datetime
import logging

## Third-party
#import matplotlib
#matplotlib.use('agg')
import matplotlib as mpl
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

mpl.rcParams['axes.linewidth'] = 1.5

def quantile(x,quantiles):
    # From DFM's triangle code
    xsorted = sorted(x)
    qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
    return zip(quantiles,qvalues)

def one_reg(ax, cropchain, model,rand_color):
#    plt.figure()
#    ax2 = plt.subplot(111)
#    rs = plot_random(cropchain, model, ax2,rand_color,False)
    rs = plot_random(cropchain, model, ax,rand_color,False)


def plot_four(bd_name,chain_filename,plot_title,
    model_filename='/home/stephanie/ldwarfs/modelSpectra/SpeX_dusty.pkl',
    mask_H=True,obs_date=None):

    medians = np.zeros(8).reshape((2,4))
    errors = np.zeros(16).reshape((2,4,2))

    # Set up brown dwarf
    bd = bdmcmc.spectra.BrownDwarf(bd_name)
    bd.get_med("sxd",obs_date=obs_date)
    good = np.where((np.isnan(bd.specs["med"]["wavelength"])==False) & 
        (np.isnan(bd.specs["med"]["flux"])==False) & 
        (np.isnan(bd.specs["med"]["unc"])==False))
    if len(good)<len(bd.specs["med"]["wavelength"]):
        bd.specs['med']['wavelength'] = bd.specs['med']['wavelength'][good]
        bd.specs['med']['flux'] = bd.specs['med']['flux'][good]
        bd.specs['med']['unc'] = bd.specs['med']['unc'][good]


    print bd.specs
    if mask_H:
        mask = mb.BandMask(bd.specs['med']['wavelength'])
        mask.mask_Hband()
        mask.make_pixel_mask()

        reverse_mask = np.delete(np.arange(len(bd.specs['med']['wavelength'])),
            mask.pixel_mask)
        mask.pixel_mask = reverse_mask

        print reverse_mask
        bd.specs['med']['wavelength'] = bd.specs['med']['wavelength'][
            mask.pixel_mask]
        bd.specs['med']['flux'] = bd.specs['med']['flux'][mask.pixel_mask]
        bd.specs['med']['unc'] = bd.specs['med']['unc'][mask.pixel_mask]


    # Set up bands
    wav = bd.specs['med']['wavelength']
    full = np.where((wav>=0.9*u.um) & (wav<2.5*u.um))[0]
    Jband = np.where((wav>=0.9*u.um) & (wav<1.4*u.um))[0]
    Hband = np.where((wav>=1.4*u.um) & (wav<1.9*u.um))[0]
    Kband = np.where((wav>=1.9*u.um) & (wav<2.5*u.um))[0]

    bands = {'J':Jband,'H':Hband,'K':Kband,'full':full}
    band_names = ['full','J','H','K']

    # get model
    am = bdmcmc.get_mod.AtmoModel(model_filename)
    print "got model"

    # Set up figure
    fig = plt.figure(figsize=(10,12))
    gs = gridspec.GridSpec(4,4) # nrow, ncol
    ax1 = plt.subplot(gs[0,:])
    ax2 = plt.subplot(gs[1,:2])
    ax3 = plt.subplot(gs[2,:2])
    ax4 = plt.subplot(gs[3,:2])
#    axes = np.append(np.append(ax1,ax2),np.append(ax3,ax4))
    axes = np.append(np.append(ax1,ax2),np.append(ax3,ax4))
    logging.debug(axes)
    print "set up axes"

    # Match colors to compare_results
    cmap = cm.get_cmap("paired")
    color_norm = Normalize(vmin=0,vmax=4)
    scalar_map = cm.ScalarMappable(norm=color_norm,cmap=cmap)
    plot_colors = [scalar_map.to_rgba(i) for i in range(4)]
#    rand_color = [plot_colors[i] for i in [3,2,0,1]] # reset from ['H','K','J','full'][i]
    rand_color = [plot_colors[i] for i in [3,0,2,1]] # reset from ['H','J','K','full'][i]
    rand_color[2] = [99.0/255,202.0/255,0,1]

#    texty = [max(bd.specs['med']['flux'].value)*0.9,
#             max(bd.specs['med']['flux'].value)*0.9,
#             bd.specs['med']['flux'][-1].value*1.1,
#             bd.specs['med']['flux'][-1].value]
    texty = [max(bd.specs['med']['flux'].value)*0.9,
             max(bd.specs['med']['flux'].value)*0.9,
             bd.specs['med']['flux'][-1].value*1.1,
             bd.specs['med']['flux'][-1].value]
#    textx = [0.95,1.01,1.51,2.01]
    textx = [1.1,1.17,1.51,2.01]
    print texty,textx

    for i in range(4):
        b = band_names[i]
        logging.debug(b)
        band_spectrum = {
            'wavelength':bd.specs['med']['wavelength'][bands[b]],
            'flux':bd.specs['med']['flux'][bands[b]],
            'unc':bd.specs['med']['unc'][bands[b]]}

        print band_spectrum["wavelength"][0:10]
        print band_spectrum["flux"][0:10]
        print band_spectrum["unc"][0:10]

        mg = bdmcmc.make_model.ModelGrid(band_spectrum,am.model,
            am.params,smooth=False)

        cpfile = open(chain_filename.replace('full',b),'rb')
        chains = cPickle.load(cpfile)
        cpfile.close()

        cropchain = chains.reshape((-1,np.shape(chains)[-1]))
        for j in range(2):
            quantiles = quantile(cropchain[:,j],[.05,.5,.95])
            medians[j][i] = quantiles[1][1]
            # -row1, +row2
            print j,quantiles
            uperr = quantiles[2][1] - quantiles[1][1]
            dnerr = quantiles[1][1] - quantiles[0][1]
            errors[j][i] = [dnerr,uperr]
        
        #one_reg(axes[i],cropchain,mg)

        #if i>1:
        #    axes[i].set_ylabel('')
        #    axes[i].tick_params(labelleft=False,labelright=False)

        ax = axes[i]
        one_reg(ax,cropchain,mg,rand_color[i])
        ax_ylim = ax.get_ylim()
        texty[i] = ax_ylim[0]+(ax_ylim[1]-ax_ylim[0])/10.

        ax.tick_params(labelsize="large")
        ax.tick_params(width=1.5,which="both")
        ax.set_yticklabels([])
        
    yl = axes[0].get_ylim()
    #ax1.set_ylim(yl)
    #ax1.set_xlim(0.92,2.45)
    axes[0].set_xlim(1.05,2.45)
    axes[0].plot((1.4,1.4),yl,':',lw=2,color="Grey")
    axes[0].plot((1.9,1.9),yl,':',lw=2,color="Grey")

    axes[0].text(textx[0],texty[0],"Data",color="k",fontsize='large')
    axes[1].text(textx[1],texty[1],"J Data",color="k",fontsize='large')
    axes[2].text(textx[2],texty[2],"H Data",color="k",fontsize='large')
    axes[3].text(textx[3],texty[3],"K Data",color="k",fontsize='large')

    axes[0].text(textx[0]+0.1,texty[0],"Models",color=rand_color[0],fontsize='large')
    axes[1].text(textx[1]+0.05,texty[1],"Models",color=rand_color[1],fontsize='large')
    axes[2].text(textx[2]+0.08,texty[2],"Models",color=rand_color[2],fontsize='large')
    axes[3].text(textx[3]+0.1,texty[3],"Models",color=rand_color[3],fontsize='large')


    plt.subplots_adjust(left=0.05,right=0.98,top=0.99,bottom=0.05,hspace=0.25)

    return rand_color,medians,errors


chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/1228-1547 Marley H 2014-05-13_chains.pkl'
results_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/1228-1547_2014-05-13_all.pkl'

chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/0036+1821 Marley H 2014-05-13_chains.pkl'
results_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/0036+1821_2014-05-13_all.pkl'

chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl_2014-08-25_med/2057-0252_full_BTSettl_2014-08-25_chains.pkl'
results_file = '/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl_2014-08-25_med/2057-0252_2014-08-25_all.pkl'

chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl10_2014-11-05_sxd/2057-0252_full_BTSettl_2014-11-05_chains.pkl'
results_file = '/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl10_2014-11-05_sxd/2057-0252_2014-11-05_all.pkl'

chain_file = '/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl10_2014-11-05_sxd/2208+2921_full_BTSettl_2014-11-05_chains.pkl'
results_file = '/home/stephanie/ldwarfs/batch_ldwarfs/BTSettl10_2014-11-05_sxd/2208+2921_2014-11-05_all.pkl'


model_file = '/home/stephanie/ldwarfs/modelSpectra/SXD_r2000_BTS13.pkl'

obj_name = chain_file.split('/')[-1].split("_")[0]
date = chain_file.split('/')[-2].split('_')[-1]

rand_color,medians,errors = plot_four(obj_name,chain_file,'',model_file,obs_date="2006-08-20")
plt.show()
plt.savefig('region_ex_{}_{}.png'.format(obj_name,date),dpi=600,bbox_inches='tight')

"""
#infile = open(results_file,'rb')
#results = cPickle.load(infile)
#infile.close()
#param_labels = ['log(g)','Fsed','Teff','ln(s)']
param_labels = ['log(g)','Teff']
K = 2
#medians = results[0][:K]
#err_temp0 = results[1][:K]
extents = [[2,6],[1200,3000]]
# swap upper and lower uncertainties so that they get plotted correctly
#errors = [[err_temp0[i][j][::-1] for j in range(np.shape(err_temp0)[1])]
#           for i in range(np.shape(err_temp0)[0])]
errors[0] = np.sqrt(errors[0]**2 + 0.25**2)
errors[1] = np.sqrt(errors[1]**2 + 50**2)
fig, axes = cr.corner(medians,errors,0.0,param_labels,['H','K','J','full'],run_colors=rand_color,extents=extents)


for i in range(K):
    axes[i,i].set_visible(False)
    axes[i,i].set_frame_on(False)
#for i in range(4):
#    if i<3:
#        axes[2,i].set_xticklabels(axes[3,i].get_xticks())
#        axes[2,i].set_xlabel(param_labels[i])
#        axes[2,i].set_ylim(1500,2050)
#        [l.set_rotation(45) for l in axes[2,i].get_xticklabels()]
#    axes[3,i].set_visible(False)
#    axes[3,i].set_frame_on(False)

#axes[1,0].set_yticks(axes[2,1].get_xticks())
fontP = FontProperties()
fontP.set_size('x-small')
axes[0,0].legend(loc='best',numpoints=1,prop = fontP)
axes[0,0].set_xlim(extents[0])
axes[0,0].set_ylim(extents[1])
axes[0,0].set_xlabel("log(g)",fontsize="xx-large")
axes[0,0].set_ylabel(r"T$eff$",fontsize="xx-large")
axes[0,0].tick_params(labelsize="large")

plt.show()
plt.savefig('region_res_{}_{}.png'.format(obj_name,date),dpi=600,bbox_inches='tight')
"""
