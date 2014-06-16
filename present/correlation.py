import datetime
import logging

## Third-party
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np
import cPickle


from bdmcmc.plotting.triangle import hist2d

basepath = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/'

def plot_corr(chain_files,param_names,x_name,y_name,extents):

    num_plots = len(chain_files)
    x_loc = np.where([param_names[i]==x_name for i in range(4)])[0][0]
    y_loc = np.where([param_names[i]==y_name for i in range(4)])[0][0]
    print num_plots, x_loc, y_loc

    K = num_plots/2
    fig = plt.figure(figsize=(8,8))
    setup_axes = [plt.subplot(K,K,i+1) for i in np.arange(num_plots)]
    axes = np.array(setup_axes)

    for i in range(num_plots):
        ax = axes[i]

        obj_name = chain_files[i].split()[0]
        band = chain_files[i].split()[2]

        infile = open(basepath+chain_files[i])
        chain = cPickle.load(infile)
        infile.close()
        cropchain = chain.reshape((-1,4))

        hist2d(cropchain[:,y_loc],cropchain[:,x_loc],ax=ax,
            extent=[extents[y_loc],extents[x_loc]])
        ax.set_title('{} {}'.format(obj_name,band))

        if i < K:
            ax.set_xticklabels([])
        else:
            [l.set_rotation(45) for l in ax.get_xticklabels()]
            ax.set_xlabel(y_name)

        if i%2==1:            
            ax.set_yticklabels([])
        else:
            [l.set_rotation(45) for l in ax.get_yticklabels()]
            ax.set_ylabel(x_name)


params = ['log(g)','Fsed','Teff','ln(s)']
extents = [[4.5,5.5],[1,4.5],[1500,2200],[0,1]]


teff_files = ['1228-1547 Marley K 2014-05-13_chains.pkl',
    '1507-1627 Marley K 2014-05-13_chains.pkl',
    '0015+3516 Marley J 2014-05-13_chains.pkl',
    '0205-1159 Marley K 2014-05-13_chains.pkl']

plot_corr(teff_files,params,'Fsed','Teff',extents)
plt.savefig('correlation_tf.png',dpi=600,bbox_inches='tight')

logg_files = ['1228-1547 Marley K 2014-05-13_chains.pkl',
    '1228-1547 Marley J 2014-05-13_chains.pkl',
#    '1507-1627 Marley full 2014-05-13_chains.pkl',
    '1506+1321 Marley K 2014-05-13_chains.pkl',
    '1506+1321 Marley full 2014-05-13_chains.pkl',
#    '0015+3516 Marley J 2014-05-13_chains.pkl',
#    '0205-1159 Marley full 2014-05-13_chains.pkl'
    ]

plot_corr(logg_files,params,'Fsed','log(g)',extents)
plt.savefig('correlation_gf.png',dpi=600,bbox_inches='tight')
