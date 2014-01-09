# Plot all versions of all the models I have
# The ORIGINAL versions, not the cropped ones from the pkl files
# Stephanie Douglas, 9 January 2014
################################################################################

import os
from datetime import date

import asciitable as at
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.colors as colors
import matplotlib.cm as cmx
from scipy.io.idl import readsav
from matplotlib.backends.backend_pdf import PdfPages

def reverse(f):
    l = len(f)
    flux = zeros(l)
    for i in arange(l):
        #print i,l-i
        flux[i]=f[l-1-i]
    return flux

clr = ['maroon','red','orangered','orange','gold','yellow','lime','green','aqua','blue','fuchsia','purple']

"""

### Marley low-res

mlowpath = '/home/stephanie/ldwarfs/summerAMNH/modelSpectra/marley_lowres/'
#Read all files by cycling through temp, grav, and fsed

temps = numpy.arange(1200,2400,100)
gravs = [100,300,1000,3000]
fseds = [1,2,3]

notthere = []
mod = []
flu = []

models = {'teff':np.array([]),'logg':np.array([]),'fsed':np.array([]),
     'wsyn':np.array([]),'fsyn':[]}

print 'starting now'
for t in temps:
    for g in gravs:
        for f in fseds:
            #print t,g,f
            filename = mlowpath+'sp_t'+str(t)+'g'+str(g)+'f'+str(f)
            if os.path.exists(filename):
                mod = at.read(filename,data_start=4)
                wav = sort(mod['col1'])
                flu1 = reverse(mod['col2'])
                flu = flu1*(3e7)/wav**2
                models['teff'] = np.append(models['teff'],t)
                models['logg'] = np.append(models['logg'],g)
                models['fsed'] = np.append(models['fsed'],f)
                models['fsyn'].append(np.asarray(flu))
                if t == 1200:
                    models['wsyn'] = np.asarray(wav)/10000.0
            else:
                models['teff'] = np.append(models['teff'],t)
                models['logg'] = np.append(models['logg'],g)
                models['fsed'] = np.append(models['fsed'],f)
                models['fsyn'].append(numpy.ones(28485))
                notthere.append((t,g,f))
pp = PdfPages('marleylowplots_{}.pdf'.format(date.isoformat(date.today())))

plt.figure(figsize=(10,7.5))
ax = plt.subplot(111)

norm_reg = np.where(abs(np.asarray(models['wsyn'])-1.27)<0.05)[0]


params = ['teff','logg','fsed']
for p in range(3):
    to_vary = params[p]
    to_cycle = np.delete(params,p)
    i_cycle = np.delete(np.arange(3),p)
    i0 = i_cycle[0]
    i1 = i_cycle[1]
    vals_0 = np.unique(models[params[i0]])
    vals_1 = np.unique(models[params[i1]])
    print to_vary,params[i0], params[i1]

    clr_skip = 1
    if len(np.unique(models[to_vary]))<len(clr):
        clr_skip = len(clr)/len(np.unique(models[to_vary]))


    for j in vals_0:
        for k in vals_1:
            static_vals = np.where((models[params[i0]]==j) & 
                (models[params[i1]]==k))[0]
            clr_count=0
            for i in static_vals:
                if len(models['fsyn'][i])==1:
                    models['fsyn'][i] = models['fsyn'][i][0]
                if sum(models['fsyn'][i])!=len(models['fsyn'][i]):
                    norm_by = np.average(models['fsyn'][i][norm_reg])
                    ax.plot(models['wsyn'],models['fsyn'][i]/norm_by,
                            color=clr[clr_count],
                            label=str(int(models[to_vary][i])))
                clr_count += clr_skip
            ax.legend(title=to_vary)
            ax.set_xlabel('Wavelength (micron)')
            ax.set_ylabel('Flux (normalized at 1.27 micron)')
            ax.set_title('{} = {}, {} = {}'.format(params[i0],j,params[i1],k))
            pp.savefig()
            ax.cla()
pp.close()

"""


### Dusty low-res
modelpath = '/home/stephanie/ldwarfs/summerAMNH/modelSpectra/'
dustylowfile = modelpath+'modelspeclowresldwarfs.save'
dl = readsav(dustylowfile)
d = len(dl.modelspec.teff)
models = {'teff':np.array([]),'logg':np.array([]),'wsyn':np.array([]),
     'fsyn':[]}

print 'starting now'
for i in numpy.arange(d):
    models['teff'] = np.append(models['teff'], dl.modelspec.teff[i])
    models['logg'] = np.append(models['logg'], dl.modelspec.logg[i])
    models['fsyn'].append(dl.modelspec.fsyn[i])
models['wsyn'] = np.asarray(dl.wsyn)/10000.
print 'models retrieved'


pp = PdfPages('dustyhighres_{}.pdf'.format(date.isoformat(date.today())))
plt.figure(figsize=(10,7.5))
ax = plt.subplot(111)

norm_reg = np.where(abs(np.asarray(models['wsyn'])-1.27)<0.05)[0]
params = ['teff','logg']
for p in range(2):
    to_vary = params[p]
    to_cycle = np.delete(params,p)
    i_cycle = np.delete(np.arange(3),p)
    i0 = i_cycle[0]
    vals_0 = np.unique(models[params[i0]])
    print to_vary,params[i0]

    my_cmap = plt.get_cmap('gist_rainbow')
    c_norm = colors.Normalize(vmin=0,vmax=len(np.unique(models[to_vary])))
    c_values = np.arange(0,len(np.unique(models[to_vary])))
    scalar_map = cmx.ScalarMappable(norm=c_norm, cmap=my_cmap)

    clr_count=0
    for j in vals_0:
        static_vals = np.where((models[params[i0]]==j))[0]
        clr_count=0
        for i in static_vals:
            if len(models['fsyn'][i])==1:
                models['fsyn'][i] = models['fsyn'][i][0]
            if sum(models['fsyn'][i])!=len(models['fsyn'][i]):
                norm_by = np.average(models['fsyn'][i][norm_reg])
                ax.plot(models['wsyn'],models['fsyn'][i]/norm_by,
                        color=scalar_map.to_rgba(c_values[clr_count]),
                        label=str(round(models[to_vary][i],p)))
            clr_count += 1
        ax.legend(title=to_vary,prop={'size':'small'})
        ax.set_xlim(0.5,3.5)
        ax.set_xlabel('Wavelength (micron)')
        ax.set_ylabel('Flux (normalized at 1.27 micron)')
        ax.set_title('{} = {}'.format(params[i0],round(j,i0)))
        pp.savefig()
        ax.cla()
pp.close()
