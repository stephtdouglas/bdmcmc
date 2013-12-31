# Script to plot spectra for my sample
# Stephanie Douglas, 31 December 2013

import numpy as np
import matplotlib.pyplot as plt
import bdmcmc.sample.fetch as fetch


def plot_low(ax,bd,offset,label,pcolor='k'):

    wav = bd.specs['low']['wavelength']
    flu = bd.specs['low']['flux']
    unc = bd.specs['low']['unc']

    if len(wav)>1:
        norm_by = np.average(flu[abs(wav-1.27)<0.005])
        flu = flu/norm_by
        unc = flu/norm_by

        ax.step(wav,flu+offset,color=pcolor)
        ax.text(2.2,0.6+offset,label,color=pcolor)
        ax.set_xlim(0.5,3)
    else:
        print 'failed', label



#ldwarfs = fetch.fetch_12()
#bds = ldwarfs.brown_dwarfs
unums = bds.keys()

offset_jump = 1
offset_start = 3*offset_jump
offset = offset_start
plt.figure(figsize=(12,9))
ax = plt.subplot(131)
ax.set_ylim(-0.1,offset_start+1.2)
for u in unums:
    plot_low(ax,bds[u],offset,u)
    if u==unums[3]:
        ax = plt.subplot(132)
        ax.set_ylim(-0.1,offset_start+1.2)
        offset = offset_start
    elif u==unums[7]:
        ax = plt.subplot(133)
        ax.set_ylim(-0.1,offset_start+1.2)
        offset = offset_start
    else:
        offset -= offset_jump
