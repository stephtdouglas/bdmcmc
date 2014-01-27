# Script to plot spectra for my sample
# Stephanie Douglas, 31 December 2013

import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.units as u

import bdmcmc.sample.fetch as fetch


def plot_low(ax,bd,offset,label,pcolor='k'):

    wav = bd.specs['low']['wavelength']
    flu = bd.specs['low']['flux']
    unc = bd.specs['low']['unc']

    if len(wav)>1:
        #print min(wav), max(wav), flu[abs(wav-1.27)<0.005]
        norm_by = np.average(flu[abs(wav-1.27*u.um)<0.005*u.um])
        flu = flu/norm_by
        unc = flu/norm_by

        ax.step(wav,flu+offset*flu.unit,color=pcolor)
        ax.text(2.,0.9+offset,label,color=pcolor)
        ax.set_xlim(0.5,3)

        #print 'done', label, np.average(flu[abs(wav-1.27)<0.005])
    else:
        print 'failed', label


def plot_high(ax,bd,order,offset,label,pcolor='k'):

    wav = bd.specs[order]['wavelength']
    flu = bd.specs[order]['flux']
    unc = bd.specs[order]['unc']

    if len(wav)>1:
        norm_by = np.average(flu)
        #print label, order, offset, min(wav), max(wav), norm_by
        flu = flu/norm_by
        unc = flu/norm_by

        ax.step(wav,flu+offset,color=pcolor)
        ax.text(plt.gca().get_xlim()[0],1.35+offset,label,color=pcolor)
        #ax.set_xlim(0.5,3)

        #print 'done', label, np.average(flu[abs(wav-1.27)<0.005])
    else:
        print 'failed', label, order


ldwarfs = fetch.fetch_12()
bds = ldwarfs.brown_dwarfs
unums = bds.keys()
num_spts = np.zeros(12)
for i in range(12):
    num_spts[i] = float(bds[unums[i]].spt[1:])
unums = np.asarray(unums)[np.argsort(num_spts)]

offset_jump = 1
offset_start = 3*offset_jump
offset = offset_start
plt.figure(figsize=(12,9))
ax = plt.subplot(131)
ax.set_ylim(-0.1,offset_start+1.2)
ax.set_ylabel('flux (normalized at 1.27 microns)')
for unum in unums:
    plot_low(ax,bds[unum],offset,'{} {}'.format(bds[unum].spt,unum))
    print unum
    if unum==unums[3]:
        ax = plt.subplot(132)
        ax.set_ylim(-0.1,offset_start+1.2)
        offset = offset_start
        ax.set_xlabel('wavelength (micron)')
    elif unum==unums[7]:
        ax = plt.subplot(133)
        ax.set_ylim(-0.1,offset_start+1.2)
        offset = offset_start
    else:
        offset -= offset_jump
plt.savefig('ldwarfs_sample12_{}.png'.format(datetime.date.isoformat(datetime.date.today())))


pp = PdfPages('ldwarfs_12_highres.pdf')
offset_jump = 1.25
offset_start = 5*offset_jump
offset = offset_start
fig = plt.figure(figsize=(12,9))
for order in np.asarray([58,59,61,62,63,64,65]):
    offset = offset_start
    fig.suptitle('Order {}'.format(order),fontsize='x-large')
    ax = plt.subplot(121)
    ax.tick_params(labelsize='large')
    ax.set_ylim(-0.2,offset_start+1.7)
    ax.set_ylabel('flux (normalized)',fontsize='x-large')
    ax.set_xlabel('wavelength (micron)',fontsize='x-large')
    for u in unums:
        plot_high(ax,bds[u],order,offset,'{} {}'.format(bds[u].spt,u))
        #print order, u
        if u==unums[5]:
            ax = plt.subplot(122)
            ax.tick_params(labelsize='large')
            ax.set_ylim(-0.3,offset_start+1.7)
            ax.set_ylabel('flux (normalized at 1.27 microns)',
                 fontsize='x-large')
            offset = offset_start
            ax.set_xlabel('wavelength (micron)')
        else:
            offset -= offset_jump
    pp.savefig()
    plt.clf()
pp.close()
plt.close()
