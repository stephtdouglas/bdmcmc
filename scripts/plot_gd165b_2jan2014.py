# Stephanie Douglas, 2 January 2014
import datetime

import numpy as np
import matplotlib.pyplot as plt

import pyfits
import asciitable as at

from bdmcmc.spectra import *


def plot_low(ax,spectrum,offset,label,pcolor='k'):

    wav = spectrum['wavelength']
    flu = spectrum['flux']
    unc = spectrum['unc']

    if len(wav)>1:
        norm_by = np.average(flu[abs(wav-1.27)<0.005])
        flu = flu/norm_by
        unc = flu/norm_by

        ax.step(wav,flu+offset,color=pcolor,label=label)
        ax.set_xlim(0.5,3)
    else:
        print 'failed', label

plt.figure(figsize=(7,9))
ax = subplot(111)


bd = BrownDwarf('U40039')
bd.get_low()
plot_low(ax,bd.specs['low'],2,'R400, 2002-05-30','r')


bd.specs['low'] = spectrum_query(bd.sid,
         '','',filename='spex_prism_G196-3B_U40039.fits')
plot_low(ax,bd.specs['low'],1,'database, no obs-date','b')


oldfile = at.read('/home/stephanie/ldwarfs/summerAMNH/LdwarfSpectra/spex_prism_gd165b_090629.txt')
old_spectrum = {'wavelength':oldfile['col1'],'flux':oldfile['col2'],
    'unc':oldfile['col3']}
plot_low(ax,old_spectrum,0,'Burgasser Unpub 2009-06-29')
ax.set_ylim(-0.1,4)
ax.set_xlabel('Wavelength (micron)')
ax.set_ylabel('Flux (normalized at 1.27 micron)')
ax.legend(loc='best')
plt.savefig('gd165b_{}.png'.format(datetime.date.isoformat(datetime.date.today())))
