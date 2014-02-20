# Script to check SpeX wavelength grids
# copied from plot_12_hilow.py
# Stephanie Douglas, 20 February 2014

import datetime

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.units as u

import bdmcmc.sample.fetch as fetch


def calc_res(bd,label):

    wav = bd.specs['low']['wavelength']
    #flu = bd.specs['low']['flux']
    #unc = bd.specs['low']['unc']
    data_wave=wav

    # Calculate resolution array
    res = np.zeros(len(wav))
    res[2:] = data_wave[2:]/(data_wave[2:]-data_wave[:-2])
    # need to fill in first 2 elements so keep same array length
    # the end elements are very noisy anyway so it shouldn't be an issue
    res[:2] = res[2:4] 
    logging.debug(str(res[:10]))

#    plt.plot(np.arange(len(res)),res,label=label)
    plt.plot(wav,res,label=label)


ldwarfs = fetch.fetch_12()
bds = ldwarfs.brown_dwarfs
unums = bds.keys()
num_spts = np.zeros(12)
for i in range(12):
    num_spts[i] = float(bds[unums[i]].spt[1:])
unums = np.asarray(unums)[np.argsort(num_spts)]

plt.figure(figsize=(10,8))

for unum in unums:
    calc_res(bds[unum],'{} {}'.format(bds[unum].spt,unum))
    
plt.legend(loc=4)
plt.xlabel('wavelength (micron)',fontsize='x-large')
plt.ylabel('R',fontsize='x-large')

plt.savefig('check_spex_lambda.png')
