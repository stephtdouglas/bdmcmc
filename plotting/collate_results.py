
import logging

import matplotlib.pyplot as plt
import numpy as np
import bdmcmc.plotting.compare_results as cr
import asciitable as at
import cPickle

result_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/results_2014-05-13.dat'

dat = at.read(result_file)
spts,filename = dat['spt'],dat['filename']
dlen = len(dat)

meds = []
errs = []

for i in range(dlen):
    infile = open(filename[i],'rb')
    med_temp, err_temp0 = cPickle.load(infile)
    infile.close()
    # swap upper and lower uncertainties so that they get plotted correctly
    err_temp = [[err_temp0[i][j][::-1] for j in range(np.shape(err_temp0)[1])]
                for i in range(np.shape(err_temp0)[0])]
    meds.append(med_temp)
    errs.append(err_temp)

all_med = [[[meds[k][i][j] for k in range(dlen)] for j in range(4)] for i in range(4)]
all_err = [[np.array([errs[k][i][j] for k in range(dlen)]).reshape((2,-1)) for j in range(4)] for i in range(4)]

cr.corner(all_med,all_err,spts,['logg','fsed','teff','ln(s)'],['H','K','J','full'])

plt.savefig('compare10_Marley.pdf')
