
import logging

import numpy as np
import bdmcmc.plotting.compare_results as cr
import asciitable as at

result_file = '/home/stephanie/ldwarfs/batch_ldwarfs/results_2014-04-16.dat'

dat = at.read(result_file)
spts,filename = dat['spt'],dat['filename']
dlen = len(dat)

meds = []
errs = []

for i in range(dlen):
    infile = open(filename[i],'rb')
    med_temp, err_temp = cPickle.load(infile)
    infile.close()
    meds.append(med_temp)
    errs.append(err_temp)

all_med = [[[meds[k][i][j] for k in range(dlen)] for j in range(4)] for i in range(3)]
all_err = [[np.array([errs[k][i][j] for k in range(dlen)]).reshape((2,-1)) for j in range(4)] for i in range(3)]

cr.corner(all_med,all_err,spts,['logg','teff','ln(s)'],['J','H','K','full'])

plt.savefig('compare10.pdf')
