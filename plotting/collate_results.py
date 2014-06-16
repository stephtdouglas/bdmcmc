
import logging

import matplotlib.pyplot as plt
import numpy as np
import bdmcmc.plotting.compare_results as cr
import asciitable as at
import cPickle

#result_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-13/results_2014-05-13.dat'
result_file = '/home/stephanie/ldwarfs/batch_ldwarfs/SpeX_2014-04-16/results_2014-04-16.dat'

dat = at.read(result_file)
spts,filename = dat['spt'],dat['filename']
dlen = len(dat)

meds = []
errs = []

for i in range(dlen):
    infile = open(filename[i],'rb')
    med_temp, err_temp0 = cPickle.load(infile)
    # med_temp has the shape (num_params, num_runs)
    # err_temp0 has the shape (num_params, num_runs, 2)
    infile.close()
    # swap upper and lower uncertainties so that they get plotted correctly
    err_temp = [[err_temp0[i][j][::-1] for j in range(np.shape(err_temp0)[1])]
                for i in range(np.shape(err_temp0)[0])]
    meds.append(med_temp)
    errs.append(err_temp)

# k gives the object count
# j gives the run count
# i gives the parameter count
# np.shape(meds) corresponds to (k, i, j)
num_params = np.shape(meds)[1]
param_counter = range(num_params) #i
num_runs = np.shape(meds)[2]
run_counter = range(num_runs) #j


all_med = [[[meds[k][i][j] for k in range(dlen)] for j in run_counter] 
           for i in param_counter]
all_err = [[np.array([errs[k][i][j] for k in range(dlen)]).reshape((2,-1)) 
            for j in run_counter] for i in param_counter]

#cr.corner(all_med,all_err,spts,['logg','fsed','teff','ln(s)'],['H','K','J','full'])

#plt.savefig('compare10_Marley.pdf')
