
import matplotlib
print matplotlib.get_backend()
print matplotlib.matplotlib_fname()
matplotlib.use('Agg')

import logging

import numpy as np
from bdmcmc.plotting.spt_plot import spt_plot
import asciitable as at
import cPickle

logging.basicConfig(level=logging.INFO)

result_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-07b/results_2014-05-07b.dat'

dat = at.read(result_file)
spts,filename = dat['spt'],dat['filename']
dlen = len(dat)
logging.debug('dlen {}'.format(dlen))

nparams = 4
nruns = 4

meds = []
errs = []

for i in range(dlen):
    infile = open("/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014-05-07b/"+filename[i],'rb')
    med_temp, err_temp = cPickle.load(infile)
    infile.close()
    meds.append(med_temp)
    errs.append(err_temp)

all_med = [[[meds[k][i][j] for k in range(dlen)] for j in range(nruns)] for i in range(nparams)]
all_err = [[np.array([errs[k][i][j] for k in range(dlen)]).reshape((2,-1)) for j in range(nruns)] for i in range(nparams)]

spt_plot(all_med,all_err,spts,['logg','fsed','teff','ln(s)'],['J','H','K','full'],"Marley",single_figure=False,extents=[[3.,6.0],[0,5],[1200,2500],[-36,-30]],model_extents=[[4.5,5.5],[0,5],[1500,2400],[-36,-30]])

######################################################
result_file = '/home/stephanie/ldwarfs/batch_ldwarfs/SpeX_2014-04-16/results_2014-04-16.dat'

dat = at.read(result_file)
spts,filename = dat['spt'],dat['filename']
dlen = len(dat)
logging.debug('dlen {}'.format(dlen))

nparams = 3
nruns = 4

meds = []
errs = []

for i in range(dlen):
    infile = open('/home/stephanie/ldwarfs/batch_ldwarfs/SpeX_2014-04-16/'+filename[i],'rb')
    med_temp, err_temp = cPickle.load(infile)
    infile.close()
    meds.append(med_temp)
    errs.append(err_temp)

all_med = [[[meds[k][i][j] for k in range(dlen)] for j in range(nruns)] for i in range(nparams)]
all_err = [[np.array([errs[k][i][j] for k in range(dlen)]).reshape((2,-1)) for j in range(nruns)] for i in range(nparams)]

spt_plot(all_med,all_err,spts,['logg','teff','ln(s)'],['J','H','K','full'],"Barman-Dusty",single_figure=False,extents=[[3.,6.0],[1200,2500],[-36,-30]],model_extents=[[3,5.5],[1400,2400],[-36,-30]])
