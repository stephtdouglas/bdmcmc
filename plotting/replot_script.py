
import logging

import numpy as np
from bdmcmc.plotting.replot import replot
import asciitable as at
import cPickle
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages

result_file = '/home/stephanie/ldwarfs/batch_ldwarfs/Marley_2014_04_29/results_2014-04-29.dat'
model_file = '/home/stephanie/ldwarfs/modelSpectra/SpeX_marley.pkl'

dat = at.read(result_file)
spts,filename = dat['spt'],dat['filename']
dlen = len(dat)

bands = ['J','H','K','full']
mask = {'J':False,'H':True,'K':False,'full':True}

for i in range(dlen):
    fsplit = filename[i].split('_')
    obj_name = fsplit[0]
    date = fsplit[1]
    pp = PdfPages('{}_{}_replot.pdf'.format(obj_name,date))
    for b in bands: 
        chain_file = "{} Marley {} {}_chains.pkl".format(obj_name,b,date)
        replot(obj_name,chain_file,'{}_{}_replot'.format(obj_name,date),
            model_file,mask_H = mask[b])
        pp.savefig()
        plt.close()
    pp.close()

