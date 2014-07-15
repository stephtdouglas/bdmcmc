import logging

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.units as u

import bdmcmc.spectra, bdmcmc.batch
from bdmcmc.config import *

def fetch_sxd():
    """
    Searches the database for all SpeX Cross-Dispersed spectra
    and returns the object names and 

    """

    # Select all spex spectra from the database
    all_spex = db.dict.execute("SELECT s.obs_date, s.filename, s.wavelength,"
        " so.shortname, s.source_id FROM spectra AS s JOIN sources AS so ON "
        "so.id=s.source_id WHERE s.telescope_id=7").fetchall()

    # only keep object names and observation_dates for those with resolution 
    # over 1000 (Prism has R~200-400, XD has R~2000)
    spex_sxd = [[i['shortname'],i['obs_date']] for i in all_spex if 
        (np.average(i['wavelength'][:-1] / np.diff(i['wavelength'])))>1000]

    return spex_sxd

def make_sxd_batch(model_name="marley",model_file="SXD_marley.pkl"):

    spex_sxd = fetch_sxd()

    h = open('submit_all_{}.sh'.format(model_name),'w')
    for name,date in spex_sxd:
        h.write('qsub run_{}_{}.sh\n'.format(model_name,name))

        f = open('run_{}_{}.py'.format(model_name,name),'w')
        f.write('import logging\n')
        f.write('from datetime import date\n\n')
        f.write('from bdmcmc.batch import OneBatch\n\n')
        f.write('logging.basicConfig(level=logging.INFO)\n\n')
        f.write("ob = OneBatch('{}',"
             "'/vega/astro/users/sd2706/modelSpectra/{}',"
             "obs_date='{}')\n".format(name,model_file,date))
        f.close()

        g = open('run_{}_{}.sh'.format(model_name,name),'w')
        g.write("#!/bin/sh\n\n")
        g.write("# Directives\n")
        g.write("#PBS -N {}{}\n".format(model_name,name))
        g.write("#PBS -W group_list=yetiastro\n")
        g.write("#PBS -l nodes=1,walltime=20:00:00,mem=4000mb\n")
        g.write("#PBS -M sd2706@columbia.edu \n")
        g.write("#PBS -m a\n")
        g.write("#PBS -V\n\n")
        g.write("# Set output and error directories\n")
        g.write("#PBS -o localhost:/vega/astro/users/sd2706/outputs/ \n")
        g.write("#PBS -e localhost:/vega/astro/users/sd2706/outputs/ \n")
        g.write("/vega/astro/users/amp2217/yt-x86_64/bin/python"
             " run_{}_{}.py\n".format(model_name,name))
        g.write("# End of script\n")
        g.close()
    h.close()

def plot_all_sxd():

    spex_sxd = fetch_sxd()

    pp = PdfPages('All_SXD.pdf')
    plt.figure()
    ax = plt.subplot(111)
    for name,date in spex_sxd:
         if (name!=None) and (name!='0559-1404') and (name!='1254-0122'):
             bd = bdmcmc.spectra.BrownDwarf(name)
             bd.get_low(obs_date=date)
             print name,date,bd.specs['low']['slit_width']
             ax.step(bd.specs['low']['wavelength'],bd.specs['low']['flux'])
             ax.set_title(name)
             ax.set_xlim(0.8,2.6)
             pp.savefig()
             ax.cla()
    pp.close()


def plot_sxd_res():
    spex_sxd = fetch_sxd()

    fig1 = plt.figure()
    ax1 = plt.subplot(111)
    for name,date in spex_sxd:
         if (name!=None) and (name!='0559-1404') and (name!='1254-0122'):
             bd = bdmcmc.spectra.BrownDwarf(name)
             bd.get_low(obs_date=date)
             print name,date,bd.specs['low']['slit_width']
             w = bd.specs['low']['wavelength']
             wdiff = w[2:]-w[:-2]
             R = w[:-2] / wdiff * (
                 0.3/bd.specs['low']['slit_width'].value)
             Rclean = R[R>1000]
             wclean = w[:-1][R>1000]#w[:-2][wdiff>0.01*u.um]
             ax1.plot(wclean, Rclean,'.',label=name)
#    ax1.legend(ncol=3)

def calc_sxd_res():
    spex_sxd = fetch_sxd()
    all_w = np.array([])
    all_R = np.array([])

    for name,date in spex_sxd:
         if (name!=None) and (name!='0559-1404') and (name!='1254-0122'):
             bd = bdmcmc.spectra.BrownDwarf(name)
             bd.get_low(obs_date=date)
             if bd.specs['low']['slit_width']==(0.5*u.um):
                  w = bd.specs['low']['wavelength']
                  wdiff = w[2:]-w[:-2]
                  R = w[:-2] / wdiff * (0.3/0.5)
                  Rclean = R[(R>1000) & (R<1375)]
                  wclean = w[:-1][(R>1000) & (R<1375)]
                  np.append(all_w,wclean)
                  np.append(all_R,Rclean)
         all_w.sort()
         all_R.sort()
         split_points = [0.93891,1.1277,1.41055,1.85]

         # there are 5 orders
         # for each order, fit a line to the R vs. lambda curve
         # Then use that point to calculate the fitted R as a
         # function of all the data wavelength points
         # Smooth the models to the right R, but resample onto this
         # horrendously oversampled wavelength grid, so that the code
         # can just interpolate

plot_sxd_res()

#plot_all_sxd()



#make_sxd_batch()
