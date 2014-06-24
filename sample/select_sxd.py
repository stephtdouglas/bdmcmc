import logging

import numpy as np

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

    # only object names and observation_dateskeep those with resolution 
    # over 1000 (Prism has R~200-400, XD has R~2000)
    spex_sxd = [[i['shortname'],i['obs_date']] for i in all_spex if 
        (np.average(i['wavelength'][:-1] / np.diff(i['wavelength'])))>1000]

    return spex_sxd

def make_sxd_batch(model_name="marley",model_file="SXD_marley.pkl"):

    spex_sxd = fetch_sxd()

    for name,date in spex_sxd:
        
    h = open('submit_all_marley.sh','w')
    for unum in unums:
        h.write('qsub run_marley_{}.sh\n'.format(unum))

        f = open('run_{}_{}.py'.format(model_name,name),'w')
        f.write('import logging\n')
        f.write('from datetime import date\n\n')
        f.write('from bdmcmc.batch import OneBatch\n\n')
        f.write('logging.basicConfig(level=logging.INFO)\n\n')
        f.write("ob = OneBatch('{}',"
             "'/vega/astro/users/sd2706/modelSpectra/{}',"
             "obs_date={})\n".format(name,model_file,date))
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
        g.write("#PBS -o localhost:/vega/astro/users/sd2706/testing/ \n")
        g.write("#PBS -e localhost:/vega/astro/users/sd2706/testing/ \n")
        g.write("/vega/astro/users/amp2217/yt-x86_64/bin/python"
             " run_{}_{}.py\n".format(model_name,name))
        g.write("# End of script\n")
        g.close()
    h.close()


