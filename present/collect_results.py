import logging, os

import cPickle
import numpy as np

from bdmcmc.spectra import BrownDwarf
from bdmcmc.plotting.spt_plot import spt_plot

#logging.basicConfig(level=logging.DEBUG)

def sort_files(chain_list):
    """ 
    determines the object names and bands that were part of the fit

    parameters
    ----------
    chain_list : string, filename
        file containing filenames and paths for chain files from 
        a set of fits.

    """

    objects = []
    band_names = []
    filenames = []

    f = open(chain_list, "r")

    check = f.readline()
    while check!="":
        this_file = check.split("/")[-1]
        this_file = this_file.replace(" ","_")
        this_file_split = this_file.split("_")
        logging.debug(this_file_split)
        ## if this isn't a chain file, skip it
        if this_file_split[-1]!="chains.pkl\n": 
            check = f.readline()
            continue

        ## If it is a chain file, pull the object and band names
        this_object = this_file_split[0]
        this_band = this_file_split[1]
        if this_band=="Marley":
            this_band = this_file_split[2]
        logging.debug("{} {}".format(this_object,this_band))
        objects.append(this_object)
        band_names.append(this_band)
        filenames.append(check[:-1])
        check = f.readline()

    f.close()

    return filenames, objects, band_names

def get_spts(objects):

    spts = np.zeros(len(objects))
    for i, obj_name in enumerate(objects):
        bd = BrownDwarf(obj_name)
        spts[i] = bd.spt

    return spts

all_extents = {'logg':[2.73,5.75],'fsed':[0.9,5.1],'teff':[1350,2450],
    'ln(s)':[-36,-30]}

model_params = {"Dusty": ["logg","teff"],
                "BTSettl": ["logg","teff"],
                "B06": ["logg","teff"],
                "Marley": ["logg", "fsed", "teff"]
                }

model_titles = {"Dusty":   "Phoenix/DUSTY (Rice et al. 2010)",
                "BTSettl": "BT-Settl (Allers et al. 2011)",
                "B06":    "Burrows et al. (2006)",
                "Marley": "Eddysed (Saumon & Marley 2008)"
                }
def quantile(x,quantiles):
    # From DFM's triangle code
    xsorted = sorted(x)
    qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
    return zip(quantiles,qvalues)

def plot_results(filenames,objects,band_names,spec_types,model_name):

    params = model_params[model_name]
    nparams = len(params)
    print params

    run_names = np.unique(band_names)
    nruns = len(run_names)
    print run_names

    unique_names = np.unique(objects)
    nobj = len(unique_names)
    print unique_names

    medians = np.zeros(nparams*nruns*nobj).reshape(nparams,nruns,nobj)
    errors = np.zeros(nparams*nruns*nobj*2).reshape(nparams,nruns,2,nobj)

    for i,filename in enumerate(filenames):
        infile = open(filename,"rb")
        chains = cPickle.load(infile)
        infile.close()
        cropchain = chains[:,:,:nparams].reshape((-1,nparams))

        quantiles = [quantile(cropchain[:,m],[.05,.16,.5,.84,.95])
                     for m in range(nparams)]
        #print quantiles

        for k, band in enumerate(run_names):
            if band==band_names[i]:
                band_ind = k

        for l, obj in enumerate(unique_names):
            if obj==objects[i]:
                obj_ind = l
        
        for j, pname in enumerate(params):
            print i, filename, obj_ind, unique_names[obj_ind]
            print j, pname, band_ind, run_names[band_ind]
            print quantiles[j]
            up_err = quantiles[j][2][1] - quantiles[j][0][1]
            dn_err = quantiles[j][4][1] - quantiles[j][2][1]
            print up_err,dn_err
            medians[j][band_ind][obj_ind] = quantiles[j][2][1]
            errors[j][band_ind][0][obj_ind] = dn_err
            errors[j][band_ind][1][obj_ind] = up_err

    extents = [all_extents[p] for p in params]

    logging.debug(medians)
    logging.debug(errors)
    logging.debug(extents)

    spt_plot(medians,errors,spec_types,params,run_names,
             model_titles[model_name],single_figure=False,
             extents=extents,model_extents=extents)


#filenames, objects, band_names = sort_files("results_and_chains.lst")
#spec_types = get_spts(np.unique(objects))

#plot_results(filenames,objects,band_names,spec_types,"BTSettl")

#plot_results(filenames,objects,band_names,spec_types,"B06")

#plot_results(filenames,objects,band_names,spec_types,"Marley")
