import logging, os

import numpy as np
import astropy.units as u
import cPickle

import bdmcmc.spectra, bdmcmc.get_mod, bdmcmc.make_model

reload(bdmcmc.spectra)
reload(bdmcmc.make_model)

lformat = "%(asctime)s %(message)s"
logging.basicConfig(level=logging.INFO,format=lformat)

mbase_path = '/vega/astro/users/sd2706/'
if os.path.exists(mbase_path)==False:
    mbase_path = '/home/stephanie/ldwarfs/'

bd = bdmcmc.spectra.BrownDwarf('U20012')
bd.get_low()

am = bdmcmc.get_mod.AtmoModel(mbase_path+'modelSpectra/SpeX_marley.pkl')

test_model_p = [4.794,2.34,2139.]
test_norm_p = [1e-16,2e-16,1.4e-16]
test_s = 5e-16
test_p = np.append(np.append(test_model_p,test_norm_p),test_s)

logging.info("testing WITH snap")
mg2 = bdmcmc.make_model.ModelGrid(bd.specs['low'],am.model,am.params,snap=True)
print mg2(test_p)

infile = open(mbase_path+"testing/2057-0252_Marley_full_2014-06-23_chains.pkl","rb")
chains = cPickle.load(infile)
infile.close()

logging.info("chain shape {}".format(np.shape(chains)))

cropchain = chains.reshape((-1,np.shape(chains)[-1]))
logging.info("cropchain shape {}".format(np.shape(cropchain)))

new_cropchain = mg2.snap_full_run(cropchain)
logging.info("new cropchain shape {}".format(np.shape(new_cropchain)))

######################################################
logging.info("testing snap with SXD")
bd = bdmcmc.spectra.BrownDwarf('0355+1133')
bd.get_low(obs_date='2007-11-13')

am = bdmcmc.get_mod.AtmoModel(mbase_path+'modelSpectra/SXD_marley.pkl')

plot_title="Test {}".format(bd.shortname)
bdsamp = bdmcmc.bdfit.BDSampler(bd.name,bd.specs['low'],am.model,
    am.params,smooth=False,plot_title=plot_title)

logging.info("set up BDSampler")
bdsamp.mcmc_go(nwalk_mult=2,nstep_mult=10)






"""
new_cropchain = np.copy(cropchain)

#This takes about 5-6 minutes for 500,000 links
logging.info("Starting to snap chains")
for j,p in enumerate(cropchain[0:20]):
#for j,p in enumerate(cropchain):
    logging.debug('starting params %s',str(p))

    # p_loc is the location in the model grid that fits all the 
    # constraints up to that point. There aren't constraints yet,
    # so it matches the full array.
    p_loc = range(len(mg2.model['fsyn']))

    for i in range(mg2.ndim):
        # find the location in the ith parameter array corresponding
        # to the ith parameter
        this_p_loc = mg2.find_nearest(mg2.model[mg2.params[i]],p[i])

        # then match that with the locations that already exist
        p_loc = np.intersect1d(p_loc,this_p_loc)

    for i in range(mg2.ndim):
        new_cropchain[j,i] = mg2.model[mg2.params[i]][p_loc]
    logging.debug("snapped params {}".format(new_cropchain[j]))
logging.info("Finished snapping chains")

# If the grid is even all around, run this for each model parameter
# This WILL NOT work when the grid isn't even
# def myround(x,base):
#    return np.round(base*np.round(x*1.0/base),2)
# myround(cropchain[:,0],0.1)


new_chains = new_cropchain.reshape(np.shape(chains))
"""
