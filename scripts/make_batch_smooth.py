import logging

import numpy as np
import asciitable as at
import cPickle

import bdmcmc.get_mod, bdmcmc.spectra

logging.basicConfig(level=logging.INFO)

def make_batch_smooth(model_name,model_file):
    modelpath = '/vega/astro/users/sd2706/modelSpectra/'
    #modelpath = '/home/stephanie/ldwarfs/modelSpectra/'
    am = bdmcmc.get_mod.AtmoModel(modelpath+model_file)
    num_jobs = len(am.model["teff"])/10 + 1

    
    h = open('submit_{}_smooth.sh'.format(model_name),'w')
    for i in range(num_jobs):
        h.write('qsub variable_smooth_{}{}.sh\n'.format(model_name,i))
    
        g = open('variable_smooth_{}{}.sh'.format(model_name,i),'w')
        g.write("#!/bin/sh\n\n")
        g.write("# Directives\n")
        g.write("#PBS -N Smooth{}{}\n".format(model_name,i))
        g.write("#PBS -W group_list=yetiastro\n")
        g.write("#PBS -l nodes=1,walltime=24:00:00,mem=3000mb\n")
        g.write("#PBS -M sd2706@columbia.edu \n")
        g.write("#PBS -m ae\n")
        g.write("#PBS -V\n\n")
        g.write("# Set output and error directories\n")
        g.write("#PBS -o localhost:/vega/astro/users/sd2706/testing/ \n")
        g.write("#PBS -e localhost:/vega/astro/users/sd2706/testing/ \n")
        g.write("/vega/astro/users/amp2217/yt-x86_64/bin/python "
                " variable_smooth_{}{}.py\n".format(model_name,i))
        g.write("# End of script\n")
        g.close()

        f = open('variable_smooth_{}{}.py'.format(model_name,i),'w')
        f.write("import logging\nfrom datetime import date\n")
        f.write("import bdmcmc.spectra,bdmcmc.get_mod,bdmcmc.smooth\n")
        f.write("import numpy as np\nimport astropy.units as u\n")
        f.write("from scipy.io.idl import readsav\nimport cPickle\n\n")
        f.write("logging.basicConfig(level=logging.INFO)\n\n")
        f.write("am = bdmcmc.get_mod.AtmoModel('{}{}')\n".format(
                modelpath,model_file))
        f.write("logging.debug(str(am.model['fsyn'][0]))\n\n")
        f.write("bd = bdmcmc.spectra.BrownDwarf('U20165')\n")
        f.write("bd.get_low()\n")
        f.write("data_wave = bd.specs['low']['wavelength']\n")
        f.write("data_flux = bd.specs['low']['flux']\n")
        f.write("logging.info('got bd')\n\n")
        f.write("r_scale = 0.3/bd.specs['low']['slit_width'].value\n")
        f.write("logging.info(str(r_scale))\n\n")
        f.write("sub_grid = am.model.copy()\n")
        f.write("for k in sub_grid.keys():\n")
        f.write("    sub_grid[k] = sub_grid[k][{}:{}]\n".format(i*10,(i+1)*10))
        f.write("new_grid = bdmcmc.smooth.smooth_grid(sub_grid,data_wave, "
                "res_scale=r_scale,incremental_outfile="
                "'{}_backup{}.pkl')\n\n".format(model_name,i))
        f.write("output = open('{}_output{}.pkl','wb')\n".format(model_name,i))
        f.write("cPickle.dump(new_grid,output)\n")
        f.write("output.close()\n")

    h.close()

def recombine_batches(filelist,original_model_file,output_file,data_wave,
    data_flux_units):

    modelpath = '/vega/astro/users/sd2706/modelSpectra/'
    #modelpath = '/home/stephanie/ldwarfs/modelSpectra/'

    sub_files0 = at.read(modelpath+"smoothfiles/"+filelist)
    sub_files = sub_files0["filenames"]
    print sub_files

    am = bdmcmc.get_mod.AtmoModel(modelpath+original_model_file)
    num_jobs = len(am.model["teff"])/10 + 1

    if num_jobs!=len(sub_files):
        print "!!! {} files but should have {}".format(num_jobs,len(sub_files))
        return
    else:
        logging.info("expected {} files, got {}!".format(num_jobs,len(sub_files)))

    orig_params = np.asarray(am.params)
    orig_params = np.delete(orig_params,np.where(orig_params=="wsyn"))
    orig_params = np.delete(orig_params,np.where(orig_params=="fsyn"))

    for i,sfile in enumerate(sub_files):
        print i, sfile
        infile = open(modelpath+"smoothfiles/"+sfile)
        sub_grid = cPickle.load(infile)
        infile.close()

        sub_params = np.asarray(sub_grid.keys())
        sub_params = np.delete(sub_params,np.where(sub_params=="wsyn"))
        sub_params = np.delete(sub_params,np.where(sub_params=="fsyn"))
        logging.info(sub_params)
        logging.info(sub_grid["teff"])

        for j in range(len(sub_grid["teff"])):
#            logging.debug("{} {}".format(i,j))
            same_params = True

            for k,param in enumerate(sub_params):
                if (am.model[param][10*i+j]!=sub_grid[param][j]):
                    same_params=False

            if same_params:
                for k,param in enumerate(sub_params):
                    logging.debug('!GOOD {} {} {} {} {}'.format(i,j,param,
                        am.model[param][10*i+j],sub_grid[param][j]))

                am.model['fsyn'][10*i + j] = sub_grid['fsyn'][j]
            else:
                for k,param in enumerate(sub_params):
                    logging.debug('OH NO {} {} {} {} {}'.format(i,j,param,
                        am.model[param][10*i+j],sub_grid[param][j]))
            
    logging.info("types fsyn {} fsyn[0] {}".format(type(am.model["fsyn"]),
        type(am.model["fsyn"][0])))
    logging.info("shapes f {} w {}".format(np.shape(am.model["fsyn"]),
        np.shape(am.model["wsyn"])))

    am.model["wsyn"] = data_wave

    new_grid = am.model.copy()

    another_fsyn = []
    for i in range(len(new_grid['logg'])):
        logging.debug(i)
        logging.debug(type(new_grid['fsyn'][i]))
        converted_fsyn = new_grid['fsyn'][i].to(data_flux_units,
            equivalencies=u.spectral_density(data_wave))
        logging.debug("converted! {}".format(converted_fsyn.unit))
        another_fsyn = np.append(another_fsyn,converted_fsyn).reshape((i+1,-1))
        logging.debug(another_fsyn[0])
    new_grid['fsyn'] = np.array(another_fsyn)*spectrum['flux'].unit

    logging.info("final w {} f {}".format(new_grid['wsyn'].unit,
                                          new_grid['fsyn'].unit))

    outfile = open(modelpath+output_file,'wb')
    cPickle.dump(new_grid,outfile)
    outfile.close()


#make_batch_smooth("Marley","marley_ldwarfs_all.pkl")
#make_batch_smooth("BTS","btsettl_wide.pkl")
#make_batch_smooth("Bur","burrows_06_cloud_expanded.pkl")

bd = bdmcmc.spectra.BrownDwarf('U20165')
bd.get_low()
data_wave = bd.specs['low']['wavelength']
data_flux_units = bd.specs['low']['flux'].unit
recombine_batches("BTS_subgrids.lst","btsettl_wide.pkl","SpeX_BTS_wide.pkl",
                  data_wave,data_flux_units)
