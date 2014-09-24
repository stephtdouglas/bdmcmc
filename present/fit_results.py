import logging, os
from datetime import date

import numpy as np
import astropy.units as u
import asciitable as at
import matplotlib.pyplot as plt
import cPickle
import astropy.units as u

import bdmcmc.get_mod, bdmcmc.spectra
import bdmcmc.make_model
import bdmcmc.mask_bands as mb
from bdmcmc.plotting.plot_random import plot_random
from bdmcmc.plotting.spt_plot import spt_plot
from bdmcmc.plotting import triangle
from bdmcmc.present import collect_results

logging.basicConfig(level=logging.INFO)

def quantile(x,quantiles):
    # From DFM's triangle code
    xsorted = sorted(x)
    qvalues = [xsorted[int(q * len(xsorted))] for q in quantiles]
    return zip(quantiles,qvalues)


figurepath = "/home/stephanie/Dropbox/paperBD/figures/"

modelpath = '/home/stephanie/ldwarfs/modelSpectra/'

all_extents = {'logg':[2.75,6.25],'fsed':[0.9,5.1],'teff':[1150,3050],
    'ln(s)':[-36,-30]}

model_params = {"Gaia-DUSTY": ["logg","teff"],
                "BT-Settl": ["logg","teff"],
                "Burrows": ["logg","teff"],
                "Marley": ["logg", "fsed", "teff"]
                }

mcolors = ["Blue","Magenta","LimeGreen","Orange"]
mod_names = ["Marley","Gaia-DUSTY","BT-Settl","Burrows"]
mod_file_names = ["Marley","Dusty","BT-Settl","B06"]
band_names = ["full","J","H","K"]
num_params = {"Marley":3,"Gaia-DUSTY":2,"BT-Settl":2,"Burrows":2}
plot_extents = {}
for i, mname in enumerate(mod_names):
    plot_extents[mname] = [all_extents[p] for p in model_params[mname]]
model_extents = {"Marley":[[4.5,5.5],[1,5],[1500,2400]],
                 "Gaia-DUSTY":[[3.0,5.5],[1400,2400]],
                 "BT-Settl":[[2.5,6.0],[1200,3000]],
                 "Burrows":[[4.5,5.5],[700,2300]]}
sxd_extents = {"Marley":[[4.,5.5],[1,5],[1200,2400]],
                 "Gaia-DUSTY":[[3.0,6.0],[1400,2400]],
                 "BT-Settl":[[2.5,6.0],[1200,3000]],
                 "Burrows":[[4.5,5.5],[700,2300]]}

fitpath = "/home/stephanie/ldwarfs/batch_ldwarfs/"
#fitfolders = ["Marley_2014-09-15_med","Dusty_2014-09-17_med/",
#              "BTSettl_2014-09-16_med",""]

def one_reg(ax, cropchain, model,rand_color):
    rs = plot_random(cropchain, model, ax,rand_color,False)

def plot_four(bd,chain_filename,plot_title,model,axes_row,
    rand_color):

    am = model


    # Set up bands
    wav = bd.specs['low']['wavelength']
    full = np.where((wav>=0.9*u.um) & (wav<2.5*u.um))[0]
    Jband = np.where((wav>=0.9*u.um) & (wav<1.4*u.um))[0]
    Hband = np.where((wav>=1.4*u.um) & (wav<1.9*u.um))[0]
    Kband = np.where((wav>=1.9*u.um) & (wav<2.5*u.um))[0]

    bands = {'J':Jband,'H':Hband,'K':Kband,'full':full}


    for i,b in enumerate(band_names):
        this_chainfile = chain_filename.replace('full',b)

        if os.path.exists(this_chainfile)==False:
            logging.info("{} does not exist! Skipping.".format(this_chainfile))
            continue

        logging.debug("NOW PLOTTING {} {}".format(bd.shortname,b))
        band_spectrum = {
            'wavelength':bd.specs['low']['wavelength'][bands[b]],
            'flux':bd.specs['low']['flux'][bands[b]],
            'unc':bd.specs['low']['unc'][bands[b]]}

        print band_spectrum["wavelength"][0:10]
        print band_spectrum["flux"][0:10]
        print band_spectrum["unc"][0:10]

        mg = bdmcmc.make_model.ModelGrid(band_spectrum,am.model,
            am.params,smooth=False,snap=True)

        cpfile = open(this_chainfile,'rb')
        chains = cPickle.load(cpfile)
        cpfile.close()

        cropchain = chains.reshape((-1,np.shape(chains)[-1]))
        
        one_reg(axes_row[i],cropchain,mg,rand_color)

def plot_all_results(object_names,sample_name,models,mcolors,names,
    mod_file_names,fitfolders,dates):

    fig = plt.figure(figsize=(8,3.5))
    for j,bd_name in enumerate(object_names):
        if bd_name=="1506+1321": continue

        mask_H=True
        # Set up brown dwarf
        bd = bdmcmc.spectra.BrownDwarf(bd_name)
        bd.get_low()#obs_date=obs_date)
        if mask_H:
            mask = mb.BandMask(bd.specs['low']['wavelength'])
            mask.mask_Hband()
            mask.make_pixel_mask()

            reverse_mask = np.delete(np.arange(len(bd.specs['low']['wavelength'])),
                mask.pixel_mask)
            mask.pixel_mask = reverse_mask
    
            bd.specs['low']['wavelength'] = bd.specs['low']['wavelength'][
                mask.pixel_mask]
            bd.specs['low']['flux'] = bd.specs['low']['flux'][mask.pixel_mask]
            bd.specs['low']['unc'] = bd.specs['low']['unc'][mask.pixel_mask]

        ## number of rows of subplots depends on the number of models
        sp_array = (len(models),5)

        for i,mod in enumerate(models):
            row_of_axes = [plt.subplot2grid(sp_array,(i,0),colspan=2),
                           plt.subplot2grid(sp_array,(i,2)),
                           plt.subplot2grid(sp_array,(i,3)),
                           plt.subplot2grid(sp_array,(i,4))]

            chain_filename = "{}{}_full_{}_{}_chains.pkl".format(fitfolders[i],
                              bd_name,mod_file_names[i],dates[i])
            if (i==0) and (sample_name=="prism"):
                chain_filename = "{}{} {} full {}_chains.pkl".format(
                              fitfolders[i],bd_name,mod_file_names[i],dates[i])
            plot_four(bd,fitpath+chain_filename,"",mod,row_of_axes,mcolors[i])


            if i==(len(models)-1):
                for ax in row_of_axes:
                    ax.set_yticklabels([])
                    ax.tick_params(labelsize="large",labelleft=False)
                    ax.set_ylabel("")
                row_of_axes[0].set_ylabel(r"F$_{\lambda}$")
                row_of_axes[0].set_xticklabels(["0.8","","1.2","","1.6","","2.0","","2.4"])
                row_of_axes[1].set_xticklabels(["","1.0","","1.2","","1.4"])
                row_of_axes[2].set_xticklabels(["","","1.6","","1.8"])
                row_of_axes[3].set_xticklabels(["","2.0","","2.2","","2.4"])
                row_of_axes[0].set_xlabel("")
                row_of_axes[2].set_xlabel("")
                row_of_axes[3].set_xlabel("")
            else:
                for ax in row_of_axes:
                    ax.set_yticklabels([])
                    ax.tick_params(labelleft=False,labelbottom=False)
                    ax.set_ylabel("")
                    ax.set_xlabel("")
                row_of_axes[0].set_ylabel(r"F$_{\lambda}$")
    
        plt.subplots_adjust(wspace=0,hspace=0,left=0.06,right=0.99,top=0.95,
            bottom=0.17)

        plt.savefig("{}_{}_results.eps".format(bd_name,sample_name),
            bbox_inches="tight")
        plt.savefig("{}_{}_results.png".format(bd_name,sample_name),
            bbox_inches="tight")
        plt.clf()

def tabulate_all_results(object_names,sample_name,names,modfilenames,
    fitfolders,dates):

    num_cols = 1
    for mname in names:
        num_cols += num_params[mname]
    cols = "l"*num_cols

    f = open("/home/stephanie/Dropbox/paperBD/{}_results.tbl".format(sample_name),
             "w")

    f.write('\\begin{deluxetable}{%s}[!t]\n' % (cols))
    f.write('%\\rotate\n')
    f.write('\\tablewidth{0pt}\n')
    f.write('\\tabletypesize{\\scriptsize}\n')
    f.write('\\tablecaption{ }\n')
    f.write('\\tablehead{\n')
    f.write("\\colhead{} & "*(num_cols-1))
    f.write("\\colhead{} \\\\\n}\n\\startdata\n")

    for j,bd_name in enumerate(object_names):
        f.write(" \\\\ \n \\tableline \\multicolumn{{{}}}{{c}}{{{}}} ".format(
                num_cols,bd_name))
        for k,b in enumerate(band_names):
            if k==0:
                f.write(" \\\\ \\tableline \n {} ".format(b))
            else:
                f.write("\\\\ \n {} ".format(b))

            for i,mod in enumerate(modfilenames):
                this_chainfile = fitpath+"{}{}_{}_{}_{}_chains.pkl".format(
                    fitfolders[i],bd_name,b,modfilenames[i],dates[i])
                if (i==0) and (sample_name=="prism"):
                   this_chainfile  = fitpath+"{}{} {} {} {}_chains.pkl".format(
                        fitfolders[i],bd_name,modfilenames[i],b,dates[i])
                if os.path.exists(this_chainfile)==False:
                    logging.info("{} does not exist! Skipping.".format(
                                 this_chainfile))
                    f.write("&  "*num_params[names[i]])
                    continue
                cpfile = open(this_chainfile,'rb')
                chains = cPickle.load(cpfile)
                cpfile.close()

                cropchain = chains.reshape((-1,np.shape(chains)[-1]))

                quantiles = [quantile(cropchain[:,m],[.05,.5,.95])
                     for m in range(num_params[names[i]])]
                for p in range(num_params[names[i]]):
                    dn_err = quantiles[p][1][1] - quantiles[p][0][1]
                    up_err = quantiles[p][2][1] - quantiles[p][1][1]
                    median = quantiles[p][1][1]
                    f.write(" & ${0:.1f}_{{ {1:.1f} }}^{{ {2:.1f} }}$ ".format(
                            median, dn_err, up_err))
#                    if (p==(num_params[names[i]]-1)) and (k==(len(band_names)-1)):
#                         f.write(" \\\\ \n")

    f.write('\\enddata\n\\end{deluxetable}\n')
    f.close()

def spt_plot_all(object_names,spts,sample_name,mcolors,names,
    modfilenames,fitfolders,dates,mod_extents):

    for i,mod in enumerate(modfilenames):
        logging.info(mod)
        nruns = len(band_names)
        nparams = num_params[names[i]]
        nobj = len(object_names)

        medians = np.zeros(nparams*nruns*nobj).reshape(nparams,nruns,nobj)
        errors = np.zeros(nparams*nruns*nobj*2).reshape(nparams,nruns,2,nobj)
#        dn_errors = np.zeros(nparams*nruns*nobj).reshape(nparams,nruns,nobj)
#        up_errors = np.zeros(nparams*nruns*nobj).reshape(nparams,nruns,nobj)
        for j,bd_name in enumerate(object_names):
            logging.info(bd_name)
            for k,b in enumerate(band_names):
                this_chainfile = fitpath+"{}{}_{}_{}_{}_chains.pkl".format(
                    fitfolders[i],bd_name,b,modfilenames[i],dates[i])
                if (i==0) and (sample_name=="prism"):
                   this_chainfile  = fitpath+"{}{} {} {} {}_chains.pkl".format(
                        fitfolders[i],bd_name,modfilenames[i],b,dates[i])
                if os.path.exists(this_chainfile)==False:
                    logging.info("{} does not exist! Skipping.".format(
                                 this_chainfile))
                    continue
                cpfile = open(this_chainfile,'rb')
                chains = cPickle.load(cpfile)
                cpfile.close()
                cropchain = chains.reshape((-1,np.shape(chains)[-1]))

                quantiles = [quantile(cropchain[:,m],[.05,.5,.95])
                     for m in range(num_params[names[i]])]
                for p in range(num_params[names[i]]):
                    dn_err = quantiles[p][1][1] - quantiles[p][0][1]
                    up_err = quantiles[p][2][1] - quantiles[p][1][1]
                    medians[p][k][j] = quantiles[p][1][1]
                    errors[p][k][0][j] = dn_err
                    errors[p][k][1][j] = up_err
                    logging.debug("{} {} {} {}".format(p,b,j,errors[p][k]))
#                    dn_errors[p][k][j] = dn_err
#                    up_errors[p][k][j] = up_err

#        errors = [[[dn_errors[pp][kk],up_errors[pp][kk]] 
#                  for kk in range(len(band_names))]
#                  for pp in range(num_params[names[i]])]

        logging.debug(medians)
        logging.debug(errors)

        spt_plot(medians,errors,spts,model_params[names[i]],band_names,
                 "{}_{}".format(sample_name,mod),single_figure=False,
                 extents=plot_extents[names[i]],
                 model_extents=mod_extents[names[i]])
        
def spt_plot_models(object_names,spts,sample_name,mcolors,names,
    modfilenames,fitfolders,dates,mod_extents):

    for k,b in enumerate(band_names):
        logging.info(b)
        nruns = len(modfilenames)
        nparams = 2
        nobj = len(object_names)

        medians = np.zeros(nparams*nruns*nobj).reshape(nparams,nruns,nobj)
        errors = np.zeros(nparams*nruns*nobj*2).reshape(nparams,nruns,2,nobj)
        for j,bd_name in enumerate(object_names):
            logging.info(bd_name)
            for i,mod in enumerate(modfilenames):
                logging.info(mod)
                this_chainfile = fitpath+"{}{}_{}_{}_{}_chains.pkl".format(
                    fitfolders[i],bd_name,b,modfilenames[i],dates[i])
                if (i==0) and (sample_name=="prism"):
                   this_chainfile  = fitpath+"{}{} {} {} {}_chains.pkl".format(
                        fitfolders[i],bd_name,modfilenames[i],b,dates[i])
                if os.path.exists(this_chainfile)==False:
                    logging.info("{} does not exist! Skipping.".format(
                                 this_chainfile))
                    continue
                cpfile = open(this_chainfile,'rb')
                chains = cPickle.load(cpfile)
                cpfile.close()
                cropchain = chains.reshape((-1,np.shape(chains)[-1]))

                quantiles = [quantile(cropchain[:,m],[.05,.5,.95])
                     for m in range(num_params[names[i]])]
                for p in range(num_params[names[i]]):
                    if (names[i]=="Marley") and (p==1): #fsed
                        continue
                    elif (names[i]=="Marley") and (p==2): #teff
                        me_loc = p-1
                        dn_err = quantiles[p][1][1] - quantiles[p][0][1]
                        up_err = quantiles[p][2][1] - quantiles[p][1][1]
                        medians[me_loc][i][j] = quantiles[p][1][1]
                        errors[me_loc][i][0][j] = dn_err
                        errors[me_loc][i][1][j] = up_err
                    else:
                        dn_err = quantiles[p][1][1] - quantiles[p][0][1]
                        up_err = quantiles[p][2][1] - quantiles[p][1][1]
                        medians[p][i][j] = quantiles[p][1][1]
                        errors[p][i][0][j] = dn_err
                        errors[p][i][1][j] = up_err
                        logging.info(medians[p])
                        logging.info(errors[p])

        logging.info("DONE WITH {}".format(b))
        logging.info(medians)
        logging.info(errors)
        spt_plot(medians,errors,spts,["logg","teff"],names,
                 "{}_{}".format(sample_name,b),single_figure=False,
                 extents=plot_extents[names[i]],
                 model_extents=[[2.5,6.0],[1200,3000]],run_colors=mcolors)
        
def plot_marley_correlations(object_names,spts,sample_name,name,
    modfilename,fitfolder,date,mod_extents):
    print "hello"

    plt.figure(figsize=(4,4))
    for j,bd_name in enumerate(object_names):
        for k,b in enumerate(band_names):
            this_chainfile = fitpath+"{}{} {} {} {}_chains.pkl".format(
                        fitfolder,bd_name,modfilename,b,date)
            if os.path.exists(this_chainfile)==False:
                logging.info("{} does not exist! Skipping.".format(
                             this_chainfile))
                continue

            cpfile = open(this_chainfile,'rb')
            chains = cPickle.load(cpfile)
            cpfile.close()
            cropchain = chains.reshape((-1,np.shape(chains)[-1]))

            bins=50
            x = cropchain[:,1] #fsed
            cmap = cm.get_cmap("gray")
            params = [r"$log(g)$","fsed",r"$T_{eff}$"]
            for i in (0,2):
                y = cropchain[:,i] #logg=0, teff=2
                extent = [[x.min(), x.max()], [y.min(), y.max()]]
                ax = plt.subplot(111)
                X = np.linspace(extent[0][0], extent[0][1], bins + 1)
                Y = np.linspace(extent[1][0], extent[1][1], bins + 1)
                H, X, Y = np.histogram2d(x.flatten(), y.flatten(), bins=(X, Y))
                V = 1.0 - np.exp(-0.5 * np.arange(0.5, 2.1, 0.5) ** 2)
                Hflat = H.flatten()
                inds = np.argsort(Hflat)[::-1]
                Hflat = Hflat[inds]
                sm = np.cumsum(Hflat)
                sm /= sm[-1]
                for ii, v0 in enumerate(V):
                    try:
                        V[ii] = Hflat[sm <= v0][-1]
                    except:
                        V[ii] = Hflat[0]

                X1, Y1 = 0.5 * (X[1:] + X[:-1]), 0.5 * (Y[1:] + Y[:-1])
                X, Y = X[:-1], Y[:-1]
                ax.pcolor(X, Y, H.max() - H.T, cmap=cmap)
                ax.contour(X1, Y1, H.T, V, colors="r", linewidths=1)
                ax.set_xlabel(r"$f_{sed}$",fontsize="x-large")
                ax.set_title("{} {}".format(bd_name,spts[j]))
                print i
                ax.set_ylabel(params[i],fontsize="x-large")
                plt.savefig("{}_{}{}_corr.eps".format(bd_name,b,i),
                            bbox_inches="tight")
                plt.savefig("{}_{}{}_corr.png".format(bd_name,b,i),
                            bbox_inches="tight")
                plt.clf()
#            break
#        break

    

############################# 

#all_filenames, all_names, all_band_names = collect_results.sort_files(
#    fitpath+"BTSettl_2014-09-15/results_and_chains.lst")

#prism_names = np.unique(all_names)
#prism_types = collect_results.get_spts(prism_names)

prism_sample = at.read("/home/stephanie/Dropbox/paperBD/prism_sample.csv",
                     delimiter="\t")
ra_order = np.argsort(prism_sample["Shortname"])
prism_names = prism_sample["Shortname"][ra_order]
prism_names[-1] = "U12128"
prism_types = prism_sample["Ntype"][ra_order]
prism_spectypes = prism_sample["Type"][ra_order]


#burrows = bdmcmc.get_mod.AtmoModel(modelpath+'SpeX_B06_wide.pkl')
#marley = bdmcmc.get_mod.AtmoModel(modelpath+'SpeX_marley.pkl')
#dusty = bdmcmc.get_mod.AtmoModel(modelpath+'SpeX_dusty_cut.pkl')
#settl = bdmcmc.get_mod.AtmoModel(modelpath+'SpeX_BTS_wide.pkl')
#prism_models=[marley,dusty,settl,burrows]

prism_fitfolders = ["Marley_2014-06-23/","Dusty_2014-09-17/",
              "BTSettl_2014-09-15/","Burrows_2014-09-16/"]
prism_dates = [prism_fitfolders[i].split("_")[1][:-1] for i in range(4)]

#plot_marley_correlations(prism_names,prism_spectypes,"prism",
#    mod_names[0],mod_file_names[0],prism_fitfolders[0],prism_dates[0],
#    model_extents["Marley"])
#plot_all_results(prism_names,"prism",prism_models,mcolors,mod_names,
#    mod_file_names,prism_fitfolders,prism_dates)
#plt.close("all")

#tabulate_all_results(prism_names,"prism",
#    mod_names,mod_file_names,prism_fitfolders,prism_dates)

#spt_plot_all(prism_names,prism_types,"prism",mcolors,
#    mod_names,mod_file_names,prism_fitfolders,prism_dates,model_extents)
#plt.close("all")

#spt_plot_models(prism_names,prism_types,"prism",mcolors,
#    mod_names,mod_file_names,prism_fitfolders,prism_dates,model_extents)
#plt.close("all")



### SXD ####

#all_filenames, all_names, all_band_names = collect_results.sort_files(
#    fitpath+"Dusty_2014-09-17_med/results_and_chains.lst")
#
#sxd_names = np.unique(all_names)
#sxd_types = collect_results.get_spts(sxd_names)


sxd_sample = at.read("/home/stephanie/Dropbox/paperBD/SXD_table.csv",
                     delimiter="\t")
ra_order = np.argsort(sxd_sample["shortname"])
sxd_names = sxd_sample["shortname"][ra_order]#[:6]
sxd_types = sxd_sample["SpT"][ra_order]#[:6]


#smarley = bdmcmc.get_mod.AtmoModel(modelpath+'SXD_r2000_Marley.pkl')
#sdusty = bdmcmc.get_mod.AtmoModel(modelpath+'SXD_r2000_Dusty.pkl')
#ssettl = bdmcmc.get_mod.AtmoModel(modelpath+'SXD_r2000_BTS.pkl')
#sxd_models=[smarley,sdusty,ssettl]

sxd_fitfolders = ["Marley_2014-09-15_med/","Dusty_2014-09-17_med/",
              "BTSettl_2014-09-16_med/"]
sxd_dates = [sxd_fitfolders[i].split("_")[1] for i in range(3)]
sxd_mod_file_names = ["Marley","Dusty","BTSettl"]

#plot_all_results(sxd_names,"SXD",sxd_models,mcolors[:-1],mod_names[:-1],
#    sxd_mod_file_names,sxd_fitfolders,sxd_dates)
#plt.close("all")

#tabulate_all_results(sxd_names,"SXD",mod_names[:-1],
#    sxd_mod_file_names,sxd_fitfolders,sxd_dates)


#spt_plot_all(sxd_names,sxd_types,"SXD",mcolors[:-1],
#    mod_names[:-1],sxd_mod_file_names,sxd_fitfolders,sxd_dates,
#    sxd_extents)
#plt.close("all")

#spt_plot_models(sxd_names,sxd_types,"SXD",mcolors[:-1],
#    mod_names[:-1],sxd_mod_file_names,sxd_fitfolders,sxd_dates,
#    sxd_extents)
#plt.close("all")


bd_name = "0033-1521"
i=2 #BT-Settl
b="K"
chain_filename = fitpath+"{}{}_{}_{}_{}_chains.pkl".format(sxd_fitfolders[i],bd_name,b,sxd_mod_file_names[i],sxd_dates[i])

cpfile = open(chain_filename,'rb')
chains = cPickle.load(cpfile)
cpfile.close()

cropchain = chains.reshape((-1,np.shape(chains)[-1]))

fig,axes = triangle.corner(cropchain,[r"$log(g)$",r"$T_{eff}$",r"$N_{K}$",r"$ln(s)$"],plot_datapoints=False,extents=[[4.8,5.7],[1800,2050],[0.986,1.016],[-30,-29.5]])#,plot_contours=True)

plt.subplots_adjust(hspace=0.05,wspace=0.05,bottom=0.12,left=0.11,top=0.98,right=0.95)
for ax in axes.flatten():
    ax.set_xlabel(ax.get_xlabel(),fontsize="x-large")
    ax.set_ylabel(ax.get_ylabel(),fontsize="x-large")

plt.savefig("corner_ex_{}_{}.png".format(bd_name,sxd_mod_file_names[i]))
plt.savefig("corner_ex_{}_{}.eps".format(bd_name,sxd_mod_file_names[i]))

