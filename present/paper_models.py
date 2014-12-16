import logging, os
from datetime import date

import bdmcmc.get_mod, bdmcmc.smooth, bdmcmc.spectra
import numpy as np
import astropy.units as u
import cPickle

logging.basicConfig(level=logging.INFO)

figurepath = "/home/stephanie/Dropbox/paperBD/figures/"

flux_unit=(u.erg / u.cm**2 / u.s / u.um)
modelpath = '/home/stephanie/ldwarfs/modelSpectra/'
#marley = bdmcmc.get_mod.AtmoModel(modelpath+'SXD_r2000_Marley.pkl')
#dusty = bdmcmc.get_mod.AtmoModel(modelpath+'SXD_r2000_Dusty.pkl')
#settl = bdmcmc.get_mod.AtmoModel(modelpath+'SXD_r2000_BTS.pkl')
#settl13 = bdmcmc.get_mod.AtmoModel(modelpath+'SXD_r2000_BTS13.pkl')
# Need to fix wavelength arrays and units


#bd = bdmcmc.spectra.BrownDwarf('U20165')
#bd.get_low()

#marley.model["wsyn"] = bd.specs['low']['wavelength']


mcolors = ["Blue","Magenta","LimeGreen","Red"]
gstyles = ["-","--",":","-."]
names = ["Marley","Gaia-DUSTY","BT-Settl10","BT-Settl13"]


for j, g in enumerate([4.5,5.0,5.5]):
    fig = plt.figure(figsize=(8,3))
    ax1 = fig.add_subplot(141)
    ax2 = fig.add_subplot(142)
    ax3 = fig.add_subplot(143)
    ax4 = fig.add_subplot(144)

    axes = [ax1,ax2,ax3,ax4]

#    for i,T in enumerate([1500,1800,2100]):
    for i,T in enumerate([1400,1700,2000,2300]):
        for k, atmo_model in enumerate([marley,dusty,settl,settl13]):
            ax = axes[k]
            model = atmo_model.model
            if k==0:
                model_loc=np.where((abs(model["teff"]-T)<1) 
                                   & (abs(model["logg"]-g)<0.05)
                                   & (model["fsed"]==1))[0]
            else:
                model_loc = np.where((abs(model["teff"]-T)<1) 
                                   & (abs(model["logg"]-g)<0.05))[0]

            if len(model_loc)==1:
                mflux = model["fsyn"][model_loc]
#                print "model wave!", model["wsyn"]
                mwave = model["wsyn"]
                while len(mflux)==1:
                    mflux = mflux[0]
                while len(mwave)==1:
                    mwave = mwave[0]
            else:
                print "UH OH {} {} {} {}".format(k,T,g,model_loc)
                ax.set_ylim(axes[0].get_ylim())
                continue
#            print "k {} lw {} lf {}".format(k, len(mwave),len(mflux))

            norm_by = np.median(mflux)
            offset = i*1.5
            ax.step(mwave,mflux/norm_by+offset,color=mcolors[k])
                #,ls=gstyles[j])
            #print "{} {} {} {}".format(k,T,g,i*5e-3)

            ax.tick_params(labelleft=False)
            ax.set_xlim(0.85,2.5)
            ax.set_xticks([1,1.5,2,2.5])
            ax.set_xticklabels(["1","1.5","2","2.5"])
            ax.set_xlabel(r"$\lambda$ (micron)",fontsize="large")

            #if k>0:
            #    #print k
            #    ax.set_ylim(axes[0].get_ylim())
            ax.set_ylim((-0.001,8))

    plt.subplots_adjust(wspace=0.0,top=0.95,bottom=0.2,left=0.1,right=0.95)
    axes[0].set_ylabel("Flux (normalized)",fontsize="large")
    axes[0].text(2,5.5,"2300 K",fontsize="small")
    axes[0].text(2,4.,"2000 K",fontsize="small")
    axes[0].text(2,2.75,"1700 K",fontsize="small")
    axes[0].text(2,1.4,r"1400 K",fontsize="small")
    axes[0].text(1.9,7,names[0],color=mcolors[0])
    axes[1].text(1.5,7,names[1],color=mcolors[1])
    axes[2].text(1.6,7,names[2],color=mcolors[2])
    axes[3].text(1.6,7,names[3],color=mcolors[3])
    plt.savefig(figurepath+"models_sxd_g{}.eps".format(g*10),bbox_inches="tight")
    #break

"""
fig = plt.figure(figsize=(8,3))
ax1 = fig.add_subplot(141)
ax2 = fig.add_subplot(142)
ax3 = fig.add_subplot(143)
ax4 = fig.add_subplot(144)

axes = [ax1,ax2,ax3,ax4]

tcolors = ["Red","Orange","Green","Blue"]

g = 5.0
model = marley.model
for i,T in enumerate([1400,1600,1800,2000]):
    for j, F in enumerate([1,2,3,4]):
        print i, T, j, F , g
        ax = axes[j]
        model_loc = np.where((abs(model["teff"]-T)<1) 
                             & (abs(model["fsed"]-F)<0.1)
                             & (abs(model["logg"]-g)<0.05))[0]
        if len(model_loc)==1:
            mflux = model["fsyn"][model_loc]
            #print model["wsyn"][model_loc]
            mwave = model["wsyn"][model_loc]
            while len(mflux)==1:
                mflux = mflux[0]
            while len(mwave)==1:
                mwave = mwave[0]
            #print "k {} lw {} lf {}".format(k, len(mwave),len(mflux))
        else:
            print "UH OH {} {} {}".format(T,g,model_loc)
            ax.set_ylim(axes[0].get_ylim())
            continue

        mflux2 = bdmcmc.smooth.smooth_model(mwave,mflux,
            bd.specs["low"]["wavelength"],
            np.average(np.diff(bd.specs["low"]["wavelength"]*2))*u.um)
        mflux = mflux2
        mwave = bd.specs["low"]["wavelength"]

        norm_by = np.sum(mflux[abs(mwave-1.27*u.um)<(0.005*u.um)])
        ax.step(mwave,mflux/norm_by+(4-i),color=tcolors[i])

        ax.tick_params(labelleft=False)
        ax.set_xticklabels(["","1","1.5","2","2.5"])
        ax.set_xlabel(r"$\lambda (\AA)$",fontsize="large")
        if j>0:
            #print k
            ax.set_ylim(axes[0].get_ylim())

plt.subplots_adjust(wspace=0.0,top=0.95,bottom=0.2,left=0.1,right=0.95)
axes[0].set_ylabel("Flux (normalized)",fontsize="large")
plt.savefig("paper_fsed_ex.eps",bbox_inches="tight")
"""
