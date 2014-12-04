import numpy as np
import matplotlib.pyplot as plt
import asciitable as at
import astropy.units as u

from bdmcmc.spectra import BrownDwarf

figurepath = "/home/stephanie/Dropbox/paperBD/"


res = at.read("/home/stephanie/Dropbox/paperBD/all_results.csv")
field = np.where((res["Grav"]=="f") & (res["SpT"]>=10))[0]
lowg = np.where((res["Grav"]=="l") & (res["SpT"]>=10))[0]

field_names0 = res["Name"][field]
field_names1 = field_names0[np.arange(0,len(field),4)]
lowg_names0 = res["Name"][lowg]
lowg_names1 = lowg_names0[np.arange(0,len(lowg),4)]

field_spt0 = res["SpT"][field]
field_spt1 = field_spt0[np.arange(0,len(field),4)]
lowg_spt0 = res["SpT"][lowg]
lowg_spt1 = lowg_spt0[np.arange(0,len(lowg),4)]

f_sorted = np.argsort(field_spt1)
l_sorted = np.argsort(lowg_spt1)

field_names = field_names1[f_sorted]
field_plotnames0 = res["shortname"][field][np.arange(0,len(field),4)][f_sorted]
field_spt = res["plotspt"][field][np.arange(0,len(field),4)][f_sorted]
field_plotnames = ["{} {}".format(field_plotnames0[i],field_spt[i]) for i in range(len(field_names))]
field_instruments = res["inst"][field][np.arange(0,len(field),4)][f_sorted]
field_obsdates = res["obs_date"][field][np.arange(0,len(field),4)][f_sorted]

lowg_names = lowg_names1[l_sorted]
lowg_plotnames0 = res["shortname"][lowg][np.arange(0,len(lowg),4)][l_sorted]
lowg_spt = res["plotspt"][lowg][np.arange(0,len(lowg),4)][l_sorted]
lowg_plotnames = ["{} {}".format(lowg_plotnames0[i],lowg_spt[i]) for i in range(len(lowg_names))]
lowg_instruments = res["inst"][lowg][np.arange(0,len(lowg),4)][l_sorted]
lowg_obsdates = res["obs_date"][lowg][np.arange(0,len(lowg),4)][l_sorted]


def plot_sample(dbnames,instruments,obs_dates,plotnames=None):
    num_spec = len(dbnames)

    if plotnames==None:
        plotnames=dbnames

    waves = []
    fluxes = []

    for i, name in enumerate(dbnames):
        bd = BrownDwarf(name)
        obs_date = None
        if obs_dates[i]!="-": obs_date = obs_dates[i]
        bd.get_med(instruments[i],obs_date=obs_date)

        waves.append(bd.specs["med"]["wavelength"])
        fluxes.append(bd.specs["med"]["flux"])


    fig = plt.figure(figsize=(10,8))
    ax = fig.add_subplot(121)
    norm_start = num_spec
    print norm_start
    for i, name in enumerate(plotnames):
        if i==(num_spec/2):
            ax1 = ax
            ax = fig.add_subplot(122)
            norm_start = num_spec * 2
            print norm_start
        if (len(fluxes[i])>1) and (len(waves[i])>1):
            good = np.where((np.isnan(fluxes[i])==False) & 
                            (np.isnan(waves[i])==False))[0]
            norm_by = np.median(fluxes[i][good])
            offset = (norm_start - i*2)
            print norm_by,offset,name,instruments[i]
            plot_flux = fluxes[i][good]/norm_by + offset
            ax.step(waves[i][good],plot_flux,color="k")
            textx = 1.9
            texty = np.min((np.max(plot_flux[(waves[i][good]>textx*u.um)]), 
                            offset + 1.5))
            ax.text(textx, texty, plotnames[i])
        else:
            print name,instruments[i]
    ax2 = ax

    ax1.set_ylim(2,num_spec+3)
    ax2.set_ylim(2,num_spec+3)
    ax1.tick_params(labelleft=False)
    ax2.tick_params(labelleft=False)
    ax1.set_xlim(0.9,2.5)
    ax2.set_xlim(0.9,2.5)
    ax1.set_xlabel(r"$\lambda\ $(micron)",fontsize="large")
    ax2.set_xlabel(r"$\lambda\ $(micron)",fontsize="large")
    ax1.set_ylabel("Flux (normalized)",fontsize="large")
    plt.subplots_adjust(wspace=0.0,top=0.95,bottom=0.1,left=0.05,right=0.95)

#plot_sample(field_names,field_instruments,field_obsdates,field_plotnames)
#plt.suptitle("Field Gravity")
#plt.savefig(figurepath+"sample_field.eps",bbox_inches="tight")
plot_sample(lowg_names,lowg_instruments,lowg_obsdates,lowg_plotnames)
plt.suptitle("Low Gravity")
plt.savefig(figurepath+"sample_lowg.eps",bbox_inches="tight")

