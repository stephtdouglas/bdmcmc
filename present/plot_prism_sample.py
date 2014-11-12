import logging

import numpy as np
import matplotlib.pyplot as plt
import asciitable as at

from bdmcmc.spectra import BrownDwarf

logging.basicConfig(level=logging.INFO)

figurepath = "/home/stephanie/Dropbox/paperBD/figures/"

prism_sample = at.read("prism_sample.csv",delimiter="\t")

waves = []
fluxes = []

for i, unum in enumerate(prism_sample["unum"]):
    bd = BrownDwarf(unum)
    bd.get_low()

    #logging.debug(bd.specs["low"])
    waves.append(bd.specs["low"]["wavelength"])
    fluxes.append(bd.specs["low"]["flux"])


fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(121)
norm_start = 4
spacing = 1.5e-3
for i, name in enumerate(prism_sample["Shortname"]):
    if i==5:
        ax1 = ax
        ax = fig.add_subplot(122)
        norm_start = 9
    norm_by = np.sum(fluxes[i])
    ax.step(waves[i],fluxes[i]/norm_by + (norm_start - i)*spacing,color="k")
    print len(waves[i]), len(fluxes[i])
    if i==7:
        ax.text(1.8,(norm_start - i)*spacing + 1.75*spacing,r"{} {}".format(name,prism_sample["Type"][i]))
    else:
        ax.text(1.8,(norm_start - i)*spacing + 1.5*spacing,r"{} {}".format(name,prism_sample["Type"][i]))
ax2 = ax

ax1.set_ylim(-4e-4,7*spacing)
ax2.set_ylim(-4e-4,7*spacing)
ax1.tick_params(labelleft=False)
ax2.tick_params(labelleft=False)
ax1.set_xlim(0.8,2.5)
ax2.set_xlim(0.8,2.5)
ax1.set_xlabel(r"$\lambda\ $(micron)",fontsize="large")
ax2.set_xlabel(r"$\lambda\ $(micron)",fontsize="large")
ax1.set_ylabel("Flux (normalized)",fontsize="large")
plt.subplots_adjust(wspace=0.0,top=0.95,bottom=0.2,left=0.1,right=0.95)
plt.savefig("Prism_sample.eps",bbox_inches="tight")
