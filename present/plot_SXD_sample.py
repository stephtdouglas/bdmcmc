import numpy as np
import matplotlib.pyplot as plt
import asciitable as at

from bdmcmc.spectra import BrownDwarf

figurepath = "/home/stephanie/Dropbox/paperBD/figures/"

sxd_sample = at.read("SXD_table.csv")

waves = []
fluxes = []

for i, unum in enumerate(sxd_sample["unum"]):
    bd = BrownDwarf(unum)
    bd.get_low(obs_date=sxd_sample["obs_date"][i])

    waves.append(bd.specs["low"]["wavelength"])
    fluxes.append(bd.specs["low"]["flux"])

fig = plt.figure(figsize=(10,8))
ax = fig.add_subplot(121)
norm_start = 7
for i, name in enumerate(sxd_sample["shortname"]):
    if i==8:
        ax1 = ax
        ax = fig.add_subplot(122)
        norm_start = 15
    norm_by = np.sum(fluxes[i])
    ax.step(waves[i],fluxes[i]/norm_by + (norm_start - i)*6e-4,color="k")
    ax.text(1.8,(norm_start - i)*6e-4 + 4.5e-4,r"{} {}".format(name,sxd_sample["SpT2"][i]))
ax2 = ax

ax1.set_ylim(0,8.1*6e-4)
ax2.set_ylim(0,8.1*6e-4)
ax1.tick_params(labelleft=False)
ax2.tick_params(labelleft=False)
ax1.set_xlim(0.8,2.5)
ax2.set_xlim(0.8,2.5)
ax1.set_xlabel(r"$\lambda\ $(micron)",fontsize="large")
ax2.set_xlabel(r"$\lambda\ $(micron)",fontsize="large")
ax1.set_ylabel("Flux (normalized)",fontsize="large")
plt.subplots_adjust(wspace=0.0,top=0.95,bottom=0.2,left=0.1,right=0.95)
plt.savefig("SXD_sample.eps",bbox_inches="tight")
