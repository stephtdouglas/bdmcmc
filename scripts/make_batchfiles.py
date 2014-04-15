import bdmcmc.sample.fetch as fetch


ldwarfs = fetch.fetch_12()
bds = ldwarfs.brown_dwarfs
unums = bds.keys()

for unum in unums:
    f = open('run_{}.py'.format(unum),'w')
    f.write('import logging\n')
    f.write('from datetime import date\n\n')
    f.write('from bdmcmc.batch import OneBatch\n\n')
    f.write('logging.basicConfig(level=logging.INFO)\n\n')
    f.write("ob = OneBatch('{}','/vega/astro/users/sd2706/modelSpectra/SpeX_dusty_cut.pkl')\n".format(unum))
    f.close()

    g = open('run_{}.sh'.format(unum),'w')
    g.write("#!/bin/sh\n\n")
    g.write("# Directives\n")
    g.write("#PBS -N Batch{}\n".format(unum))
    g.write("#PBS -W group_list=yetiastro\n")
    g.write("#PBS -l nodes=1,walltime=8:00:00,mem=3000mb\n")
    g.write("#PBS -M sd2706@columbia.edu \n")
    g.write("#PBS -m ae\n")
    g.write("#PBS -V\n\n")
    g.write("# Set output and error directories\n")
    g.write("#PBS -o localhost:/vega/astro/users/sd2706/testing/ \n")
    g.write("#PBS -e localhost:/vega/astro/users/sd2706/testing/ \n")
    g.write("/vega/astro/users/amp2217/yt-x86_64/bin/python run_{}.py\n".format(unum))
    g.write("# End of script\n")
    g.close()
