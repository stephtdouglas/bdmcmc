import bdmcmc.sample.fetch as fetch

model_name = "BT-Settl"
model_file = "SpeX_BTS_wide.pkl"
modelpath = '/vega/astro/users/sd2706/modelSpectra/'
#modelpath = '/home/stephanie/ldwarfs/modelSpectra/'


ldwarfs = fetch.fetch_12()
bds = ldwarfs.brown_dwarfs
unums = bds.keys()

h = open('submit_all_{}.sh'.format(model_name),'w')
for unum in unums:
    h.write('qsub run_{}_{}.sh\n'.format(model_name,unum))

    f = open('run_{}_{}.py'.format(model_name,unum),'w')
    f.write('import logging\n')
    f.write('from datetime import date\n\n')
    f.write('from bdmcmc.batch import OneBatch\n\n')
    f.write('logging.basicConfig(level=logging.INFO)\n\n')
    f.write("ob = OneBatch('{}','/vega/astro/users/sd2706/modelSpectra/'"
            " )\n".format(model_file,unum))
    f.close()

    g = open('run_{}_{}.sh'.format(model_name,unum),'w')
    g.write("#!/bin/sh\n\n")
    g.write("# Directives\n")
    g.write("#PBS -N {}{}\n".format(model_name,unum))
    g.write("#PBS -W group_list=yetiastro\n")
    g.write("#PBS -l nodes=1,walltime=8:00:00,mem=4000mb\n")
    g.write("#PBS -M sd2706@columbia.edu \n")
    g.write("#PBS -m a\n")
    g.write("#PBS -V\n\n")
    g.write("# Set output and error directories\n")
    g.write("#PBS -o localhost:/vega/astro/users/sd2706/outputs/ \n")
    g.write("#PBS -e localhost:/vega/astro/users/sd2706/outputs/ \n")
    g.write("/vega/astro/users/amp2217/yt-x86_64/bin/python "
            " run_{}_{}.py\n".format(model_name,unum))
    g.write("# End of script\n")
    g.close()
h.close()
