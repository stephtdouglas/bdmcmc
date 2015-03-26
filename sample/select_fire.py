import logging

import numpy as np
import matplotlib
matplotlib.use('agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import astropy.units as u

import bdmcmc.spectra, bdmcmc.batch
from bdmcmc.config import *

logging.basicConfig(level=logging.DEBUG)

def fetch_fire():
    """
    Searches the database for all MAGELLAN FIRE spectra
    and returns the object names and other info
    """

    # Select all spex spectra from the database
    all_fire = db.dict.execute("SELECT s.obs_date, s.filename, s.wavelength,"
        " so.shortname, s.source_id FROM spectra AS s JOIN sources AS so ON "
        "so.id=s.source_id WHERE s.instrument_id=11").fetchall()

    mag_fire = [[i['shortname'],i['filename'],i['source_id']] 
                 for i in all_fire]

    return mag_fire

def make_fire_batch(model_name="marley",model_file="SXD_Marley.pkl"):

    logging.debug("model {} file {}".format(model_name,model_file))

    mag_fire = fetch_fire()
    logging.debug("Got Sample")

    h = open('submit_all_{}.sh'.format(model_name),'w')
    for name,filename,sid in mag_fire:
        if name==None:
            name = db.query.execute("SELECT names FROM sources WHERE id={}".format(sid)
                                    ).fetchone()[0]

        logging.debug("name {} obs on {} sid {}".format(name,filename,sid))
        h.write('qsub run_{}_{}.sh\n'.format(model_name,name))

        f = open('run_{}_{}.py'.format(model_name,name),'w')
        f.write('import logging\n')
        f.write('from datetime import date\n\n')
        f.write('from bdmcmc.batch import OneBatch\n\n')
        f.write('logging.basicConfig(level=logging.INFO)\n\n')
        f.write("ob = OneBatch('{}',"
             "'/vega/astro/users/sd2706/modelSpectra/{}',"
             "'{}',filename='{}',spectrograph='{}')\n".format(
             name,model_file,model_name,filename,"FIRE"))
        f.close()

        g = open('run_{}_{}.sh'.format(model_name,name),'w')
        g.write("#!/bin/sh\n\n")
        g.write("# Directives\n")
        g.write("#PBS -N {}{}\n".format(model_name,name))
        g.write("#PBS -W group_list=yetiastro\n")
        g.write("#PBS -l nodes=1,walltime=24:00:00,mem=2500mb\n")
        g.write("#PBS -M sd2706@columbia.edu \n")
        g.write("#PBS -m a\n")
        g.write("#PBS -V\n\n")
        g.write("# Set output and error directories\n")
        g.write("#PBS -o localhost:/vega/astro/users/sd2706/outputs/ \n")
        g.write("#PBS -e localhost:/vega/astro/users/sd2706/outputs/ \n")
        g.write("/vega/astro/users/amp2217/yt-x86_64/bin/python"
             " run_{}_{}.py\n".format(model_name,name))
        g.write("# End of script\n")
        g.close()
    h.close()

def plot_all_fire():

    mag_fire = fetch_fire()

    pp = PdfPages('All_FIRE.pdf')
    plt.figure()
    ax = plt.subplot(111)
    for name,filename,sid in mag_fire:
        if name==None:
            name = db.query.execute("SELECT names FROM sources WHERE id={}".format(sid)
                                    ).fetchone()[0]
        if (name!=None) and (name!='0559-1404') and (name!='1254-0122'):
            bd = bdmcmc.spectra.BrownDwarf(name)
            bd.get_low(filename=filename)
            print name,filename,bd.specs['low']['slit_width']
            ax.step(bd.specs['low']['wavelength'],bd.specs['low']['flux'])
            ax.set_title(name)
            ax.set_xlim(0.8,2.6)
            pp.savefig()
            ax.cla()
    pp.close()


def plot_fire_res():
    mag_fire = fetch_fire()

    fig1 = plt.figure()
    ax1 = plt.subplot(111)
    for name,filename,sid in mag_fire:
        if name==None:
            name = db.query.execute("SELECT names FROM sources WHERE id={}".format(sid)
                                    ).fetchone()[0]
        if (name!=None) and (name!='0559-1404') and (name!='1254-0122'):
            bd = bdmcmc.spectra.BrownDwarf(name)
            bd.get_low(filename=filename)
            print name,filename,bd.specs['low']['slit_width']
            w = bd.specs['low']['wavelength']
            wdiff = w[2:]-w[:-2]
            R = w[:-2] / wdiff * (
                0.3/bd.specs['low']['slit_width'].value)
            Rclean = R[R>1000]
            wclean = w[:-1][R>1000]#w[:-2][wdiff>0.01*u.um]
            ax1.plot(wclean, Rclean,'.',label=name)
#    ax1.legend(ncol=3)

def calc_fire_res():
    mag_fire = fetch_fire()
    all_w = np.array([])
    all_R = np.array([])

    for name,filename,sid in mag_fire:
         if (name!=None) and (name!='0559-1404') and (name!='1254-0122'):
             bd = bdmcmc.spectra.BrownDwarf(name)
             bd.get_low(filename=filename)
             if bd.specs['low']['slit_width']==(0.5*u.um):
                  w = bd.specs['low']['wavelength']
                  wdiff = w[2:]-w[:-2]
                  R = w[:-2] / wdiff * (0.3/0.5)
                  Rclean = R[(R>1000) & (R<1375)]
                  wclean = w[:-1][(R>1000) & (R<1375)]
                  np.append(all_w,wclean)
                  np.append(all_R,Rclean)
         all_w.sort()
         all_R.sort()
         split_points = [0.93891,1.1277,1.41055,1.85]

         # there are 5 orders
         # for each order, fit a line to the R vs. lambda curve
         # Then use that point to calculate the fitted R as a
         # function of all the data wavelength points
         # Smooth the models to the right R, but resample onto this
         # horrendously oversampled wavelength grid, so that the code
         # can just interpolate

def get_source_info():
    mag_fire = fetch_fire()
    # Info I want:
    # sources: ra, dec, publication_id, shortname, components
    # parallaxes: parallax, parallax_unc, publication_id
    # spectral_types: spectral_type, gravity, publication_id, regime, adopted
    # spectra: publication_id
    # (tie all above to publications: id, shortname, eventually bibtex)

    outfile = "/home/stephanie/ldwarfs/FIRE_info.dat"
    f = open(outfile,"w")
    f.write("shortname\tcomponents\t")
    f.write("spec_pub\t")
    f.write("SpT\tGrav\tRegime\n")
    for name,filename,sid in mag_fire:
        if name==None:
            name = db.query.execute("SELECT names FROM sources WHERE id={}".format(sid)
                                    ).fetchone()[0]
        result = db.dict.execute("SELECT s.ra, s.dec, "#s.publication_id, "
            "s.shortname, s.components FROM sources AS s WHERE s.id={}".format(
            sid)).fetchall()
        #print len(result)
        f.write("{}\t{}\t".format(result[0]["shortname"],
            result[0]["components"]))

        result = db.dict.execute("SELECT p.shortname "
            " FROM spectra AS s JOIN publications AS p ON "
            " s.publication_id=p.id WHERE s.source_id={}".format(
            sid)).fetchall()
        if len(result)>0:
            f.write("{}\t".format(result[0]["shortname"]))
        else:
            f.write("unknown\t")

        result = db.dict.execute("SELECT p.shortname, "
            " s.spectral_type, s.gravity, s.regime "
            " FROM spectral_types AS s JOIN publications AS p ON "
            " s.publication_id=p.id WHERE s.source_id={}".format(
            sid)).fetchall()
        optical_found=False
        for i,res_line in enumerate(result):
            if res_line["regime"]=="OPT":
                f.write("{}\t{}\t{}\t".format(res_line["spectral_type"],
                     res_line["gravity"],res_line["regime"]))
                #print "{}\t{}\t{}\t".format(res_line["spectral_type"],
                #     res_line["gravity"],res_line["regime"])
                optical_found=True
                break
        #print "spt {} optical {}".format(len(result),optical_found)
        if (optical_found==False) and (len(result)>0):
            res_line = result[0]
            f.write("{}\t{}\t{}\n".format(res_line["spectral_type"],
                 res_line["gravity"],res_line["regime"]))
        else:
            f.write("n/a\tn/a\tn/a\n")
    f.close()

#        result = db.dict.execute("SELECT "
#            " FROM  WHERE .source_id={}".format(
#            sid))
#        f.write("{}\t{}\t".format())


#get_source_info()

#plot_fire_res()

#plot_all_fire()



#make_fire_batch("Marley","SXD_r2000_Marley.pkl")
make_fire_batch("BTSettl","SXD_r2000_BTS.pkl")
#make_fire_batch("B06_","B06_cloud_exp_r1200.pkl")
#make_fire_batch("Dusty","SXD_r2000_Dusty.pkl")
