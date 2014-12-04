# Module to set up a few things for bdmcmc
# 2 December 2013, Stephanie Douglas
################################################################################

import os

import BDdb

base_path = '/vega/astro/users/sd2706/'
if os.path.exists(base_path)==False:
    base_path = '/home/stephanie/Dropbox/'

#db = BDdb.get_db(base_path+'BDNYCdb/BDNYC.db')
db = BDdb.get_db(base_path+'BDNYC_new.db')

#db = BDdb.get_db(base_path+"paperBD/BDNYC_Yeti.db")
