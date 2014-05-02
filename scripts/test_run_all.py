
import logging
from datetime import date

from bdmcmc.batch import OneBatch

logging.basicConfig(level=logging.DEBUG)

ob = OneBatch('U10668','/vega/astro/users/sd2706/modelSpectra/SpeX_marley.pkl')
