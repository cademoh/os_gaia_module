#running the query required to build the comoving stars catalogue

from astroquery.gaia import Gaia
from astroquery.utils.tap.core import TapPlus

import astropy.table as table
from astropy.table import Table



#running the query in the paper with columns described in find_binaries_edr3.py
#to change to DR3 context, change gaiaedr3.gaia_source to gaiadr3.gaia_source
print('../common/gaia_credentials.txt')

#log in to Gaia Server - Can change to different credentials file for a different user
Gaia.login(credentials_file='../common/gaia_credentials.txt')

#get username from credentials file for query
with open('../common/gaia_credentials.txt', 'r') as file:
    username = file.readline()
print("inputting job query, this can take a very long time. go get a coffee or something")
#submit query
job = Gaia.launch_job_async("select source_id, ra, dec, parallax, parallax_error, pmra, pmdec, pmra_error, pmdec_error, phot_g_mean_mag "
                            "from gaiaedr3.gaia_source "
                            "where parallax > 1 "
                            "and parallax_over_error > 5 "
                            "and parallax_error < 2 "
                            "and phot_g_mean_mag is not null")

# storing results of the query in Astropy Table
data = job.get_results()
print("we got the job")
#writing results to a csv
data.write('/edr3_parallax_snr5_goodG.csv')

#Deleting job from Gaia ESA server so we don't clog the memory
Gaia.remove_jobs(job.jobid)

Gaia.logout()


