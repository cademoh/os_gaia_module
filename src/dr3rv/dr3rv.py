# %
# PROCESS THE GAIA DR3 RADIAL VELOCITY CATALOG:
# https://gea.esac.esa.int/archive/
#
#
# ZACK REEVES
# CREATED: 2024
#
# VERSIONS:
#  1.1  MAR 2024 CREATE JUPYTER NOTEBOOK

# %
# **STEPS TO RUN THIS CODE**

# There are two queries you can run to generate this catalog
# The first is a query to grab every star in Gaia with a radial velocity (and a gmag) ~33 million stars
# The second is a query to grab every Gaia star with a radial velocity AND which pass some data quality tests (parallax error) ~29 million stars

# First choose a query, then run it.  
# If you have a decent internet connection and can run the code for hours, the code will execute the query and download the table
# If the code times out, the query results should still be available on your gaia archive account

# Download the data either from the code or from the gaia archive, then run the rest of the processing code.  This can take several hours
# Recommended to download as a vot.gz - this is the smallest file size the gaia server offers
# If any errors occur, consider slicing the data down to the first 1000 rows (data[:1000]) to debug
# Can also add "select TOP 1000" to the query to grab 1000 stars for testing

# %
#This code pulls a catalogue from Gaia DR3 of each star that has a reliable radial velocity
#Based on https://www.aanda.org/articles/aa/full_html/2023/06/aa44220-22/aa44220-22.html#S14,
#We correct radial velocities for stars of high magnitude grvs_mag>11 by Katz et al.
#We also correct stars of high effective temperature 14500>rv_template_teff>8500 and 6>grvs_mag>12 by Blomme et al. as recommended by Katz^

# %
import pandas as pd
import numpy as np
import sys
import os
import collections

import astropy.units as u
import astropy.coordinates
from astropy.table import Table
from astropy.io import ascii
from astropy.io import fits

from astroquery.gaia import Gaia

sys.path.insert(0, '..')
from common import file_functions, calculations, gaia_functions

from matplotlib import pyplot as plt, colors

# %
# Define the metadata for the data set. 
metadata = {}

metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
metadata['sub_project'] = 'Gaia DR3 Radial Velocities'

metadata['catalog'] = 'Gaia Data Release 3: Properties and validation of the radial velocities (Katz et al., 2023)'
metadata['catalog_author'] = 'Katz et al.'
metadata['catalog_year'] = '2023'
metadata['prepared_by'] = 'Zack Reeves (AMNH)'
metadata['version'] = '1.1'

metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
metadata['raw_data_dir'] = ''

metadata['data_group_title'] = 'RadialVelocityStars'
metadata['data_group_desc'] = 'Gaia DR3 Radial Velocity'
metadata['data_group_desc_long'] = 'Gaia DR3 Radial Velocity' #need to expand

metadata['fileroot'] = 'gdr3rv'

file_functions.generate_asset_file(metadata)
file_functions.generate_license_file(metadata)

# %
#query the catalogue from Gaia
#https://gea.esac.esa.int/archive/
#The query pulls the source id, positional data, velocity data as well as teff and magnitude for correction purposes
#corrective data included in the query was informed by https://www.aanda.org/articles/aa/full_html/2023/06/aa44220-22/aa44220-22.html#S14
print("Which query would you like to run? type '1' or '2' for the queries, or '3' for a downloaded file")
querytype = int(input())
# %
if querytype == 1:
    #QUERY #1 - 33 million stars

    #log in to Gaia Server - Can change to different credentials file for a different user
    #query runs in a little over an hour
    #file is 3.2 gigabytes, 33,653,049 objects
    Gaia.login(credentials_file='../common/gaia_credentials.txt')

    #Query Gaia DR3 source for parallaxes
    job = Gaia.launch_job_async("select a.source_id, a.ra, a.dec, a.pmra, a.pmdec, a.parallax, a.parallax_error, a.phot_g_mean_mag, a.bp_g, a.radial_velocity, a.radial_velocity_error, a.grvs_mag, a.rv_template_teff, "
                                "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo "
                                "from gaiadr3.gaia_source a left join external.gaiaedr3_distance bj on a.source_id = bj.source_id "
                                "where a.radial_velocity is not null and a.phot_g_mean_mag > 0 and parallax > 0",
                                dump_to_file=False)

    #Put the resulting table into a Table
    data = job.get_results()
        
    #Gaia.remove_jobs(job.jobid) # UNCOMMENT THIS LINE IF YOU WANT TO PURELY READ THE DATA FROM NOTEBOOK CODE - otherwise remember to delete the job from the gaia archive to not clog your memory

    Gaia.logout()

# %
elif querytype == 2:
    #QUERY #2 - 29 million stars

    #log in to Gaia Server - Can change to different credentials file for a different user
    #query runs in a little over an hour
    #file is 2.8 gigabytes, 29,946,388 objects
    Gaia.login(credentials_file='../common/gaia_credentials.txt')

    #Query Gaia DR3 source for parallaxes
    job = Gaia.launch_job_async("select a.source_id, a.ra, a.dec, a.pmra, a.pmdec, a.parallax, a.parallax_error, a.phot_g_mean_mag, a.bp_g, a.radial_velocity, a.radial_velocity_error, a.grvs_mag, a.rv_template_teff, "
                                "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo "
                                "from gaiadr3.gaia_source a left join external.gaiaedr3_distance bj on a.source_id = bj.source_id "
                                "where a.radial_velocity is not null and a.phot_g_mean_mag > 0 and parallax > 0 and a.parallax / a.parallax_error > 5",
                                dump_to_file=False)

    #Put the resulting table into a Table
    data = job.get_results()
        
    Gaia.remove_jobs(job.jobid)

    Gaia.logout()

elif querytype == 3:
    #uncomment this if you are using downloaded data, keep it commented if you are trying to run from the query
    data = Table.read('dr3rv_query_1_result.fits')

print("data acquired")

# %
data

# %
gaia_functions.set_bj_distance(data)

# %
data

# %
#calculating distance in light years and parsecs
calculations.get_distance(data, dist='bj_distance', use='distance')
print("calculations done")
# %
# data quality check
data.remove_rows(np.where(data['dist_pc']<=0)[0])

# %
gaia_functions.get_magnitudes(data)

# %
gaia_functions.get_luminosity(data)

# %
gaia_functions.get_bp_g_color(data) #may want to change to bp_rp

print("functions done")
# %
# data check on the G mag
x = data['phot_g_mean_mag']
q25, q75 = np.percentile(x, [25, 75])
bin_width = 2 * (q75 - q25) * len(x) ** (-1/3)
bins = round((x.max() - x.min()) / bin_width)
print("Freedman–Diaconis number of bins:", bins)
plt.hist(x, bins=bins);

# %
#applying corrections from papers
data['radial_velocity_correction'] = [0.0]*len(data)

#Katz correction
katz_indexes = np.where(data['grvs_mag']>11)[0]
data['radial_velocity_correction'][katz_indexes] = [(0.02755*data['grvs_mag'][i]**2 - 0.55863*data['grvs_mag'][i] + 2.81129) for i in katz_indexes]

#Blomme correction
blomme_indexes = np.where((data['grvs_mag']>11)&(data['rv_template_teff']>8500)&(data['rv_template_teff']<14500))[0]
data['radial_velocity_correction'][blomme_indexes] = [7.98 - 1.135*data['grvs_mag'][i] for i in blomme_indexes]

data['corrected_radial_velocity'] = np.subtract(data['radial_velocity'], data['radial_velocity_correction'])
data['corrected_radial_velocity'].unit=u.km/u.s

# %
#calculating cartesian coordinates
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='corrected_radial_velocity', frame='icrs')

# %
data

print("all calculations done")
# %
#2D Visualization
fig, ax = plt.subplots(1, 2)

#XY Plane
ax[0].scatter(data['x'], data['y'])
ax[0].set_title('XY Plane')

#XZ Plane
ax[1].scatter(data['x'], data['z'])
ax[1].set_title('XZ Plane')

#set good spacing
fig.tight_layout()
fig.set_size_inches(10, 4, forward=True)
plt.show

# %
#2D Density Visualization
fig, ax = plt.subplots(1, 2)

#XY Plane
ax[0].hist2d(data['x'], data['y'], 
           bins = 200,  
           norm = colors.LogNorm(),  
           cmap = "RdYlGn_r",) 
ax[0].set_title('XY Plane')

#XZ Plane
ax[1].hist2d(data['x'], data['z'], 
           bins = 200,  
           norm = colors.LogNorm(),  
           cmap = "RdYlGn_r",) 
ax[1].set_title('XZ Plane')

#set good spacing
fig.tight_layout()
fig.set_size_inches(10, 4, forward=True)
#plt.show

# %
#construct a speck comment column
data['speck_label'] = data.Column(data=['#__'+str(name) for name in data['SOURCE_ID']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                  description='Gaia DR3 Source ID')

#construct a label column
data['label'] = ['GaiaDR3_'+ str(source) for source in data['SOURCE_ID']]  #leaving for now in case we want to add other labels

# %
#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')

# %
#Getting the column metadata
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'speck_label'])
columns

print("plots done")
# %
# Print the csv file using the to_csv function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(data), columns)

# %
# Print the speck file using the to_speck function in file_functions
file_functions.to_speck(metadata, Table.to_pandas(data), columns)

# %
# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))

# %



