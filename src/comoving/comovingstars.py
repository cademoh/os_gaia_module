# 
# PROCESS THE COMOVING STARS CATAOLOG:
# https://zenodo.org/records/4435257
# https://academic.oup.com/mnras/article/506/2/2269/6131876#
#
# ZACK REEVES
# CREATED: 2024
# CADE MOHRHARDT
# UPDATED: 2025
# VERSIONS:
#  1.1  JAN 2024 CREATE JUPYTER NOTEBOOK
#  Python 3.12.12 OCT 2025

# 
# **STEPS TO RUN THIS CODE**

# There are two ways you can generate this catalog:

# The first (more straightforward) is to go to the above zenodo link to download the catalog
# Download 'all_columns_catalog.fits.gz' - the csv.gz seems to be corrupt

# The second way can be useful for updates with future data releases, but there are some kinks
# Run el_badry_query.py, then run num_neighbors_edr3.py, then run find_binaries_edr3.py, which should generate the catalog

# Once the catalog is obtained, run the rest of the processing code.  There should be ~1.8 million binary systems or >3 million stars, so it can be slow

# If any errors occur, consider slicing the data down to the first 1000 rows (data[:1000]) to debug
# Can also add "select TOP 1000" to the query to grab 1000 stars for testing
# 
import pandas as pd
import numpy as np
import sys
import collections

from matplotlib import pyplot as plt, colors

from astroquery.gaia import Gaia
from astroquery.utils.tap.core import TapPlus

from astropy.io import fits
import astropy.table as table
from astropy.table import Table, vstack

sys.path.insert(0, '..')
from common import file_functions, calculations, gaia_functions
# 
# Define the metadata for the data set. 
metadata = {}

metadata['project'] = 'Digital Universe Atlas'
metadata['sub_project'] = 'Comoving Stars'

metadata['catalog'] = 'A million binaries from Gaia eDR3: sample selection and validation of Gaia parallax uncertainties (El-Badry+, 2021)'
metadata['catalog_author'] = 'El-Badry+'
metadata['catalog_year'] = '2021'
metadata['catalog_doi'] = 'https://doi.org/10.1093/mnras/stab323'
metadata['catalog_bibcode'] = '2021MNRAS.506.2269E'


metadata['prepared_by'] = 'Zack Reeves (AMNH), Cade Mohrhardt (AMNH)'
metadata['version'] = '1.1'

metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
metadata['raw_data_dir'] = ''

metadata['data_group_title'] = 'ComovingStars'
metadata['data_group_desc'] = 'Comoving Star catalog'
metadata['data_group_desc_long'] = 'Comoving Star catalog'
metadata['fileroot'] = 'comov'

file_functions.generate_license_file(metadata)
#file_functions.generate_asset_file(metadata, data_display_type="GaiaRenderable")
file_functions.generate_asset_file(metadata)

# 
# The following block is method 2 for generating the catalog - uncomment if you want to do it that way

# #running code with -i allows us to run the .py files in the same namespace as our .ipynb
#%run -i el_badry_query.py
#from comoving import el_badry_query
#running the query required to build the comoving stars catalogue

# #running code with -i allows us to run the .py files in the same namespace as our .ipynb
# %run -i num_neighbors_edr3.py

# %run -i find_binaries_edr3.py
#from comoving import find_binaries_edr3

# Method 1 - download 'all_columns_catalog.fits.gz' from https://zenodo.org/records/4435257
# Make sure you grab the FITS.gz!!

#Uncomment these lines if you want to see all of the columns in the table
#binaries_table = Table.read('raw_data/all_columns_catalog.fits.gz')
#binaries = binaries_table[['source_id1', 'source_id2', 'ra1', 'ra2', 'dec1', 'dec2', 'parallax1', 'parallax2', 'parallax_error1', 'parallax_error2', 'pmra1', 'pmra2', 'pmdec1', 'pmdec2', 'dr2_radial_velocity1', 'dr2_radial_velocity2', 'phot_g_mean_mag1', 'phot_g_mean_mag2', 'bp_rp1', 'bp_rp2']]

#pares the table down to just the columns we want
#this line currently works with a downloaded .fits.gz file, it may need to be adjusted once the query is working
binaries = Table.read('all_columns_catalog.fits.gz')[['source_id1', 'source_id2', 'ra1', 'ra2', 'dec1', 'dec2', 'parallax1', 'parallax2', 'parallax_error1', 'parallax_error2', 'pmra1', 'pmra2', 'pmdec1', 'pmdec2', 'dr2_radial_velocity1', 'dr2_radial_velocity2', 'phot_g_mean_mag1', 'phot_g_mean_mag2', 'bp_rp1', 'bp_rp2']]
#This line may work when pulling from the el_badry code
#binaries = Table.read('binary_catalog.fits')[['source_id1', 'source_id2', 'ra1', 'ra2', 'dec1', 'dec2', 'parallax1', 'parallax2', 'parallax_error1', 'parallax_error2', 'pmra1', 'pmra2', 'pmdec1', 'pmdec2', 'dr2_radial_velocity1', 'dr2_radial_velocity2', 'phot_g_mean_mag1', 'phot_g_mean_mag2', 'bp_rp1', 'bp_rp2']]

#only uncomment this for testing and debugging, we need the whole dataset 
#binaries = binaries[:1000]

binaries
print("houston we have the data")
# 
#creating a table for the primary stars in the system
data_1 = binaries[['source_id1', 'ra1', 'dec1', 'parallax1', 'parallax_error1', 'pmra1', 'pmdec1', 'dr2_radial_velocity1', 'phot_g_mean_mag1', 'bp_rp1']]

# 
#creating a table for the secondary stars in the system
data_2 = binaries[['source_id2', 'ra2', 'dec2', 'parallax2', 'parallax_error2', 'pmra2', 'pmdec2', 'dr2_radial_velocity2', 'phot_g_mean_mag2', 'bp_rp2']]

# 
#renaming columns to remove the 1s and 2s to facilitate table stacking
data_1.rename_columns(['source_id1', 'ra1', 'dec1', 'parallax1', 'parallax_error1', 'pmra1', 'pmdec1', 'dr2_radial_velocity1', 'phot_g_mean_mag1', 'bp_rp1'], ['source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'dr2_radial_velocity', 'phot_g_mean_mag', 'bp_rp']),
data_2.rename_columns(['source_id2', 'ra2', 'dec2', 'parallax2', 'parallax_error2', 'pmra2', 'pmdec2', 'dr2_radial_velocity2', 'phot_g_mean_mag2', 'bp_rp2'], ['source_id', 'ra', 'dec', 'parallax', 'parallax_error', 'pmra', 'pmdec', 'dr2_radial_velocity', 'phot_g_mean_mag', 'bp_rp'])

# 
#adding a column to indicate which stars are primary or secondary
data_1['binary_index']=[str(i)+'_1' for i in range(len(binaries))]
data_2['binary_index']=[str(i)+'_2' for i in range(len(binaries))]

# 
#stacking the tables together for one whole dataset
data = vstack([data_1, data_2])

# 
data
print("about to start calculating. buckle up")
# 
#calculating distance in light years and parsecs
#this dataset only uses gaia parallaxes to calculate distance to avoid the cpmutational expense of uploading >3 million stars to grab BJ distances
calculations.get_distance(data, parallax='parallax', use='parallax')

#setting metadata for dcalc
data['dcalc'] = data.Column([3]*len(data),
                            meta=collections.OrderedDict([('ucd', 'meta.dcalc')]),
                            description='Distance Indicator: 1 indicates a Bailer-Jones photogeometric distance; 2 indicates a Bailer-Jones geometric distance; 3 indicates a Gaia parallax-based distance')

# 
#calculating cartesian coordinates
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='dr2_radial_velocity', frame='icrs')
print("about to start plotting")
# 
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

# 
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
print("last step, making some files now")
# 
# data quality test
data['error_over_parallax']=[data['parallax_error'][i]/data['parallax'][i] for i in range(len(data))]
len(data[data['error_over_parallax']>0.15])

# 
gaia_functions.get_magnitudes(data)
gaia_functions.get_luminosity(data)
gaia_functions.get_bp_g_color(data, color='bp_rp')

# 
#construct a speck comment column
data['speck_label'] = data.Column(data=['#__'+str(name) for name in data['source_id']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                  description='Gaia EDR3 Source ID')

#construct a label column
data['label'] = ['GaiaEDR3_'+ str(source) for source in data['source_id']]  #leaving for now in case we want to add other labels

# 
#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')

# 
#Getting the column metadata
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'speck_label'])
columns

# 
# Print the speck file using the to_speck function in file_functions
file_functions.to_speck(metadata, Table.to_pandas(data), columns)

# 
# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))

# 
# Print the csv file using the to_csv function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(data), columns)


print("done")