# %
# PROCESS THE GAIA CATALOG OF RGB Stars:
# https://ui.adsabs.harvard.edu/abs/2023ApJS..267....8A/abstract
# https://zenodo.org/records/7945154
#
#
# ZACK REEVES
# CREATED: 2024
#
# VERSIONS:
#  1.1  JUN 2024 CREATE JUPYTER NOTEBOOK

# %
import pandas as pd
import numpy as np
import sys
import os
import collections

import astropy.units as u
import astropy.coordinates
from astropy.table import Table, join, vstack
from astropy.io import ascii

from astroquery.gaia import Gaia

sys.path.insert(0, '..')
from common import file_functions, calculations, gaia_functions

from matplotlib import pyplot as plt, colors

# %
# Define the metadata for the data set.  NEED TO EDIT
#https://ui.adsabs.harvard.edu/abs/2023A%26A...674A..39G/abstract
metadata = {}

metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
metadata['sub_project'] = 'Red Giant Branch Stars'

metadata['catalog'] = 'Robust Data-driven Metallicities for 175 Million Stars from Gaia XP Spectra (Andrae, 2023)'
metadata['catalog_author'] = 'Andrae+'
metadata['catalog_year'] = '2023'
metadata['catalog_doi'] = 'doi:10.3847/1538-4365/acd53e'
metadata['catalog_bibcode'] = '2023ApJS..267....8A'

metadata['prepared_by'] = 'Brian Abbott, Zack Reeves'
metadata['version'] = '1.1'

metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
metadata['raw_data_dir'] = ''

metadata['data_group_title'] = 'Giants'
metadata['data_group_desc'] = 'Red Giant Branch Stars' #need to fix
metadata['data_group_desc_long'] = 'The Sun is the reference point in much of stellar astronomy and astrophysics. Solar analogues are stars that resemble the Sun in terms of a restricted set of parameters. In contrast to the Sun, they can be observed in the night sky and with the very same instruments used to study stars in the Milky Way.'
metadata['fileroot'] = 'giant'

file_functions.generate_license_file(metadata)
file_functions.generate_asset_file(metadata)
print("asset file created")
# %
#download the data from https://zenodo.org/records/7945154 
#~12 million stars
data = Table.read('table_2_catwise.fits.gz')
data
print("data is read")
# %
#calculating distance in light years and parsecs
#this dataset only uses gaia parallaxes to calculate distance to avoid the cpmutational expense of uploading >3 million stars to grab BJ distances

data['parallax'].unit=u.mas
calculations.get_distance(data, parallax='parallax')

#setting metadata for dcalc
data['dcalc'] = data.Column([3]*len(data),
                            meta=collections.OrderedDict([('ucd', 'meta.dcalc')]),
                            description='Distance Indicator: 1 indicates a Bailer-Jones photogeometric distance; 2 indicates a Bailer-Jones geometric distance; 3 indicates a Gaia parallax-based distance')

print("partway through calculations (7 min)")
# %
#setting necessary units and calculating galactic cartesian XYZ
data['ra'].unit=u.deg
data['dec'].unit=u.deg
data['pmra'].unit=u.mas/u.yr
data['pmdec'].unit=u.mas/u.yr
data['radial_velocity'].unit=u.km/u.s
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='radial_velocity', frame='icrs')

print("halfway through calulations (7 min)")
# %
#setting necessary units
data['phot_g_mean_mag'].unit=u.mag
data['phot_bp_mean_mag'].unit=u.mag
data['phot_rp_mean_mag'].unit=u.mag

# %
#calculating absolute and apparent magnitudes, luminosity, and color
gaia_functions.get_magnitudes(data)
gaia_functions.get_luminosity(data)
data['bp_rp'] = [data['phot_bp_mean_mag'][i]-data['phot_rp_mean_mag'][i] for i in range(len(data))]
gaia_functions.get_bp_g_color(data, color='bp_rp')

print("calculations done (4 min)")
# %
data

# %
plt.hist(data['bp_rp'], bins=250);

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
data['speck_label'] = data.Column(data=['#__'+str(name) for name in data['source_id']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                  description='Gaia EDR3 Source ID')

#construct a label column
data['label'] = ['GaiaEDR3_'+ str(source) for source in data['source_id']]  #leaving for now in case we want to add other labels

# %
#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')

# %
#Getting the column metadata
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'speck_label'])
columns

print("plots done (6 min)")
# %
# Print the csv file using the to_csv function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(data), columns)

# %
# Print the speck file using the to_speck function in file_functions
file_functions.to_speck(metadata, Table.to_pandas(data), columns)

# %
# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))


print("done")