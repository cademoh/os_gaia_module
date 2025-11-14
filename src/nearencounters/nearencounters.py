# %
# PROCESS THE BAILER-JONES NEAR ENCOUNTERS CATALOG:
# https://cdsarc.cds.unistra.fr/viz-bin/cat/J/ApJ/935/L9
#
#
# ZACK REEVES
# CREATED: 2024
#
# VERSIONS:
#  1.1  MAR 2024 CREATE JUPYTER NOTEBOOK

# %
import pandas as pd
import numpy as np
import sys
import os
import collections

import astropy.units as u
import astropy.coordinates
from astropy.table import Table, join

from astroquery.vizier import Vizier
from astroquery.gaia import Gaia

sys.path.insert(0, '..')
from common import file_functions, calculations

from matplotlib import pyplot as plt, colors

# %
# Define the metadata for the data set. 
metadata = {}

metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
metadata['sub_project'] = 'Close encounters to the Sun'

metadata['catalog'] = 'Close encounters to the Sun in Gaia DR3 (Bailer-Jones, 2022)'
metadata['catalog_author'] = 'Bailer-Jones'
metadata['catalog_year'] = '2022'
metadata['catalog_doi'] = 'doi:10.1051/0004-6361/202039498'  #need to fix
metadata['catalog_bibcode'] = '2021A&A...649A...6G' #need to fix

metadata['prepared_by'] = 'Zack Reeves (AMNH)'
metadata['version'] = '1.1'

metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
metadata['raw_data_dir'] = ''

metadata['data_group_title'] = 'Near Encounters to the Sun'
metadata['data_group_desc'] = 'Near Encounters to the Sun'
metadata['data_group_desc_long'] = '' #need to write

metadata['fileroot'] = 'near_encounters'
file_functions.generate_license_file(metadata)
file_functions.generate_asset_file(metadata)

# %
#reading in the catalogue
catalog = Vizier(catalog='J/ApJ/935/L9/table12', columns=['**'], row_limit=-1).query_constraints()
catalog[0]

# %
#reducing data down to the necessary columns
data = catalog[0][['GaiaDR3', 'tphmed', 'dphmed', 'vphmed', 'Plx', 'e_Plx', 'RV', 'Gmag', 'GMAG', 'GLON', 'GLAT']]

# %
#Query Gaia ESA ADQL server using Gaia EDR3 IDs to obtain proper motion and RA/DEC to calculate uvw

#log in to Gaia Server - Can change to different credentials file for a different user
Gaia.login(credentials_file='../common/gaia_credentials.txt')

#grab username from file
file = open('../common/gaia_credentials.txt', 'r')
username = file.readline().strip()

#Upload table (table name will be forced to lowercase)
job = Gaia.upload_table(upload_resource=data[['GaiaDR3']], table_name="near_encounters", format="csv")

#Query Gaia DR3 source for parallaxes
#Potentially want Bailer Jones distances pending figuring out the parallax error issue
job = Gaia.launch_job_async("select a.GaiaDR3, b.ra, b.dec, b.pmra, b.pmdec, bp_g "
                            "from user_"+username+".near_encounters a left join gaiadr3.gaia_source b on a.GaiaDR3 = b.source_id ",
                            dump_to_file=False)

#put the resulting table into a dataframe and drop the unnecessary index column
data = join(data, job.get_results(), keys='GaiaDR3', join_type='left')
# data.remove_column('xhip_main_oid')
#Deleting table and job from Gaia ESA server so we don't clog the memory
Gaia.delete_user_table('near_encounters')
Gaia.remove_jobs(job.jobid)

Gaia.logout()

# %
data

# %
#calculating distance in light years and parsecs
calculations.get_distance(data, parallax='Plx', use='parallax')

# %
#calculating cartesian coordinates
#calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='RV', frame='icrs')
calculations.get_cartesian(data, glon='GLON', glat='GLAT', pmra='pmra', pmde='pmdec', radial_velocity='RV', frame='icrs')

# %
data

# %
#setting dcalc
#setting metadata for dcalc
#since we calculate distance only using parallax in this dataset, dcalc is always 2
data['dcalc'] = data.Column(data=[2]*len(data),
                            meta=collections.OrderedDict([('ucd', 'meta.dcalc')]),
                            description='Distance Indicator: 1 indicates a Bailer-Jones photogeometric distance; 2 indicates a Bailer-Jones geometric distance')

# %
#calculating absolute magnitudes
#calculate absolute V mag based on apparent magnitude and distance
data['appmag'] = data.MaskedColumn(data=data['Gmag'],
                             unit=u.mag,
                             meta=collections.OrderedDict([('ucd', 'phot.mag;em.opt.G')]),
                             format='{:.6f}',
                             description='Apparent magnitude in Gaia G-band')
data['absmag'] = data.MaskedColumn(data=[data['appmag'][i]+5-5*np.log10(data['dist_pc'][i]) for i in range(len(data))],
                             unit=u.mag,
                             meta=collections.OrderedDict([('ucd', 'phot.magAbs;em.opt.G')]),
                             format='{:.6f}',
                             description='Absolute magnitude in Gaia G-band')

# %
#calculate luminosity based on absolute magnitude
data['lum'] = [10**(1.89 - 0.4*data['absmag'][i]) for i in range(len(data))]
small_luminosities = np.where((data['lum']>0.0) & (data['lum']<0.001))[0]
data['lum'][small_luminosities] = [0.001]*len(small_luminosities)

data['lum'] = data.MaskedColumn(data=data['lum'],
                             unit=u.solLum,
                             meta=collections.OrderedDict([('ucd', 'phys.luminosity')]),
                             format='{:.6f}',
                             description='Stellar Luminosity')

# %
#setting color and visualizing
data['color'] = data.MaskedColumn(data=data['bp_g'],
                             unit=u.solLum,
                             meta=collections.OrderedDict([('ucd', 'phys.color')]),
                             format='{:.2f}',
                             description='Gaia BP-G color')
plt.hist(data['color'], bins=10)

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
data['error_over_parallax']=[data['e_Plx'][i]/data['Plx'][i] for i in range(len(data))]

# %
len(data[data['error_over_parallax']>0.2])

# %
#construct a speck comment column
data['speck_label'] = data.Column(data=['#__'+str(name) for name in data['GaiaDR3']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                  description='Gaia DR3 Source ID')

#construct a label column
data['label'] = ['GaiaDR3_'+ str(source) for source in data['GaiaDR3']]  #leaving for now in case we want to add other labels

# %
#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')

# %
#Getting the column metadata
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'tphmed', 'dphmed', 'vphmed', 'speck_label'])
columns

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
from astropy.units import imperial

# %
min(data['dphmed'].quantity).to(imperial.mile)

# %
data[data['dphmed']<0.07]

# %
data

print("done")
