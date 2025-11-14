# %%
# PROCESS THE YOUNG STARS CATAOLOG:
# https://cdsarc.cds.unistra.fr/viz-bin/cat/J/A+A/620/A172
#
#
# ZACK REEVES
# CREATED: 2023
#
# VERSIONS:
#  1.1  OCT 2023 CREATE JUPYTER NOTEBOOK

# %%
import pandas as pd
import numpy as np
import sys
import collections

from astropy.io import ascii
import astropy.units as u
import astropy.coordinates
from astropy.coordinates import Angle
from astropy.table import unique, vstack, Table, join

from astroquery.vizier import Vizier

sys.path.insert(0, '..')
from common import file_functions, calculations, get_bailer_jones, gaia_functions

from matplotlib import pyplot as plt, colors

# %%
# Define the metadata for the data set.  NEED TO EDIT
#https://www.aanda.org/articles/aa/full_html/2023/06/aa43964-22/aa43964-22.html
metadata = {}

metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
metadata['sub_project'] = 'Young Stars'

metadata['catalog'] = '3D mapping of young stars in the solar neighbourhood with Gaia DR2 (Zari+, 2023)'  #need to edit
metadata['catalog_author'] = 'Zari+'
metadata['catalog_year'] = '2023'
metadata['catalog_doi'] = 'doi:10.1051/0004-6361/202039498' #need to fix
metadata['catalog_bibcode'] = '2021A&A...649A...6G' #need to fix

metadata['prepared_by'] = 'Brian Abbott, Zack Reeves'
metadata['version'] = '1.1'

metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
metadata['raw_data_dir'] = ''

metadata['data_group_title'] = 'Young'
metadata['data_group_desc'] = 'Young Stars'
metadata['data_group_desc_long'] = 'Young Stars in the Milky Way mapped by Gaia'
metadata['fileroot'] = 'young_stars'

file_functions.generate_asset_file(metadata)
file_functions.generate_license_file(metadata)

# %%
#Reading in the catalog with Vizier
#We specify the row limit to make sure we get all the stars in the catalog
#We place constraints on the Parallax as a preliminary thresh
#We specify columns = ['**'] to get all of the columns, not just the default ones
catalog = Vizier(catalog='J/A+A/620/A172', columns=['**'], row_limit=-1).query_constraints(Plx='> 0.0')

# %%
#This catalog comes with 4 tables:
# - Pre main sequence (has SIMBAD column)
# - Upper main sequence (has SIMBAD column)
# - Pre main sequence S=2 tangential velocity
# - Pre main sequence S=3 tangential velocity

#We first label each object with the table it came from
catalog[0]['table'] = catalog[0].Column(data=['Pre-main sequence']*len(catalog[0]),
                                        meta = collections.OrderedDict([('ucd', 'meta.table')]),
                                        description='Catalog Table')

catalog[1]['table'] = catalog[1].Column(data=['Upper main sequence']*len(catalog[1]),
                                        meta = collections.OrderedDict([('ucd', 'meta.table')]),
                                        description='Catalog Table')

catalog[2]['table'] = catalog[2].Column(data=['Pre-main sequence S=2']*len(catalog[2]),
                                        meta = collections.OrderedDict([('ucd', 'meta.table')]),
                                        description='Catalog Table')

catalog[3]['table'] = catalog[3].Column(data=['Pre-main sequence S=3']*len(catalog[3]),
                                        meta = collections.OrderedDict([('ucd', 'meta.table')]),
                                        description='Catalog Table')


#We concatenate these tables into one for a full catalog
#Some stars exist in multiple tables and present as duplicate objects
#We remove duplicate objects using the unique function
data = unique(vstack([catalog[0], catalog[1], catalog[2], catalog[3]], 
              metadata_conflicts='silent'), keys='Source', keep='first')

# %%
data

# %%
#querying Gaia for bailer-jones distances
distances = get_bailer_jones.get_bj_distances(data, source_id='Source')
distances

# %%
data = join(data, distances, keys='Source', join_type='inner')

# %%
gaia_functions.set_bj_distance(data)

# %%
# #fixing parallax units (Vizier labels it as a magnitude, probably meant milliarcseconds (mag versus mas))
# data['Plx'].unit=u.mas

# %%
# #thresh on parallax error (cutting on >10% error removes 1614 stars)
# data['parallax_over_error'] = [data['Plx'][i] / data['e_Plx'][i] for i in range(len(data))]
# data.remove_rows(np.where(data['parallax_over_error']<10)[0])

# %%
#calculating distance in light years and parsecs
#calculations.get_distance(data, parallax='Plx', use='parallax')

# %%
calculations.get_distance(data, dist='bj_distance', use='distance')

# %%
len(data)

# %%
#threshing on distance
data.remove_rows(np.where(data['dist_pc']<0.10)[0])

# %%
data

# %%
#calculating cartesian coordinates
calculations.get_cartesian(data, glon='GLON', glat='GLAT', pmglon='pmGLON', pmglat='pmGLAT', 
                           radial_velocity='RV', frame='galactic')

# %%
gaia_functions.get_magnitudes(data, gmag='Gmag')
gaia_functions.get_luminosity(data)
data['bp_rp'] = [data['BPmag'][i] - data['RPmag'][i] for i in range(len(data))]
gaia_functions.get_bp_g_color(data, color='bp_rp')

# %%
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

# %%
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
plt.show

# %%
#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')

# %%
#construct a speck comment column
data['speck_label'] = data.Column(data=['#__'+str(name) for name in data['Source']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                 description='Gaia DR2 Source ID')

#construct a label column
data['label'] = ['GaiaDR2_'+ str(source) for source in data['Source']]  #leaving for now in case we want to add other labels

# %%
#construct a metadata table
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'appmag', 'absmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'speck_label'])
columns

# %%
# Print the speck file using the to_speck function in file_functions
file_functions.to_speck(metadata, Table.to_pandas(data), columns)

# %%
# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))

# %%
# Print the csv file using the to_label function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(data), columns)

# %%



