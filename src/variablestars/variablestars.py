# %
# PROCESS THE GAIA CATALOG OF VARIABLE STARS:
# Gaia has a table for RR Lyrae stars and Cepheid variables
# gaiadr3.vari_rrlyrae & gaiadr3.vari_cepheid
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
from astropy.table import Table, join, vstack
from astropy.io import ascii

from astroquery.gaia import Gaia

sys.path.insert(0, '..')
from common import file_functions, calculations

from matplotlib import pyplot as plt, colors

# %
# Define the metadata for the data set.  NEED TO EDIT
#https://www.aanda.org/articles/aa/full_html/2023/06/aa43964-22/aa43964-22.html
metadata = {}

metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
metadata['sub_project'] = 'Cepheid and RRLyrae Variable Stars'

metadata['catalog'] = 'The Gaia Catalogue of Nearby Stars (Gaia Collaboration, 2021)'  #need to edit
metadata['catalog_author'] = 'Gaia Collaboration'
metadata['catalog_year'] = '2021'
metadata['catalog_doi'] = 'doi:10.1051/0004-6361/202039498'
metadata['catalog_bibcode'] = '2021A&A...649A...6G'

metadata['prepared_by'] = 'Brian Abbott, Zack Reeves'
metadata['version'] = '1.1'

metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
metadata['raw_data_dir'] = ''

metadata['data_group_title'] = 'Variable stars'
metadata['data_group_desc'] = 'Variable stars in the Milky Way mapped by Gaia'
metadata['data_group_desc_long'] = 'RR Lyrae stars are excellent tracers of the oldest stars (ages greater than or equal to 9/10 Gyr) and standard candles for measuring the distance to stellar systems that are mainly composed of an old stellar population. The Gaia Third Data Release (DR3) publishes a catalogue of full sky RR Lyrae stars observed during the initial 34 months of science operations. They were processed through the Specific Object Study (SOS) pipeline, which was developed to validate and characterise Cepheids and RR Lyrae stars (SOS Cep&RRL) observed by Gaia.'
metadata['fileroot'] = 'variable_stars'

file_functions.generate_license_file(metadata)
file_functions.generate_asset_file(metadata)
# %
#query cepheid table

#log in to Gaia Server - Can change to different credentials file for a different user
#query runs in less than a minute
#file is small, 15021 objects
Gaia.login(credentials_file='../common/gaia_credentials.txt')

#Query Gaia DR3 **we can add more params later, maybe get metallicity, periods, etc

job = Gaia.launch_job_async("select a.source_id, "
                            "b.ra, b.dec, b.pmra, b.pmdec, b.parallax, b.parallax_error, b.phot_g_mean_mag, b.bp_g, b.radial_velocity, b.radial_velocity_error, b.grvs_mag, b.rv_template_teff, "
                            "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo "
                            "from gaiadr3.vari_cepheid a inner join gaiadr3.gaia_source b on a.source_id = b.source_id "
                            "left join external.gaiaedr3_distance bj on a.source_id = bj.source_id",
                            dump_to_file=False)

#Put the resulting table into a Table
cepheids = job.get_results()
cepheids['variable_type'] = [1]*len(cepheids)

Gaia.remove_jobs(job.jobid)

Gaia.logout()

# %
#query RR Lyrae table

#log in to Gaia Server - Can change to different credentials file for a different user
#query runs in a few minutes
#file is smallish, 271779 objects
Gaia.login(credentials_file='../common/gaia_credentials.txt')

#Query Gaia DR3 **we can add more params later, maybe get metallicity, periods, etc

job = Gaia.launch_job_async("select a.source_id, "
                            "b.ra, b.dec, b.pmra, b.pmdec, b.parallax, b.parallax_error, b.phot_g_mean_mag, b.bp_g, b.radial_velocity, b.radial_velocity_error, b.grvs_mag, b.rv_template_teff, "
                            "bj.r_med_geo, bj.r_hi_geo, bj.r_lo_geo, bj.r_med_photogeo, bj.r_hi_photogeo, bj.r_lo_photogeo "
                            "from gaiadr3.vari_rrlyrae a inner join gaiadr3.gaia_source b on a.source_id = b.source_id "
                            "left join external.gaiaedr3_distance bj on a.source_id = bj.source_id",
                            dump_to_file=False)

#Put the resulting table into a Table
rrls = job.get_results()
rrls['variable_type'] = [2]*len(rrls)

Gaia.remove_jobs(job.jobid)

Gaia.logout()

# %
rrls

# %
data = vstack([cepheids, rrls])
data['variable_type'] = data.Column(data['variable_type'],
                            meta=collections.OrderedDict([('ucd', 'meta.vari_type')]),
                            description='Type of Variable Star: 1 indicates Cepheid, 2 indicates RR Lyrae')
data

# %
#setting dcalc based on r_med_geo (if>500pc and photogeo exists, we choose photogeo and set dcalc to 1, else geo and dcalc to 2)
data['dcalc'] = [1 if((not(np.ma.is_masked(data['r_med_photogeo'][i])))and(data['r_med_geo'][i]>500)) else 2 for i in range(len(data))]

#setting metadata for dcalc
data['dcalc'] = data.Column(data['dcalc'],
                            meta=collections.OrderedDict([('ucd', 'meta.dcalc')]),
                            description='Distance Indicator: 1 indicates a Bailer-Jones photogeometric distance; 2 indicates a Bailer-Jones geometric distance')

#Choosing distance based on dcalc
data['bj_distance'] = [data['r_med_photogeo'][i] if data['dcalc'][i]==1 else data['r_med_geo'][i] for i in range(len(data))]
data['bj_distance'].unit=u.pc

#Choosing and calculating distance error based on the distance we chose
data['e_bj_dist'] = [((data['r_hi_photogeo'][i]-data['r_lo_photogeo'][i])/2)*u.pc if((not(np.ma.is_masked(data['r_med_photogeo'][i])))and(data['r_med_geo'][i]>500)) else ((data['r_hi_geo'][i]-data['r_lo_geo'][i])/2)*u.pc for i in range(len(data))]

# %
#calculating distance in light years and parsecs
calculations.get_distance(data, dist='bj_distance', use='distance')

# %
#dropping stars with null distance
data.remove_rows([i for i in range(len(data)) if np.ma.is_masked(data['bj_distance'][i])])

# %
#calculating cartesian coordinates
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', radial_velocity='radial_velocity', frame='icrs')

# %
#calculating absolute magnitudes
#calculate absolute V mag based on apparent magnitude and distance
data['appmag'] = data.MaskedColumn(data=data['phot_g_mean_mag'],
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
plt.hist(data['color'], bins=250)

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
#2D Visualization
fig, ax = plt.subplots(1, 2)

#XY Plane
ax[0].scatter(data[data['variable_type']=='RRLyrae']['x'], data[data['variable_type']=='RRLyrae']['y'], color='red')
ax[0].scatter(data[data['variable_type']=='cepheid']['x'], data[data['variable_type']=='cepheid']['y'], color='blue')

ax[0].set_title('XY Plane')

#XZ Plane
ax[1].scatter(data[data['variable_type']=='RRLyrae']['x'], data[data['variable_type']=='RRLyrae']['z'], color='red')
ax[1].scatter(data[data['variable_type']=='cepheid']['x'], data[data['variable_type']=='cepheid']['z'], color='blue')
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
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'variable_type', 'speck_label'])
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
data[data['SOURCE_ID']==4685634433183799680]

# %


print("done")
