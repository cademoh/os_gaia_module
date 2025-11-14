# %
# PROCESS THE GLOBULAR CLUSTER STARS CATAOLOG:
# https://zenodo.org/records/4891252
#
#
# ZACK REEVES
# CREATED: 2024
#
# VERSIONS:
#  1.1  JAN 2024 CREATE JUPYTER NOTEBOOK

# %
import pandas as pd
import numpy as np
import sys
import os
import collections

from tqdm import tqdm

from astropy.io import ascii
import astropy.units as u
import astropy.coordinates
from astropy.table import unique, vstack, Table, join

from astroquery.gaia import Gaia

sys.path.insert(0, '..')
from common import file_functions, calculations, gaia_functions, get_bailer_jones

from matplotlib import pyplot as plt, colors

# %
# Define the metadata for the data set. #FIX LATER
metadata = {}

metadata['project'] = 'Digital Universe Atlas Gaia Subsets'
metadata['sub_project'] = 'Gaia Stars in Globular Clusters'

metadata['catalog'] = 'Catalogue of stars in Milky Way globular clusters from Gaia EDR3 (Vasiliev+, 2021)'
metadata['catalog_author'] = 'Vasiliev+'
metadata['catalog_year'] = '2021'
metadata['catalog_doi'] = 'https://doi.org/10.5281/zenodo.4891252'
metadata['catalog_bibcode'] = '10.5281/zenodo.4891252'

metadata['prepared_by'] = 'Brian Abbott, Zack Reeves'
metadata['version'] = '1.1'

metadata['dir'] = metadata['sub_project'].replace(' ', '_').lower()
metadata['raw_data_dir'] = ''

metadata['data_group_title'] = 'Gaia Globular Cluster Stars'
metadata['data_group_desc'] = 'Stars in the Milky Way identified to be members of globular clusters'
metadata['data_group_desc_long'] = ''
metadata['fileroot'] = 'globstars'

file_functions.generate_license_file(metadata)
file_functions.generate_asset_file(metadata)

# %
#download data from https://zenodo.org/records/4891252

# %
#reading in the data

#data are downloaded in a .zip file.  Once extracted, the stars associated with each cluster are stored in folders
#named by the cluster.  We combine each of these folders into one table with an appended column describing the
#cluster

#iterate through files in catalogues directory and stack

directory_str = 'clusters/catalogues/'
directory = os.fsencode(directory_str)
    
tables=[]
for file in tqdm(os.listdir(directory)):
    filename = os.fsdecode(file)
    if filename.endswith(".txt"):
        #reading in the table
        table = ascii.read(directory_str+filename)
        #adding a column with the cluster name
        table['cluster_name'] = table.Column(data=[filename[:len(filename)-4]]*len(table),
                                             meta = collections.OrderedDict([('ucd', 'meta.name.cluster')]),
                                             description='Name of associated Globular Cluster')
        #adding table to array for stacking
        tables.append(table)

# %
#combining tables in tables array
data = vstack(tables)
data

# %
data['source_id'] = [int(i) for i in data['source_id']]

# %
#renaming columns 'x' and 'y' to not be confused with cartesian x and y and adding clarification
data.rename_column('x', 'cluster_x')
data['cluster_x'].description ='X coordinate centered on cluster'
data.rename_column('y', 'cluster_y')
data['cluster_y'].description ='Y coordinate centered on cluster'

#removing proper motion columns because we'll get them from Gaia DR3
data.remove_column('pmra')
data.remove_column('pmdec')

#removing rows with parallax <=0.0 and rows with memberprob < 0.5
data.remove_rows(np.where(data['plx']<=0.0)[0])
data.remove_rows(np.where(data['memberprob']<0.5)[0])

# %
data.write('gc.fits', overwrite=True)

# %
#data = Table.read('gc.fits')

# %
distances = get_bailer_jones.get_bj_distances(data, get_motion=True)

# %
data = join(data, distances, keys='source_id', join_type='left')
data

# %
data['cluster_name'] = table.Column(data=data['cluster_name'],
                                         meta = collections.OrderedDict([('ucd', 'meta.name.cluster')]),
                                         description='Name of associated Globular Cluster')

# %
#fixing parallax units
data['plx'].unit=u.mas

#fixing RA/Dec units
data['ra'].unit=u.deg
data['dec'].unit=u.deg

# %
#calculating distance in light years and parsecs
calculations.get_distance(data, dist='bj_distance', use='distance')

# %
#calculating cartesian coordinates
calculations.get_cartesian(data, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', frame='icrs')

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
gaia_functions.get_magnitudes(data, gmag='g_mag')
gaia_functions.get_luminosity(data)
gaia_functions.get_bp_g_color(data, color='bp_rp')

# %
#construct a speck comment column
data['speck_label'] = data.Column(data=['#__'+str(name) for name in data['source_id']], 
                                  meta=collections.OrderedDict([('ucd', 'meta.id')]),
                                  description='Gaia EDR3 Source ID')

#construct a label column
data['label'] = ['GaiaEDR3_'+ str(source) for source in data['source_id']]  #leaving for now in case we want to add other labels

#setting texture number column
data['texnum'] = data.Column(data=[1]*len(data), 
                                  meta=collections.OrderedDict([('ucd', 'meta.texnum')]),
                                  description='Texture Number')

# %
#Getting the column metadata
columns = file_functions.get_metadata(data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'cluster_name', 'speck_label'])
columns

# %
# Print the csv file using the to_speck function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(data), columns)

# Print the speck file using the to_speck function in file_functions
file_functions.to_speck(metadata, Table.to_pandas(data), columns)

# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(data))

# %
cluster_names = unique(data, keys='cluster_name')['cluster_name']

# %
#for i in cluster_names:
#    print(i)

# %
#the distances given for the stars in globular clusters are inaccurate due to local motions within the clusters
#we disect the dataset back into its constituent clusters to artificially constrain the stars in each cluster to the hypothetical cluster radius

# %
from astroquery.vizier import Vizier
# #reading in the catalogue
#catalog = Vizier(catalog='J/MNRAS/505/5978', columns=['**'], row_limit=-1).query_constraints()
#clusters = catalog[0]
#clusters

# %
# The table read in the commented out code in the previous cell usefully provides the mean parallaxes of each cluster
# The cluster names in that table don't exactly match the ones we have, so the following table is a custom correlation table
#clustable = Table.to_pandas(clusters)
#clustable.to_csv("globclusters.csv", index=False, encoding='utf-8')

clusters = Table.read('globclusters.csv')

# %
def angle_radius(rscale_theta, distance):
    return distance*np.tan(rscale_theta)

# %
import scipy.stats as stats

# %
cluster_dataframes=[]
for i in tqdm(range(len(cluster_names))):
    cluster = cluster_names[i]
    df = data[data['cluster_name']==cluster]
    
    clusters['plx'].unit=u.mas
    
    cluster_parallax = clusters['plx'][clusters['bonus_column']==cluster]
    if(cluster_parallax>0):
        cluster_distance = (clusters['plx'][clusters['bonus_column']==cluster]).to(u.pc, equivalencies=u.parallax())[0]
    else:
        cluster_distance = sum(df['bj_distance']) / len(df)
        
    cluster_radius = angle_radius(clusters['Rscale'][clusters['bonus_column']==cluster], cluster_distance)
    
    #calculating percentiles of the distances of each star in the cluster; we will maintain these percentilesas we narrow the distribution down
    df['dist_percentile'] = [stats.percentileofscore(df['bj_distance'], df['bj_distance'][i], kind='weak') for i in range(len(df))]
    
    #dropping stars which are egregiously far from the cluster (higher than 5 sigma)
    #df.remove_rows(np.where((df['dist_percentile']<(1-0.999999426696856))|(df['dist_percentile']>0.999999426696856))[0])
    
    #mu is the desired distance (the center of the cluster) and sigma will be the adjusted stdev
    mu = cluster_distance
    #we want 95% of our stars within our calculated diameter, so 4*sigma=diameter -> sigma=radius / 2
    sigma = angle_radius(rscale_theta=14.62*u.arcmin, distance=5181.347*u.pc) / 2

    adjusted_distance_distribution = stats.norm(loc=mu, scale=sigma)

    # df['adjusted_distances'] = [adjusted_distance_distribution.ppf(i/100) for i in df['dist_percentile']]
    df['adjusted_distances'] = [adjusted_distance_distribution.ppf(df['dist_percentile'][i]/100) if((df['dist_percentile'][i]<100.0)&(df['dist_percentile'][i]>0.0)) else df['bj_distance'][i] for i in range(len(df))]

    df['adjusted_distances'].unit=u.pc
    
    calculations.get_distance(df, dist='adjusted_distances', use='distance')
    calculations.get_cartesian(df, ra='ra', dec='dec', pmra='pmra', pmde='pmdec', frame='icrs')
    cluster_dataframes.append(df)

# %
adjusted_data = vstack(cluster_dataframes)

# %
#2D Visualization
fig, ax = plt.subplots(1, 2)

#XY Plane
ax[0].scatter(adjusted_data['x'], adjusted_data['y'])
ax[0].set_title('XY Plane')

#XZ Plane
ax[1].scatter(adjusted_data['x'], adjusted_data['z'])
ax[1].set_title('XZ Plane')

#set good spacing
fig.tight_layout()
fig.set_size_inches(10, 4, forward=True)
plt.show

# %
#2D Density Visualization
fig, ax = plt.subplots(1, 2)

#XY Plane
ax[0].hist2d(adjusted_data['x'], adjusted_data['y'], 
           bins = 200,  
           norm = colors.LogNorm(),  
           cmap = "RdYlGn_r",) 
ax[0].set_title('XY Plane')

#XZ Plane
ax[1].hist2d(adjusted_data['x'], adjusted_data['z'], 
           bins = 200,  
           norm = colors.LogNorm(),  
           cmap = "RdYlGn_r",) 
ax[1].set_title('XZ Plane')

#set good spacing
fig.tight_layout()
fig.set_size_inches(10, 4, forward=True)
#plt.show

# %
#Getting the column metadata
columns = file_functions.get_metadata(adjusted_data, columns=['x', 'y', 'z', 'color', 'lum', 'absmag', 'appmag', 'texnum', 'dist_ly', 'dcalc', 'u', 'v', 'w', 'speed', 'cluster_name', 'speck_label'])
columns

# %
adjusted_data

# %
# Print the csv file using the to_speck function in file_functions
file_functions.to_csv(metadata, Table.to_pandas(adjusted_data), columns)

# Print the speck file using the to_speck function in file_functions
file_functions.to_speck(metadata, Table.to_pandas(adjusted_data), columns)

# Print the label file using the to_label function in file_functions
file_functions.to_label(metadata, Table.to_pandas(adjusted_data))

# %
max(adjusted_data['dist_pc'])

# %
adjusted_data[adjusted_data['dist_pc']>300000]

# %
1000/0.003

# %
clusters[clusters['bonus_column']=='NGC_2419']

# %


print("done")
