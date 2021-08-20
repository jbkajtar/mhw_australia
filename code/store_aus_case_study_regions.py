# JBK 2020-12-17
# store area averaged SST timeseries from observational data, for further processing
# updated to use regionmask module
# 1. read SST data
# 2. create masks for index regions
# 3. compute area average and store SST timeseries
# OUTPUT: sst_ts.aus.NOAA_OISST.AVHRR.v2-1_modified.nc

# load required modules
import numpy as np
import xarray as xr
import pandas
import glob
from datetime import date
import time
import cftime
import scipy.io as io
import regionmask

# set paths and filenames
sourcepath = '' # location of observational data
outpath = ''    # output directory

# define bounds of stored region
wa = np.array([[110, -32.5], [110, -22], [116, -22], [116, -32.5]])
tor = np.array([[132, -18], [132, -9], [148, -9], [148, -18]])
gbr = np.array([[142, -25], [142, -10], [155, -10], [155, -25]])
tas = np.array([[147, -46], [147, -37], [157, -37], [157, -46]])
sab = np.array([[132, -35], [132, -31], [143, -31], [143, -40], [138, -40]])

# define headers and index names
idx_head = "Australian MHW case study regions"
idx_name = ["Western Australia", "Torres Strait", "Great Barrier Reef", "Tasman Sea", "South Australian Basin"]
idx_code = ["wa", "tor", "gbr", "tas", "sab"]
reg_str = 'aus'

# create regionmask
index_masks = regionmask.Regions([wa, tor, gbr, tas, sab], names=idx_name, abbrevs=idx_code, name=idx_head)

print('PROCESSING:')
print(sourcepath)

# split sourcepath names to retrieve dataset parameters
filecode = sourcepath.split("/")
filecode = '.'.join(filecode[4:7])

# read netcdf filenames from sourcepath
infiles = sorted(glob.glob(sourcepath + '*.nc'))

# load data
print('Loading data... ')
start = time.time()
ds = xr.open_mfdataset(infiles,combine='nested',concat_dim='time')
end = time.time()
print(end - start)

# create region masks for dataset
mask = index_masks.mask_3D(ds)

# trim time, and store start and end of time array as labels
# start at 1 Jan 1982, end at 31 Dec 2020
ds = ds.where(ds.time >= np.datetime64('1982-01-01T00:00:00'), drop=True)
ds = ds.where(ds.time <= np.datetime64('2020-12-31T23:59:59'), drop=True)

tlab1 = np.datetime_as_string(ds.time.values[0], unit='Y')
tlab2 = np.datetime_as_string(ds.time.values[-1], unit='Y')
date_start = np.datetime_as_string(ds.time.values[0])
date_end = np.datetime_as_string(ds.time.values[-1])

# year span label
tlab = tlab1 + '-' + tlab2
print(tlab)
print('Start: ' + date_start)
print('End: ' + date_end)

# clear all variables except 'sst'
varlist = list(ds.data_vars)  # first get the variable list
varlist.remove('sst')         # remove 'sst' from the list, which will be kept
ds = ds.drop_vars(varlist)

# remove zlev parameter from OBS data
ds = ds.sel(zlev=0, drop=True)

print('SST dimensions: ')
print(ds.sst.dims)

# create weights for area averaging
# NOTE: It is better to use a model’s original grid cell area (e.g. areacella). 
# cos(lat) works reasonably well for regular lat/lon grids. For irregular 
# grids (regional models, ocean models, …) it is not appropriate.
# See: https://regionmask.readthedocs.io/en/stable/notebooks/mask_3D.html
weights = np.cos(np.deg2rad(ds.lat))

# compute SST indices
print('Computing indices... ')
start = time.time()
sst_ts = ds.weighted(mask * weights).mean(dim=("lat", "lon"))
end = time.time()
print(end - start)

# create output filename
outfile = outpath + 'sst_ts.' + reg_str + '.' + filecode + '.nc'

# set global attributes
sst_ts.attrs['source_code'] = "https://github.com/jbkajtar/mhw_australia"
sst_ts.attrs['title'] = "Area-averaged SST for Australian case study regions"
sst_ts.attrs['summary'] = "Data generated for Kajtar et al., 'A catalogue of marine heatwave characteristics and trends for the Australian region', (2021)"
sst_ts.attrs['source_data'] = filecode
sst_ts.attrs['keywords'] = "marine heatwave; extreme event; impact; ocean warming; Australia; observational"
sst_ts.attrs['Conventions'] = "ACDD-1.3"

# save data to file
print('Saving data... ')
start = time.time()
sst_ts.to_netcdf(outfile, encoding={'abbrevs':{'dtype':'S1'}, 'names':{'dtype':'S1'}})
end = time.time()
print(end - start)

