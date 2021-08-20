# JBK 2021-06-18
# process marine heatwave statistics for Australian region
# 1. read SST from observational dataset
# 2. trim to required region
# 3. compute MHW statistics for region
# 4. store required data to netCDF
# OUTPUT: mhw_cats.aus.NOAA_OISST.AVHRR.v2-1_modified.nc

# load required modules
import numpy as np
import xarray as xr
import pandas
import glob
from datetime import date
import time
import cftime
import scipy.io as io

# marine heatwave modules
import marineHeatWaves_jbk as mhw

# set paths and filenames
sourcepath = '' # location of observational data
outpath = ''    # output directory

# define bounds of stored region
reg_lab = 'aus'
reg_bnds = [100, 170, -50, 0]   # 100E-170E, 50S-0

# define climatological baseline period
clim_b = [1983,2012]   # climatology used by Hobday et al. (2018)

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

# trim to required region
da = ds.where((ds.lon > reg_bnds[0]) & (ds.lon < reg_bnds[1])
    & (ds.lat > reg_bnds[2]) & (ds.lat < reg_bnds[3]), drop=True)
 
# trim time, and store start and end of time array as labels
# start at 1 Jan 1982, end at 31 Dec 2020
da = da.where(da.time >= np.datetime64('1982-01-01T00:00:00'), drop=True)
da = da.where(da.time <= np.datetime64('2020-12-31T23:59:59'), drop=True)

tlab1 = np.datetime_as_string(da.time.values[0], unit='Y')
tlab2 = np.datetime_as_string(da.time.values[-1], unit='Y')
date_start = np.datetime_as_string(da.time.values[0])
date_end = np.datetime_as_string(da.time.values[-1])

# year span label
tlab = tlab1 + '-' + tlab2
print(tlab)
print('Start: ' + date_start)
print('End: ' + date_end)
 
# clear all variables except 'sst'
varlist = list(ds.data_vars)  # first get the variable list
varlist.remove('sst')         # remove 'sst' from the list, which will be kept
da = da.drop_vars(varlist)
 
# remove zlev parameter from OBS data
da = da.sel(zlev=0, drop=True)
 
# create output filename
outfile = outpath + 'mhw_cats.' + reg_lab + '.' + filecode
 
# load required data to memory
print('Reading data to memory... ')
start = time.time()
sst = da.sst.values
time_arr = da.time.values
lat = da.lat.values
lon = da.lon.values
end = time.time()
print(end - start)
 
sst.shape   # check dimensions of array [time,lat,lon]
 
# convert time array to ordinal time array
Ly_set = False   # required for mhw.detect
yr1 = pandas.to_datetime(time_arr[0]).year
yr2 = pandas.to_datetime(time_arr[-1]).year
time_o = [0] * len(time_arr)
for i in range(0,len(time_arr)):
  time_o[i] = pandas.Timestamp(time_arr[i]).toordinal() # convert to pandas time

# convert time array (list) to a numpy array (necessary for 'detect')
time_o = np.array(time_o)
 
# initialise variables for MHW statistics
time_yr = list(range(yr1,yr2+1))
cat = list(range(1,5))
i_which = range(sst.shape[2])
j_which = range(sst.shape[1])
mhw_total_events = np.NaN*np.zeros((len(j_which), len(i_which)))
mhw_cats_dpy = np.NaN*np.zeros((len(time_yr),len(cat),len(j_which), len(i_which)))
mhw_count = np.NaN*np.zeros((len(time_yr),len(j_which), len(i_which)))
mhw_intensity = np.NaN*np.zeros((len(time_yr),len(j_which), len(i_which)))
mhw_duration = np.NaN*np.zeros((len(time_yr),len(j_which), len(i_which)))
mhw_count_tr = np.NaN*np.zeros((len(j_which), len(i_which)))
mhw_intensity_tr = np.NaN*np.zeros((len(j_which), len(i_which)))
mhw_duration_tr = np.NaN*np.zeros((len(j_which), len(i_which)))
ev_max_max = np.NaN*np.zeros((len(j_which), len(i_which)))
ev_max_dur = np.NaN*np.zeros((len(j_which), len(i_which)))
ev_max_start = np.empty((len(j_which), len(i_which)), dtype="<U10")
ev_dur_max = np.NaN*np.zeros((len(j_which), len(i_which)))
ev_dur_dur = np.NaN*np.zeros((len(j_which), len(i_which)))
ev_dur_start = np.empty((len(j_which), len(i_which)), dtype="<U10")

# loop over every gridpoint, and compute MHW statistics
for i in i_which:
 start = time.time()
 print(i, 'of', len(i_which)-1)
 for j in j_which:
  # process single SST timeseries
  sst1 = sst[:,j,i]
  # skip cells with land, ice
  if np.logical_not(np.isfinite(sst1.sum())) or (sst1<-1).sum()>0:
   mhw_total_events[j,i] = 0
  else:
   # detect MHWs
   mhws, clim = mhw.detect(time_o, sst1, climatologyPeriod=clim_b, pctile=90,Ly=Ly_set)
   # perform annual averaging of statistics
   mhwBlock = mhw.blockAverage(time_o, mhws, clim, temp=sst1)
   # store total MHW counts
   mhw_total_events[j,i] = mhwBlock['count'].sum()
   # store days per year (dpy) in each MHW category
   mhw_cats_dpy[:,0,j,i] = mhwBlock['moderate_days']
   mhw_cats_dpy[:,1,j,i] = mhwBlock['strong_days']
   mhw_cats_dpy[:,2,j,i] = mhwBlock['severe_days']
   mhw_cats_dpy[:,3,j,i] = mhwBlock['extreme_days']
   # store store additional MHW variables
   mhw_count[:,j,i] = mhwBlock['count']
   mhw_intensity[:,j,i] = mhwBlock['intensity_max']  # annual mean of max MHW intensity
   mhw_duration[:,j,i] = mhwBlock['duration']
   # computes means and trends
   mean, trend, dtrend = mhw.meanTrend(mhwBlock)
   # store trend data
   mhw_count_tr[j,i] = trend['count']
   mhw_intensity_tr[j,i] = trend['intensity_max']
   mhw_duration_tr[j,i] = trend['duration']
   # store start dates of strongest/longest events
   ev_m = np.argmax(mhws['intensity_max'])         # find strongest (intensity_max) event
   ev_d = np.argmax(mhws['duration'])              # find longest (duration) event
   # store statistics for strongest (intensity_max) event
   ev_max_max[j,i] = mhws['intensity_max'][ev_m]
   ev_max_dur[j,i] = mhws['duration'][ev_m]
   ev_max_start[j,i] = mhws['date_start'][ev_m].strftime("%Y-%m-%d")
   # store statistics for longest (duration) event
   ev_dur_max[j,i] = mhws['intensity_max'][ev_d]
   ev_dur_dur[j,i] = mhws['duration'][ev_d]
   ev_dur_start[j,i] = mhws['date_start'][ev_d].strftime("%Y-%m-%d")
 
 end = time.time()
 print(end - start)

# create xarray dataset for processed results
ds_out = da.drop_vars({'sst'})   # follow format of input dataset, without sst
dim = da.sst.dims                # read dimension names
ds_out.attrs = {}                # clear attributes
ds_out['time'] = time_yr         # set new time coordinate
ds_out = ds_out.assign_coords({"cat" : cat})    # set category coordinate
 
# store new variables in dataset
ds_out['mhw_total_events'] = ((dim[-2], dim[-1]), mhw_total_events)
ds_out['mhw_cats_dpy'] = (('time', 'cat', dim[-2], dim[-1]), mhw_cats_dpy)
ds_out['mhw_count'] = (('time', dim[-2], dim[-1]), mhw_count)
ds_out['mhw_intensity'] = (('time', dim[-2], dim[-1]), mhw_intensity)
ds_out['mhw_duration'] = (('time', dim[-2], dim[-1]), mhw_duration)
ds_out['mhw_count_tr'] = ((dim[-2], dim[-1]), mhw_count_tr)
ds_out['mhw_intensity_tr'] = ((dim[-2], dim[-1]), mhw_intensity_tr)
ds_out['mhw_duration_tr'] = ((dim[-2], dim[-1]), mhw_duration_tr)

# store statistics for strongest and longest events
ds_out['ev_max_max'] = ((dim[-2], dim[-1]), ev_max_max)
ds_out['ev_max_dur'] = ((dim[-2], dim[-1]), ev_max_dur)
ds_out['ev_max_start'] = ((dim[-2], dim[-1]), ev_max_start)
ds_out['ev_dur_max'] = ((dim[-2], dim[-1]), ev_dur_max)
ds_out['ev_dur_dur'] = ((dim[-2], dim[-1]), ev_dur_dur)
ds_out['ev_dur_start'] = ((dim[-2], dim[-1]), ev_dur_start)
 
# set coverage_content_type: https://wiki.esipfed.org/Concepts_Glossary
cct = "physicalMeasurement"

# set attributes of variables
ds_out['time'] = ds_out.time.assign_attrs(units="years", standard_name="time", long_name="calendar year", axis = "T", calendar="proleptic_gregorian")
ds_out['cat'] = ds_out.cat.assign_attrs(units="1", long_name="Marine heatwave category following Hobday et al. (2018) definition", mhw_categories="1: Moderate, 2: Strong, 3: Severe, 4: Extreme")
ds_out['mhw_total_events'] = ds_out.mhw_total_events.assign_attrs(units="1", standard_name="n/a", long_name="Total number of marine heatwaves detected", coverage_content_type="auxiliaryInformation")
ds_out['mhw_cats_dpy'] = ds_out.mhw_cats_dpy.assign_attrs(units="1", standard_name="n/a", long_name="Count of days per year (dpy) in each marine heatwave category", coverage_content_type="auxiliaryInformation")
ds_out['mhw_count'] = ds_out.mhw_count.assign_attrs(units="1", standard_name="n/a", long_name="Count of marine heatwave events in each year", coverage_content_type="auxiliaryInformation")
ds_out['mhw_intensity'] = ds_out.mhw_intensity.assign_attrs(units=da.sst.attrs['units'], standard_name="n/a", long_name="Annual mean of maximum marine heatwave intensities in each year (as an anomaly w.r.t. seasonal climatology)", coverage_content_type="auxiliaryInformation")
ds_out['mhw_duration'] = ds_out.mhw_duration.assign_attrs(units="1", standard_name="n/a", long_name="Mean duration (in days) of marine heatwave events in each year", coverage_content_type="auxiliaryInformation")
ds_out['mhw_count_tr'] = ds_out.mhw_count_tr.assign_attrs(units="1", standard_name="n/a", long_name="Trend in annual counts of marine heatwave events (events/year)", coverage_content_type="auxiliaryInformation")
ds_out['mhw_intensity_tr'] = ds_out.mhw_intensity_tr.assign_attrs(units=da.sst.attrs['units'], standard_name="n/a", long_name="Trend in annual mean of maximum of marine heatwave intensities (degree_C/year)", coverage_content_type="auxiliaryInformation")
ds_out['mhw_duration_tr'] = ds_out.mhw_duration_tr.assign_attrs(units="1", standard_name="n/a", long_name="Trend in annual mean duration of marine heatwaves (days/year)", coverage_content_type="auxiliaryInformation")

# set attributes for strongest and longest event variables
ds_out['ev_max_max'] = ds_out.ev_max_max.assign_attrs(units=da.sst.attrs['units'], standard_name="n/a", long_name="Maximum intensity (as an anomaly w.r.t. seasonal climatology) of grid-point largest maximum intensity marine heatwave", coverage_content_type="auxiliaryInformation")
ds_out['ev_max_dur'] = ds_out.ev_max_dur.assign_attrs(units="1", standard_name="n/a", long_name="Duration (in days) of grid-point largest maximum intensity marine heatwave", coverage_content_type="auxiliaryInformation")
ds_out['ev_max_start'] = ds_out.ev_max_start.assign_attrs(units="date", standard_name="n/a", long_name="Start date of grid-point largest maximum intensity marine heatwave", coverage_content_type="auxiliaryInformation")
ds_out['ev_dur_max'] = ds_out.ev_dur_max.assign_attrs(units=da.sst.attrs['units'], standard_name="n/a", long_name="Maximum intensity (as an anomaly w.r.t. seasonal climatology) of grid-point longest duration marine heatwave", coverage_content_type="auxiliaryInformation")
ds_out['ev_dur_dur'] = ds_out.ev_dur_dur.assign_attrs(units="1", standard_name="n/a", long_name="Duration (in days) of grid-point longest duration marine heatwave", coverage_content_type="auxiliaryInformation")
ds_out['ev_dur_start'] = ds_out.ev_dur_start.assign_attrs(units="date", standard_name="n/a", long_name="Start date of grid-point longest duration marine heatwave", coverage_content_type="auxiliaryInformation")

# set global attributes
ds_out.attrs['source_code'] = "https://github.com/jbkajtar/mhw_australia"
ds_out.attrs['title'] = "Marine heatwave statistics for the Australian region (100E-170E, 50S-0)"
ds_out.attrs['summary'] = "Data generated for Kajtar et al., 'A catalogue of marine heatwave characteristics and trends for the Australian region', (2021)"
ds_out.attrs['source_data'] = filecode
ds_out.attrs['keywords'] = "marine heatwave; extreme event; impact; ocean warming; Australia; observational"
ds_out.attrs['reference climatology'] = str(clim_b[0]) + '-' + str(clim_b[1])
ds_out.attrs['Conventions'] = "ACDD-1.3"
 
# save dataset to netcdf file
print('Saving data... ')
start = time.time()
ds_out.to_netcdf(outfile + '_uncompressed' + '.nc')
end = time.time()
print(end - start)

# save compressed dataset to netcdf file
comp = dict(zlib=True, complevel=5)
encoding = {var: comp for var in ds_out.data_vars}
print('Saving data... ')
start = time.time()
ds_out.to_netcdf(outfile + '.nc', encoding=encoding)
end = time.time()
print(end - start)
