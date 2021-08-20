# JBK 2021-05-01
# process marine heatwave statistics for Australian case study regions
# 1. read observational dataset of area-averaged SST timeseries
# 2. compute MHW statistics for selected region
# 3. create plot

# load required modules
import numpy as np
import numpy.matlib
import scipy.io
import pandas as pd
import xarray as xr
from datetime import date, timedelta

# MHW modules
import marineHeatWaves_jbk as mhw

# plotting
import matplotlib.pyplot as plt
import matplotlib.dates as dt
import matplotlib.ticker as ticker
import matplotlib.patches as mpatches

inpath = ''   # location of input data
outpath = ''  # output path of figure

# set desired region, options: tas, wa, sab, tor, gbr
reg_code = 'tor'

# define climatological baseline period
clim_b = [1983,2012]   # climatology used by Hobday et al. (2018)

# read ONI
fd = open(inpath + 'oni.data.txt','r')
d = np.loadtxt(fd,skiprows=4,max_rows=71,dtype={'names':('oni_yr','m1','m2','m3','m4',\
    'm5','m6','m7','m8','m9','m10','m11','m12'), 'formats': ('i4','f4','f4','f4','f4',\
    'f4','f4','f4','f4','f4','f4','f4','f4')})
fd.close()

# reshape monthly ONI date into single 1D array
oni_yr=d['oni_yr']
oni=np.stack((d['m1'],d['m2'],d['m3'],d['m4'],d['m5'],d['m6'],d['m7'],d['m8'],d['m9'],d['m10'],d['m11'],d['m12']), axis=1)
oni=oni.flatten()

# create monthly time array for ONI
oni_yr=np.matlib.repmat(oni_yr,12,1)
oni_mths=np.matlib.repmat(1/12*np.arange(0.5,12.5,1),1,oni_yr.shape[1])
oni_yr=np.reshape(oni_yr,oni_mths.shape,order='F')
oni_time=oni_yr+oni_mths
oni_time=oni_time.T
oni_time=oni_time[:,0]

# create 2D ONI array for shading
oni_d=np.arange(-80,1,1)
oni=np.matlib.repmat(oni,len(oni_d),1).T

# open observed SST timeseries data
infile = inpath + 'sst_ts.aus.NOAA_OISST.AVHRR.v2-1_modified.nc'
ds = xr.open_dataset(infile)

# read data
time = ds.time.values
ds = ds.where(ds.abbrevs==reg_code,drop=True)
sst = ds.sst.values[:,0]
reg_name = ds.names.values[0]

# convert time to numpy datetime64 array
time = time.astype('datetime64',copy=False)

# convert time array to ordinal time array
# see http://stackoverflow.com/a/13753918/4279
Ly_set = False   # required for mhw.detect, set to 360 for Datetime360Day
if np.issubdtype(time.dtype, np.datetime64):
    # handle time arrays that are numpy.datetime64
    yr1 = pd.to_datetime(time[0]).year
    yr2 = pd.to_datetime(time[-1]).year
    time_o = [0] * len(time)
    for i in range(0,len(time)):
        time_o[i] = pd.Timestamp(time[i]).toordinal() # convert to pandas time

# convert time array (list) to a numpy array (necessary for 'detect')
time_o = np.array(time_o)

# redefine some variables for use with Eric Oliver's code
t = time_o
dates = [date.fromordinal(tt.astype(int)) for tt in t]
time_yr = list(range(yr1,yr2+1))

# detect MHWs
mhws, clim = mhw.detect(time_o, sst, climatologyPeriod=clim_b, pctile=90,Ly=Ly_set)
mhwBlock = mhw.blockAverage(time_o, mhws, clim, temp=sst)
rank, returnPeriod = mhw.rank(time_o, mhws)

# set cat 1, 2, 3, and 4 levels
cat_diff = clim['thresh'] - clim['seas']
thresh_cat1 = clim['seas'] + cat_diff
thresh_cat2 = clim['seas'] + 2*cat_diff
thresh_cat3 = clim['seas'] + 3*cat_diff
thresh_cat4 = clim['seas'] + 4*cat_diff

# create table of strongest and longest events
for i in range(mhws['n_events']):
 if rank['duration'][i] < 4 or rank['intensity_max'][i] < 4 \
 or rank['intensity_max_norm'][i] < 4 or rank['intensity_cumulative'][i] < 4:
  print(mhws['date_start'][i].strftime("%d/%m/%Y"), '-', mhws['date_end'][i].strftime("%d/%m/%Y") \
  + ',', mhws['duration'][i], '(' + str(rank['duration'][i]) + '),', \
  "{:.2f}".format(mhws['intensity_max'][i]), '(' + str(rank['intensity_max'][i]) + '),', \
  "{:.2f}".format(mhws['intensity_max_norm'][i]), '(' + str(rank['intensity_max_norm'][i]) + '),', \
  "{:.1f}".format(mhws['intensity_cumulative'][i]), '(' + str(rank['intensity_cumulative'][i]) + '),', i)

# create event list for plotting
if reg_code == 'wa':
    ev_list=[27, 39, 45]
elif reg_code == 'tas':
    ev_list=[76, 79, 81]
elif reg_code == 'tor':
    ev_list=[69, 84, 90]
elif reg_code == 'gbr':
    ev_list=[14, 51, 67, 80]
elif reg_code == 'sab':
    ev_list=[40, 54, 70]

# list of years of target events
ev_yrs=[mhws['date_peak'][x].year for x in ev_list]

# adjust years for Tasman Sea events
if reg_code == 'tas':
    ev_yrs[1] = ev_yrs[1] - 0.49
    ev_yrs[2] = ev_yrs[2] + 0.49

# create figure
plt.rcParams.update({'font.size': 20})
plt.rcParams["font.family"] = "Helvetica"
plt.rc('axes', titleweight='bold')

fig, axs = plt.subplots(len(ev_list)+1, 1, figsize=(14, 4*(len(ev_list)+1)))
fig.set_facecolor('white')

ax = axs[0]

# shading of ONI (ENSO)
# rasterized=True significantly reduces the EPS file size
ax.pcolor(oni_time, oni_d, oni[:-1,:-1].T, cmap='RdBu_r', vmin=-4, vmax=4, rasterized=True)

# stacked bar plot of MHW days in each category
xvals = [x+0.5 for x in time_yr]
ax.bar(xvals, mhwBlock['moderate_days'], 0.8, color=[1, 0.81, 0.41])
ax.bar(xvals, mhwBlock['strong_days'], 0.8, bottom=mhwBlock['moderate_days'], color=[1, 0.44, 0])
ax.bar(xvals, mhwBlock['severe_days'], 0.8, \
       bottom=mhwBlock['moderate_days']+mhwBlock['strong_days'], color=[0.66, 0, 0])
ax.bar(xvals, mhwBlock['extreme_days'], 0.8, \
       bottom=mhwBlock['moderate_days']+mhwBlock['strong_days']+mhwBlock['severe_days'], \
       color=[0.2, 0, 0])
    
# plot legend
modb = mpatches.Patch(color=[1, 0.81, 0.41], label='Moderate')
strb = mpatches.Patch(color=[1, 0.44, 0], label='Strong')
sevb = mpatches.Patch(color=[0.66, 0, 0], label='Severe')
extb = mpatches.Patch(color=[0.2, 0, 0], label='Extreme')
ax.legend(handles=[modb, strb, sevb, extb],bbox_to_anchor=(1, 0.67), loc='upper left', \
          borderaxespad=0.,frameon=False, prop={'size':18,'family':'Helvetica'})

# other plot parameters
t_min = time_yr[0]
t_max = time_yr[-1]+1
ax.set_xlim(t_min, t_max)
ax.set_ylim(-49, 300)
ax.set_ylabel(r'Number of MHW days')
ax.xaxis.set_major_locator(ticker.MultipleLocator(5))
ax.xaxis.set_minor_locator(ticker.MultipleLocator(1))
ax.yaxis.set_minor_locator(ticker.MultipleLocator(50))

# plot other subplot labels
for kk in range(len(ev_yrs)):
    i1 = np.where(mhwBlock['years_centre'] == round(ev_yrs[kk]))[0][0]
    y_days = mhwBlock['total_days'][i1]
    ax.text(ev_yrs[kk]+0.5,y_days+10,'(' + chr(98+kk) + ')',ha='center',fontsize=16, fontname="Helvetica")
    
# subplot title
t_str = '(a) Annual marine heatwave days: ' + reg_name
ax.set_title(t_str,{'horizontalalignment':'left'},fontname="Helvetica",loc='left')

# legend header
ax.text(1.01, 0.7, 'MHW Categories', fontsize=18, fontname="Helvetica", fontweight="bold", transform=ax.transAxes)
    
ax.tick_params(axis='both', which='major', labelsize=18)

# plot individual MHW events
for nn, ax in enumerate(axs[1:]):
    ev = ev_list[nn]

    # set time-spans and thresholds
    t1 = np.where(t==mhws['time_start'][ev])[0][0]
    t2 = np.where(t==mhws['time_end'][ev])[0][0]
    x0 = dates[t1:t2+1]
    y0 = sst[t1:t2+1]
    y1 = thresh_cat1[t1:t2+1]
    y2 = thresh_cat2[t1:t2+1]
    y3 = thresh_cat3[t1:t2+1]
    y4 = thresh_cat4[t1:t2+1]

    # Find indices for all ten MHWs before and after event of interest and shade accordingly
    emin = ev-10
    emax = ev+10
    if emax > mhws['n_events']:
        emax = mhws['n_events']

    for ev0 in np.arange(emin, emax, 1):
        t1 = np.where(t==mhws['time_start'][ev0])[0][0]
        t2 = np.where(t==mhws['time_end'][ev0])[0][0]
        ax.fill_between(dates[t1:t2+1], sst[t1:t2+1], clim['thresh'][t1:t2+1], \
                         color=(1,0.85,0.85))

    # shade main event according to category
    ax.fill_between(x0, y0, y1, where = (y1 <= y0), color=[1, 0.81, 0.41], interpolate=True)
    ax.fill_between(x0, y0, y2, where = (y2 <= y0), color=[1, 0.44, 0], interpolate=True)
    ax.fill_between(x0, y0, y3, where = (y3 <= y0), color=[0.66, 0, 0], interpolate=True)
    ax.fill_between(x0, y0, y4, where = (y4 <= y0), color=[0.2, 0, 0], interpolate=True)

    # Plot SST, seasonal cycle, threshold, shade MHWs with main event in red
    ax.plot(dates, sst, 'k-', linewidth=3)
    ax.plot(dates, clim['thresh'], ':', color=[1, 0.81, 0.41], linewidth=1)
    ax.plot(dates, thresh_cat2, ':', color=[1, 0.44, 0], linewidth=1)
    ax.plot(dates, thresh_cat3, ':', color=[0.66, 0, 0], linewidth=1)
    ax.plot(dates, thresh_cat4, ':', color=[0.2, 0, 0], linewidth=1)
    ax.plot(dates, clim['seas'], 'b-', linewidth=2)

    # other plot parameters
    t_centre = mhws['date_start'][ev] + (mhws['date_end'][ev] - mhws['date_start'][ev]) / 2
    t_min = t_centre - timedelta(days=183)
    t_max = t_centre + timedelta(days=182)
    ax.set_xlim(t_min, t_max)
    ax.set_ylim(clim['seas'].min()-0.8, thresh_cat4.max()+0.8)
    ax.set_ylabel(r'SST ($^o$C)')
    ax.xaxis.set_major_locator(dt.MonthLocator(interval=2))
    ax.xaxis.set_major_formatter(dt.DateFormatter("%b %Y"))
    ax.xaxis.set_minor_locator(dt.MonthLocator(interval=1))
    ax.xaxis.set_minor_formatter(ticker.NullFormatter())

    # plot title
    t_str = '(' + chr(98+nn) + ') ' + mhws['date_start'][ev].strftime("%d/%m/%Y") + ' - ' + \
          mhws['date_end'][ev].strftime("%d/%m/%Y" + ' (' + str(mhws['duration'][ev]) + ' days, ' + \
          "{:.1f}".format(mhws['intensity_max'][ev]) + '$^o$C max. intensity)')
    ax.set_title(t_str,{'horizontalalignment':'left'},fontname="Helvetica",loc='left')

    # Text about ranking
    ax.text(1.01, 0.7, 'Event Rank', fontsize=18, fontname="Helvetica", fontweight="bold", transform=ax.transAxes)
    ax.text(1.01, 0.6, 'Duration: ' + str(rank['duration'][ev]), fontsize=18, transform=ax.transAxes)
    ax.text(1.01, 0.5, 'Max. Intensity: ' + str(rank['intensity_max'][ev]), fontsize=18, transform=ax.transAxes)
    ax.text(1.01, 0.4, 'Max. Severity: ' + str(rank['intensity_max_norm'][ev]), fontsize=18, transform=ax.transAxes)
    ax.text(1.01, 0.3, 'Cumul. Intensity: ' + str(rank['intensity_cumulative'][ev]), fontsize=18, transform=ax.transAxes)
    
    ax.tick_params(axis='both', which='major', labelsize=18)
    
fig.tight_layout(rect=[0,0,1,1])

# Save the plot by calling plt.savefig() BEFORE plt.show()
outfile = outpath + 'fig_case_study_' + reg_code
plt.savefig(outfile + '.png', facecolor=fig.get_facecolor(), edgecolor='none')
plt.savefig(outfile + '.eps', edgecolor='none')
plt.show()


