% JBK 2020-06-24
% post-process data for plotting MHW statistics in the Australian region
% from NOAA OISST observations
% Start dates of strongest events by:
% - max intensity
% - duration
% Statistics of each of the strongest events
% OUTPUT: mhw_strongest.processed.aus.NOAA_OISST.AVHRR.v2-1_modified.nc

clear all;

% region bounds
reg_lab='aus';
rb=[100, 170, -50, 0];

% set paths, sourcepath -> location of mhw_cats.aus.NOAA_OISST.AVHRR.v2-1_modified.nc
sourcepath='';
inpath=sourcepath;
outpath=sourcepath;

% read data
infile=[inpath 'mhw_cats.' reg_lab '.NOAA_OISST.AVHRR.v2-1_modified.nc'];
lon=ncread(infile,'lon');
lat=ncread(infile,'lat');
time=ncread(infile,'time');
ds.ev_max_start=h5read(infile,'/ev_max_start');
ds.ev_max_max=ncread(infile,'ev_max_max');
ds.ev_max_dur=ncread(infile,'ev_max_dur');
ds.ev_dur_start=h5read(infile,'/ev_dur_start');
ds.ev_dur_max=ncread(infile,'ev_dur_max');
ds.ev_dur_dur=ncread(infile,'ev_dur_dur');

S=[length(lon),length(lat)];

% convert start date strings into datenums
ev_max_start=datenum(datetime(ds.ev_max_start));
ev_dur_start=datenum(datetime(ds.ev_dur_start));

% store required data in single array
% store in order: start date, max intensity, duration

% store start dates
mhw_str_ev(:,:,1,1)=ev_max_start;
mhw_str_ev(:,:,2,1)=ev_dur_start;

% store max intensities of each of the strongest MHWs
mhw_str_ev(:,:,1,2)=ds.ev_max_max;
mhw_str_ev(:,:,2,2)=ds.ev_dur_max;

% store durations of each of the strongest MHWs
mhw_str_ev(:,:,1,3)=ds.ev_max_dur;
mhw_str_ev(:,:,2,3)=ds.ev_dur_dur;

% write to netCDF
f1=[outpath 'mhw_strongest.processed.' reg_lab '.NOAA_OISST.AVHRR.v2-1_modified.nc'];
fmt='netcdf4_classic';

% save data to netcdf
nccreate(f1,'time', 'Dimensions',{'time',length(time)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'time',time);
ncwriteatt(f1,'time','units','years')
ncwriteatt(f1,'time','standard_name','time');
ncwriteatt(f1,'time','long_name','calendar year');
ncwriteatt(f1,'time','axis','T');
ncwriteatt(f1,'time','calendar','proleptic_gregorian');

nccreate(f1,'lon', 'Dimensions',{'lon',length(lon)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'lon',lon);
ncwriteatt(f1,'lon','standard_name','longitude');
ncwriteatt(f1,'lon','long_name','Longitude');
ncwriteatt(f1,'lon','units','degrees_east');
ncwriteatt(f1,'lon','axis','X');

nccreate(f1,'lat', 'Dimensions',{'lat',length(lat)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'lat',lat);
ncwriteatt(f1,'lat','standard_name','latitude');
ncwriteatt(f1,'lat','long_name','Latitude');
ncwriteatt(f1,'lat','units','degrees_north');
ncwriteatt(f1,'lat','axis','Y');

nccreate(f1,'mhw_str_ev', 'Dimensions',{'lon',length(lon),'lat',length(lat),'metric1',length(1:2),'metric2',length(1:3)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'mhw_str_ev',mhw_str_ev);
ncwriteatt(f1,'mhw_str_ev','units','1');
ncwriteatt(f1,'mhw_str_ev','standard_name','n/a');
ncwriteatt(f1,'mhw_str_ev','long_name','Start dates, maximum intensities, and durations of longest and strongest marine heatwaves');
ncwriteatt(f1,'mhw_str_ev','coverage_content_type','auxiliaryInformation');


