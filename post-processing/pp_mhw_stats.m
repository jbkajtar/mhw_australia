% JBK 2020-06-11
% post-process data for plotting MHW statistics in the Australian region
% from NOAA OISST observations
% Mean MHW frequency, intensity, and duration
% Trends in MHW frequency, intensity, and duration
% Time-series of MHW stats
% Days in each MHW category
% Time-series of MHW categories
% OUTPUT: mhw_stats.processed.aus.NOAA_OISST.AVHRR.v2-1_modified.nc

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
ds.mhw_cats_dpy=ncread(infile,'mhw_cats_dpy');
ds.mhw_count=ncread(infile,'mhw_count');
ds.mhw_intensity=ncread(infile,'mhw_intensity');
ds.mhw_duration=ncread(infile,'mhw_duration');
ds.mhw_count_tr=ncread(infile,'mhw_count_tr');
ds.mhw_intensity_tr=ncread(infile,'mhw_intensity_tr');
ds.mhw_duration_tr=ncread(infile,'mhw_duration_tr');

% store data in global array, keep only 1982-2020
[k1,k2]=findrange(time,1982,2020);
mhw_stats_fld(:,:,1)=nanmean(ds.mhw_count(:,:,k1:k2),3);
mhw_stats_fld(:,:,2)=nanmean(ds.mhw_intensity(:,:,k1:k2),3);
mhw_stats_fld(:,:,3)=nanmean(ds.mhw_duration(:,:,k1:k2),3);
mhw_stats_fld(:,:,4)=10*ds.mhw_count_tr;    % trends/yr -> trends/decade
mhw_stats_fld(:,:,5)=10*ds.mhw_intensity_tr;
mhw_stats_fld(:,:,6)=10*ds.mhw_duration_tr;
time=time(k1:k2);

% trim to required region
[i1,i2]=findrange(lon,rb(1),rb(2));
[j1,j2]=findrange(lat,rb(3),rb(4));

% create AREA array
AREA=repmat(area_make(lon(i1:i2),lat(j1:j2)),[1 1 length(time)]);

% compute area averages, note: lat is in first position in data source
mhw_stats_ts(:,1)=area_average(ds.mhw_count(i1:i2,j1:j2,k1:k2),AREA,1);
mhw_stats_ts(:,2)=area_average(ds.mhw_intensity(i1:i2,j1:j2,k1:k2),AREA,1);
mhw_stats_ts(:,3)=area_average(ds.mhw_duration(i1:i2,j1:j2,k1:k2),AREA,1);

% set periods of category blocks
yr_blk=[1982,2001; 2001,2020];
clear yr_str;
for jj=1:size(yr_blk,1)
 yr_str{jj}=[num2str(yr_blk(jj,1)) '-' num2str(yr_blk(jj,2))];
end

% store sum of days in different catergories, in each block
for ii=1:2
 [i1,i2]=findrange(time,yr_blk(ii,1),yr_blk(ii,2));
 for kk=1:4
  mhw_cats(:,:,kk,ii)=mean(sum(ds.mhw_cats_dpy(:,:,kk:4,i1:i2),3),4);
 end
end

% store the number of ocean cells for the region
mask_ocean=~isnan(nanmean(ds.mhw_cats_dpy(:,:,1,:),4));
mask_land=isnan(nanmean(ds.mhw_cats_dpy(:,:,1,:),4));
num_oc=numel(find(mask_ocean));   % number of ocean cells
num_lc=numel(find(mask_land));  % number of land cells

% store the MHW maximum category for each year
mhw_max_cat=double(ds.mhw_cats_dpy>0);

% remove higher categories events from lower ones, but then need to set any '-1' values to zero,
% since sometimes a higher category MHW is recorded without passing through all lower categories
for kk=3:-1:1
 mhw1=mhw_max_cat(:,:,kk,:)-sum(mhw_max_cat(:,:,kk+1:4,:),3);
 mhw1(mhw1<0)=0;
 mhw_max_cat(:,:,kk,:)=mhw1;
end

% regions where no MHWs were recorded
mhw_max_cat(:,:,5,:)=double(sum(ds.mhw_cats_dpy,3)==0);

% apply ocean mask
om=double(mask_ocean);
om(om==0)=nan;
om=repmat(om,[1 1 size(mhw_max_cat,3) size(mhw_max_cat,4)]);
mhw_max_cat=mhw_max_cat.*om;
clear om;

% create AREA array
AREA=area_make(lon,lat);
AREA=repmat(AREA,[1 1 size(mhw_max_cat,3) size(mhw_max_cat,4)]);

% take area average of MHW max catergory array
mhw_frac_ts=area_average(mhw_max_cat,AREA,1);
mhw_frac_ts(6,:)=sum(mhw_frac_ts(1:5,:),1);

% reshape array
mhw_frac_ts=mhw_frac_ts';

% additional variables for netCDF storage
mhw_categories=[1:4];
time_block=yr_blk;

% write to netCDF
f1=[outpath 'mhw_stats.processed.' reg_lab '.NOAA_OISST.AVHRR.v2-1_modified.nc'];
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

nccreate(f1,'cat', 'Dimensions',{'cat',length(mhw_categories)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'cat',mhw_categories);
ncwriteatt(f1,'cat','units','1');
ncwriteatt(f1,'cat','long_name','Marine heatwave category following Hobday et al. (2018) definition');
ncwriteatt(f1,'cat','mhw_categories','1: Moderate, 2: Strong, 3: Severe, 4: Extreme');

nccreate(f1,'time_block', 'Dimensions',{'time_block',size(time_block,1),'time_bnds',size(time_block,2)}, 'Datatype','single','Format',fmt);
ncwrite(f1,'time_block',time_block);
ncwriteatt(f1,'time_block','units','years')
ncwriteatt(f1,'time_block','standard_name','time');
ncwriteatt(f1,'time_block','long_name','Bounds of time averaging blocks');
ncwriteatt(f1,'time_block','axis','T');

nccreate(f1,'mhw_cats', 'Dimensions',{'lon',length(lon),'lat',length(lat),'cat',length(mhw_categories), 'time_block',size(time_block,1)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'mhw_cats',mhw_cats);
ncwriteatt(f1,'mhw_cats','units','1');
ncwriteatt(f1,'mhw_cats','standard_name','n/a');
ncwriteatt(f1,'mhw_cats','long_name','Count of days per year (dpy) in each marine heatwave category, averaged in each time block');
ncwriteatt(f1,'mhw_cats','coverage_content_type','auxiliaryInformation');

nccreate(f1,'mhw_stats_fld', 'Dimensions',{'lon',length(lon),'lat',length(lat),'metric',length(1:6)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'mhw_stats_fld',mhw_stats_fld);
ncwriteatt(f1,'mhw_stats_fld','units','1');
ncwriteatt(f1,'mhw_stats_fld','standard_name','n/a');
ncwriteatt(f1,'mhw_stats_fld','long_name','Long-term means and trends of marine heatwave metrics');
ncwriteatt(f1,'mhw_stats_fld','coverage_content_type','auxiliaryInformation');

nccreate(f1,'mhw_stats_ts', 'Dimensions',{'time',length(time),'metric2',length(1:3)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'mhw_stats_ts',mhw_stats_ts);
ncwriteatt(f1,'mhw_stats_ts','units','1');
ncwriteatt(f1,'mhw_stats_ts','standard_name','n/a');
ncwriteatt(f1,'mhw_stats_ts','long_name','Annual region-wide spatial means of marine heatwave metrics');
ncwriteatt(f1,'mhw_stats_ts','coverage_content_type','auxiliaryInformation');

nccreate(f1,'mhw_frac_ts', 'Dimensions',{'time',length(time),'metric3',length(1:6)}, 'Datatype','single', 'Format',fmt, 'DeflateLevel',2);
ncwrite(f1,'mhw_frac_ts',mhw_frac_ts);
ncwriteatt(f1,'mhw_frac_ts','units','1');
ncwriteatt(f1,'mhw_frac_ts','standard_name','n/a');
ncwriteatt(f1,'mhw_frac_ts','long_name','Annual fraction of ocean at maximum marine heatwave category');
ncwriteatt(f1,'mhw_frac_ts','coverage_content_type','auxiliaryInformation');

ncwriteatt(f1,'/','source_code','https://github.com/jbkajtar/mhw_australian');
ncwriteatt(f1,'/','title','Marine heatwave statistics for the Australian region (100E-170E, 50S-0)');
ncwriteatt(f1,'/','summary','Data generated for Kajtar et al., ''A catalogue of marine heatwave characteristics and trends for the Australian region'', (2021)');
ncwriteatt(f1,'/','keywords','marine heatwave; extreme event; impact; ocean warming; Australia; observational');
ncwriteatt(f1,'/','Conventions','ACDD-1.3');


