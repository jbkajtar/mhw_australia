function area=area_make(lon,lat)
%
% JBK 24/07/19
% generate area mesh from lon, lat
% based on area_make by A. Santoso, but does not lose 1 grid-point at edge
% and updated to use Matlab's areaquad function
% assumes the use of referenceSphere('earth','m'): spherical Earth and units in metres

% ensure x and y are column vectors, and areaquad requires doubles
lon=double(lon(:));
lat=double(lat(:));

% create lon bounds array from lon array
lonb=interp1(1:length(lon),lon,0.5:length(lon)+0.5,'linear','extrap');

% create lat bounds array from lat array
latb=interp1(1:length(lat),lat,0.5:length(lat)+0.5,'linear','extrap');
% round end points of lat array towards zero, if either end > 90
% works for flipped lat arrays too
if abs(latb(1))>90
 latb(1)=fix(latb(1));
end
if abs(latb(end))>90
 latb(end)=fix(latb(end));
end

% make 2D arrays of bound arrays
[latn,lonn]=meshgrid(latb,lonb);

% set reference sphere
rs=referenceSphere('earth','m');

% compute areas of cells with areaquad function
area=areaquad(latn(1:end-1,1:end-1),lonn(1:end-1,1:end-1),latn(2:end,2:end),lonn(2:end,2:end),rs);

return
