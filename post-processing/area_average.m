function var1=area_average(vari,A,mode)
% 
% area average for 3D variable 
% mode = 1 : average over longitude and latitude range
% mode = 2 : average over latitude range (meridional average)
% mode = 3 : average over longitude range (zonal average)
% A. Santoso
% updated by JBK: do not need to set nans in A before using function

% ensure that A array has same nan values as vari
A(isnan(vari))=nan;

% compute area average
if(mode==1)
 var1=squeeze(nansum(nansum(vari.*A,1),2))./ ...
 squeeze(nansum(nansum(A,1),2));
elseif(mode==2)
 var1=squeeze(nansum(vari.*A,2))./ ...
 squeeze(nansum(A,2));
elseif(mode==3)
 var1=squeeze(nansum(vari.*A,1))./ ...
 squeeze(nansum(A,1));
end

