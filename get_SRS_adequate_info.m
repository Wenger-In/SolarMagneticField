clear; close all;
save_or_not = 0;
%% PART 1: Eliminate unreasonable data
srs_dir = 'E:\Research\Data\NOAA\';
data_srs_v0 = importdata([srs_dir,'gather_SRS.csv']); % original data including area=0 or lon_ext=0
srs_v1 = data_srs_v0; % eliminate  area=0 and lon_ext=0
% data profile: YYYY; MM; DD; No.; lat; lon_obs; lon_CR; area; lon_extent
i_srs = 1;
while i_srs < length(srs_v1)
    srs_sub = srs_v1(i_srs,:);
    if srs_sub(8) == 0 || srs_sub(9) == 0
        srs_v1(i_srs,:) = [];
    else
        i_srs = i_srs + 1;
    end
end
clear i_srs srs_sub
%% PART 1


%% PART 2: Separate by Carrington Rotation
cr_dir = 'E:\Research\Data\Sunspot\';
data_cr = importdata([cr_dir, 'sn_interp.dat']);
cr_lst_full = data_cr(:,1);
cr_frac_lst_full = data_cr(:,2);
year2day = 365.2422; % one year
CR2day = 27.2753; % one Carrington Rotation
% time period with SRS data
cr_lst = cr_lst_full(263:end);
cr_frac_lst = cr_frac_lst_full(263:end);
% pseudo following CR for labeled
cr_lst = [cr_lst; cr_lst(end)+1];
cr_frac_lst = [cr_frac_lst; cr_frac_lst(end)+CR2day./year2day];
% transform SRS date into fraction
srs_frac_lst = zeros(length(srs_v1),1);
for i_srs = 1 : length(srs_v1)
    srs_sub = srs_v1(i_srs,:);
    date_sub = datetime(srs_sub(1),srs_sub(2),srs_sub(3));
    doy_sub = day(date_sub,'dayofyear');
    srs_frac_lst(i_srs) = srs_sub(1) + doy_sub ./ year2day; % srs date in fraction of year
end
clear i_srs
% label SRS with CR
srs_cr_lst = zeros(length(srs_frac_lst),1);
for i_cr = 1 : length(cr_frac_lst)-1
   for i_srs = 1 : length(srs_frac_lst)
       if srs_frac_lst(i_srs)>=cr_frac_lst(i_cr) && srs_frac_lst(i_srs)<cr_frac_lst(i_cr+1)
           srs_cr_lst(i_srs) = cr_lst(i_cr);
       else
           continue
       end
   end
end
clear i_srs
%% PART 2


%% PART 3: Calculate latitudinal extent
No_lst = srs_v1(:,4);
No_lst(No_lst<4000) = No_lst(No_lst<4000) + 10000; % renamed whose No. excess 10000
lat_lst = srs_v1(:,5); % [deg.]
lon_obs_lst = srs_v1(:,6); % [deg.]
lon_lst = srs_v1(:,7); % [deg.]
area_lst = srs_v1(:,8); % [uSH]
lon_ext_lst = srs_v1(:,9); % [deg.]
lat_ext_lst = cal_lat_ext(lat_lst,area_lst,lon_ext_lst); % [deg.]
%% PART 3


%% PART 4: Organize and save data
save_dir = 'E:\Research\Data\NOAA\';
save_var = [srs_cr_lst,No_lst,lat_lst,lon_obs_lst,lon_lst,lat_ext_lst,lon_ext_lst];
if save_or_not == 1
    csvwrite([save_dir,'SRS_adequate_info.csv'],save_var);
    disp('Save NOAA SRS Adequate Infomation Successfully');
end
%% PART 4


%% functions
function lat_ext = cal_lat_ext(lat,area,lon_ext)
    deg2rad = pi./180;
    % normalize unit as [r] & [r^2]
    area_norm = area .*(2*pi*1e-6); % [r^2]
    lon_ext_norm = cosd(lat) .* lon_ext .* deg2rad; % [r]
    lat_ext_norm = area_norm ./ lon_ext_norm; % [r]
    lat_ext = lat_ext_norm ./ deg2rad;
end