clear; close all;
save_or_not = 0;
store_dir = 'E:\Research\Data\NOAA\SRS\';
%% PART 0: Specify date of SRS
date_beg = datenum(1996,01,02);
date_end = datenum(2022,09,01);
date_gap = 1;
date_lst = datestr(date_beg:date_gap:date_end,'yyyymmdd');
date_num = length(date_lst);
%% PART 0


%% PART 1: organize data
ar_info = [];
for i_date = 1 : date_num
    date = date_lst(i_date,:);
    file_name = [store_dir,date,'SRS.txt'];
    % check whether file for this date exists
    if exist(file_name,'file') == 0
    else
        % import data
        file_id = fopen(file_name);
        data_scan = textscan(file_id,'%s');
        fclose(file_id);
        % extract data
        ar_item = 8;
        ar_info_sub = [];
        data_cell = data_scan{1,1};
        for i_find = 55 : 70
            data_find = data_cell{i_find,1};
            if strcmp(data_find,'Nmbr') || strcmp(data_find,'NMBR')
                i_row = i_find + 9; % where AR data begins
            end
        end
        while i_row
            beg_str = data_cell{i_row,1};
            if isempty(str2num(beg_str))
                break
            else
                year = str2num(date(1:4));
                month = str2num(date(5:6));
                day = str2num(date(7:8));
                number = str2num(data_cell{i_row,1});
                lat_lon_obs = data_cell{i_row+1,1};
                [lat_obs,lon_obs] = symbol_lat_lon(lat_lon_obs);
                lon = str2num(data_cell{i_row+2,1});
                area = str2num(data_cell{i_row+3,1});
                lon_ext = str2num(data_cell{i_row+5,1});
                ar_info_subb = [year,month,day,number,lat_obs,lon_obs,lon,area,lon_ext];
            end
            ar_info_sub = [ar_info_sub;ar_info_subb];
            i_row = i_row + ar_item;
        end
    end
    ar_info = [ar_info;ar_info_sub];
end
%% PART 1


%% PART 2: save data
save_dir = 'E:\Research\Data\NOAA\';
save_var = ar_info;
if save_or_not == 1
%     save_file = [save_dir,'gather_SRS.dat'];
%     save(save_file,'save_var','-ascii');
    csvwrite([save_dir,'gather_SRS.csv'],save_var);
    disp('Neaten NOAA Solar Region Summary Table Successfully');
end
%% PART 2


%% functions
function [lat,lon] = symbol_lat_lon(lat_lon_str)
    % symbolize latitude data
    lat_str = lat_lon_str(1:3);
    if lat_str(1) == 'N'
        lat = str2num(lat_str(2:3));
    elseif lat_str(1) == 'S'
        lat = -str2num(lat_str(2:3));
    end
    % symbolize longitude data
    lon_str = lat_lon_str(4:6);
    if lon_str(1) == 'E'
        lon = str2num(lon_str(2:3));
    elseif lon_str(1) == 'W'
        lon = 360 - str2num(lon_str(2:3));
    end
end