clear; close all;
save_or_not = 0;
%% import data
data_dir = 'E:\Research\Data\';
save_dir = 'E:\Research\Work\[else]\synoptic_map_comparison\';
% data path in order of WSO, GONG, MDI, HMI, SOLIS
obs_title = {'WSO','GONG','MDI','HMI','SOLIS'};
obs_dir = {'WSO\field\','GONG\fits\','MDI\','HMI\','SOLIS\'};
obs_name_front = {'cr','mrzqs_c','synop_Mr_0.','hmi.Synoptic_Mr.','kbv7g101210t0911c'};
obs_name_behnd = {'.dat','.fits','.fits','.fits','_000_int-mas_dim-900.fits'};
lon_res = [73,360,3600,3600,1800];
lat_res = [30,180,1080,1440, 900];
%% raw synoptic map
for cr = 2096 : 2096
    close all;
    Br_raw = cell(1,5); % [G]
    grid_raw = cell(2,5); % [deg.]

    figure();
    red_white_blue = generate_colorbar();
    LineWidth = 2;
    FontSize = 15;

    for i_obs = 1 : 4
        % read data
        obs_file = [data_dir,obs_dir{i_obs},obs_name_front{i_obs},num2str(cr),obs_name_behnd{i_obs}];
        if i_obs == 1 % WSO [.dat]
            Br_sub = importdata(obs_file); % [uT]
            Br_sub = Br_sub / 100; % [G]
            lon_sub = linspace(360,0,73);
            lat_sin = linspace(-14.5/15,14.5/15,30);
            clim = 10;
        else
            Br_sub = fitsread(obs_file); % [G]
            lon_sub = linspace(0,360,lon_res(i_obs));
            lat_sin = linspace(-1,1,lat_res(i_obs));
            clim = 30;
        end
        lat_sub = asind(lat_sin);
        [grid_raw{1, i_obs}, ~] = meshgrid(lon_sub,lat_sub);
        [~, grid_raw{2, i_obs}] = meshgrid(lon_sub,lat_sub);
        % plot figure
        subplot('Position',[0.55-0.5*mod(i_obs,2),0.5-0.5*floor((i_obs-0.1)/2),0.45,0.5])
        h = pcolor(lon_sub,lat_sub,Br_sub);
        set(h,'LineStyle','none');
        shading interp
        axis equal
        xlim([0 360])
        ylim([-90 90])

        cb = colorbar;
        colormap(red_white_blue);
        title(obs_title{i_obs},'FontSize',FontSize)
        set(gca,'CLim',[-clim clim],'TickDir','out','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
    end

    sgtitle(['CR',num2str(cr)],'FontSize',FontSize)
    annotation('textbox', [0.8, 0.8, 0.2, 0.2], 'String', 'plotted by synoptic\_map\_comparison.m', 'FitBoxToText', 'on');

    if save_or_not == 1
        save_name = ['CR',num2str(cr),'.png'];
        saveas(gca,[save_dir,save_name]);
    end
end
%% plot figure


%% function
function  red_white_blue = generate_colorbar()
    color_red   = [1,0,0];
    color_white = [1,1,1];
    color_blue  = [0,0,1];
    n1 = 100;
    n2 = 100;
    R_comp = [linspace(color_red(1),color_white(1),n1),linspace(color_white(1),color_blue(1),n2)];
    G_comp = [linspace(color_red(2),color_white(2),n1),linspace(color_white(2),color_blue(2),n2)];
    B_comp = [linspace(color_red(3),color_white(3),n1),linspace(color_white(3),color_blue(3),n2)];
    red_white_blue = [R_comp',G_comp',B_comp'];
end