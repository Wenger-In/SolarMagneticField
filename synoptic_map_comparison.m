clear; close all;
save_or_not = 0;
%% import data
data_dir = 'E:\Research\Data\';
save_dir = 'E:\Research\Work\[else]\synoptic_map_comparison\';
% data path in order of WSO, GONG, SOLIS, MDI, HMI
obs_title = {'WSO','GONG','SOLIS','MDI','HMI'};
obs_dir = {'WSO\field\','GONG\fits\','SOLIS\','MDI\','HMI\'};
obs_name_front = {'cr','mrzqs_c','cr','synop_Mr_0.','hmi.Synoptic_Mr.'};
obs_name_behnd = {'.dat','.fits','.fits','.fits','.fits'};
lon_res = [73,360,1800,3600,3600];
lat_res = [30,180, 900,1080,1440];

red_white_blue = generate_colorbar();
LineWidth = 2;
FontSize = 15;
figure_width = 1600;
figure_height = 900;

for cr = 2097 : 2104
    close all;

    %% raw synoptic map
    Br_raw = cell(1,5); % [G]
    grid_raw = cell(2,5); % [deg.]

    % read data
    for i_obs = 1 : 5
        obs_file = [data_dir,obs_dir{i_obs},obs_name_front{i_obs},num2str(cr),obs_name_behnd{i_obs}];
        if i_obs == 1 % WSO [.dat]
            Br_sub = importdata(obs_file); % [uT]
            Br_sub = Br_sub / 100; % [G]
            Br_sub = flipud(Br_sub); % filp up-down
            lon_sub = linspace(360,0,73);
            lat_sin = linspace(-14.5/15,14.5/15,30);
        else
            Br_sub = fitsread(obs_file); % [G]
            lon_sub = linspace(0,360,lon_res(i_obs));
            lat_sin = linspace(-1,1,lat_res(i_obs));
        end
        lat_sub = asind(lat_sin);
        Br_raw{i_obs} = Br_sub;
        [grid_raw{1, i_obs}, ~] = meshgrid(lon_sub,lat_sub);
        [~, grid_raw{2, i_obs}] = meshgrid(lon_sub,lat_sub);
    end

    % plot figure
    figure();
    for i_obs = 1 : 5
        if i_obs == 1
            clim = 10;
        else
            clim = 30;
        end
        subplot('Position',[0.03+0.33*mod(i_obs-1,3),0.5-0.5*floor((i_obs-0.1)/3),0.3,0.5])
%         subplot(2,3,i_obs)
        h = pcolor(grid_raw{1, i_obs},grid_raw{2, i_obs},Br_raw{i_obs});
        set(h,'LineStyle','none');
        shading interp
        axis equal
        xlim([0 360])
        ylim([-90 90])

        colorbar('Location', 'SouthOutside')
        colormap(red_white_blue);
        title(['CR',num2str(cr),'-',obs_title{i_obs}],'FontSize',FontSize)
        set(gca,'CLim',[-clim clim],'TickDir','out','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
    end

    annotation('textbox', [0.8, 0.8, 0.2, 0.2], 'String', 'plotted by synoptic\_map\_comparison.m', 'FitBoxToText', 'on');
    set(gcf,'Position',[10,10,figure_width,figure_height])

    if save_or_not == 1
        save_name = ['CR',num2str(cr),'_raw.png'];
        saveas(gca,[save_dir,save_name]);
    end

    %% interpolated synoptic map
    % standard grid
    hmi_lon = grid_raw{1,5};
    hmi_lat = grid_raw{2,5};
    std_index = find(abs(sind(hmi_lat))<=14.5/15);
    std_lon = reshape(hmi_lon(std_index), [], 3600);
    std_lat = reshape(hmi_lat(std_index), [], 3600);

    % interpolate
    Br_interp = cell(1,5);
    for i_obs = 1 : 5
        Br_interp_sub = interp2(grid_raw{1,i_obs}, grid_raw{2,i_obs}, Br_raw{i_obs}, std_lon, std_lat, 'linear');
        Br_interp{i_obs} = Br_interp_sub;
    end

    % plot figure
    figure();
    for i_obs = 1 : 5
        if i_obs == 1 % WSO
            clim = 10;
        else
            clim = 30;
        end
        subplot('Position',[0.03+0.33*mod(i_obs-1,3),0.5-0.5*floor((i_obs-0.1)/3),0.3,0.5])
%         subplot(2,3,i_obs)
        h = pcolor(std_lon,std_lat,Br_interp{i_obs});
        set(h,'LineStyle','none');
        shading interp
        axis equal
        xlim([0 360])
        ylim([-90 90])

        colorbar('Location', 'SouthOutside')
        colormap(red_white_blue);
        title([obs_title{i_obs},' interp','-','CR',num2str(cr)],'FontSize',FontSize)
        set(gca,'CLim',[-clim clim],'TickDir','out','XminorTick','on','YminorTick','on','LineWidth',LineWidth,'FontSize',FontSize);
    end

    annotation('textbox', [0.8, 0.8, 0.2, 0.2], 'String', 'plotted by synoptic\_map\_comparison.m', 'FitBoxToText', 'on');
    set(gcf,'Position',[10,10,figure_width,figure_height])

    if save_or_not == 1
        save_name = ['CR',num2str(cr),'_interp.png'];
        saveas(gca,[save_dir,save_name]);
    end

    %% comparison to HMI
    figure();
    for i_obs = 1 : 4
        % plot scatter
        subplot(2,2,i_obs)
        Br_interp_sub = Br_interp{i_obs}(:);
        Br_interp_hmi = Br_interp{5}(:);
        scatter(Br_interp_hmi,Br_interp_sub,2,'filled')
        grid on; hold on

        % regression
        nan_index = isnan(Br_interp_hmi) | isnan(Br_interp_sub);
        Br_interp_hmi(nan_index) = 0;
        Br_interp_sub(nan_index) = 0;
        fit = polyfit(Br_interp_hmi, Br_interp_sub, 1);
        Br_fit_sub = polyval(fit, Br_interp_hmi);
        plot(Br_interp_hmi,Br_fit_sub,'r','LineWidth',LineWidth)
        % calculate RMSE
        cc = corrcoef(Br_interp_sub, Br_fit_sub);
        legend(['cc=',num2str(cc(2))],['k=',num2str(fit(1))],'Location','southeast')

        xlabel('HMI [G]')
        ylabel([obs_title{i_obs},'[G]'])
        sgtitle(['CR',num2str(cr)],'FontSize',FontSize)
        set(gca,'LineWidth',LineWidth,'FontSize',FontSize)
    end

    annotation('textbox', [0.8, 0.8, 0.2, 0.2], 'String', 'plotted by synoptic\_map\_comparison.m', 'FitBoxToText', 'on');
    set(gcf,'Position',[10,10,figure_width,figure_height])

    if save_or_not == 1
        save_name = ['CR',num2str(cr),'_comparison.png'];
        saveas(gca,[save_dir,save_name]);
    end
end


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