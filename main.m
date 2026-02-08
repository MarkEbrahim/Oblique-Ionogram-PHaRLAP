%
% Name :
%   main.m
%
% Purpose :
%   1. Generates oblique ionograms between a receiver site in 
%       Auburn, AL and transmitters in Chesapeake, VA and Corpus 
%       Christi, TX at 6:00 AM on April 15, 2002.
%   2. Repeats (1) at 6:00 PM on April 15, 2002.
%   3. 
%
% Latitude and Longitude:
% Auburn, AL            32.606759   -85.480299
% Chesapeake, VA        36.764497   -76.287469
% Corpus Christi, TX    27.754056   -97.429634

% Bearings and distance
% From Chesapeake to Auburn        958.6 km @ 243.866389° (from Auburn to Chesapeake is 58.623889°)
% From Corpus Christi to Auburn    1268 km @ 61.908611° (from Auburn to CC is 247.9375°)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% LOCATION CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
auburn_lat = 32.606759;          % latitude of the start point of ray
auburn_long = -85.480299;         % longitude of the start point of ray

ches_lat = 36.764497;          % latitude of the start point of ray
ches_long = 76.287469;         % longitude of the start point of ray

cc_lat = 27.754056;
cc_long = -97.429634;

ches_to_au_bearing = 243.866;
ches_to_au_distance = 958.6;

cc_to_au_bearing = 61.9086;
cc_to_au_distance = 1268;

t1 = [2026 4 15 6 0];        % UT - year, month, day, hour, minute
t2 = [2026 4 15 18 0];

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



    % iospheric environment parameters

    % max_range = 10000;      % maximum range for sampling the ionosphere (km)
    % num_range = 201;        % number of ranges (must be < 2000)
    % range_inc = max_range ./ (num_range - 1);  % range cell size (km)
    % 
    % start_height = 0 ;      % start height for ionospheric grid (km)
    % height_inc = 3;         % height increment (km)
    % num_heights = 200;      % number of  heights (must be < 2000)


fprintf( ['\n' ...
  'Example of 2D numerical raytracing for a fan of rays for a WGS84 ellipsoidal' ...
  ' Earth\n\n'])

function [all_freqs, all_heights, grid, ray_path_datas]=f1(origin_lat, origin_long, bearing, distance, datetime)

    tol = [1e-7 .01 10];         % ODE tolerance and min/max step sizes
    nhops = 1;                   % number of hops to raytrace
    doppler_flag = 1;            % generate ionosphere 5 minutes later so that
                                 % Doppler shift can be calculated
    irregs_flag = 0;             % no irregularities - not interested in 
                                 % Doppler spread or field aligned irregularities
    kp = 0;                      % kp not used as irregs_flag = 0. Set it to a 
                                 % dummy value 
    speed_of_light = 2.99792458e8;




    elevs = [2:.05:80];            % initial ray elevation
    % elevs = [2:2:60];            % initial ray elevation
    M = length(elevs);
    frequencies = [1:.1:10];
    % frequencies = [1:.1:2];

    % grids = struct(size(frequencies));


    cfg = config();


    
    all_heights_temp = zeros(size(frequencies));
    idxs = length(frequencies);

    clear iri_options
    iri_options.Ne_B0B1_model = 'Bil-2000'; % this is a non-standard setting for 
                                                % IRI but is used as an example
    fprintf('Generating ionospheric grid... ')

    tic
    %
    % generate ionospheric, geomagnetic and irregularity grids
    %
    [iono_pf_grid, iono_pf_grid_5, collision_freq, irreg] = ...
        gen_iono_grid_2d(origin_lat, origin_long, cfg.R12, datetime, bearing, ...
                         cfg.max_range, cfg.num_range, cfg.range_inc, cfg.start_height, ...
	             cfg.height_inc, cfg.num_heights, kp, doppler_flag, 'iri2020', ...
	             iri_options);
    toc
    grid = iono_pf_grid;

    C = {};
    % freq is in MHz
    for idx = 1:idxs
        freq = frequencies(idx);


                
        % convert plasma frequency grid to  electron density in electrons/cm^3
        iono_en_grid = iono_pf_grid.^2 / 80.6164e-6;
        iono_en_grid_5 = iono_pf_grid_5.^2 / 80.6164e-6;
        
        

        freqs = freq.*ones(size(elevs));
        % call raytrace for a fan of rays
        % first call to raytrace so pass in the ionospheric and geomagnetic grids 
        % fprintf('Generating %d 2D NRT rays ...', M);
        tic
        [ray_data, ray_path_data] = ...
           raytrace_2d(origin_lat, origin_long, elevs, bearing, freqs, nhops, ...
                       tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
 	               collision_freq, cfg.start_height, cfg.height_inc, cfg.range_inc, irreg);
        toc;
        
        

        % find frequency that travels closest to Auburn
        % Identify the frequency that travels closest to Auburn
        mask  = arrayfun(@(s) ~isnan(s.virtual_height), ray_data);
        ground_ranges  = arrayfun(@(s) s.ground_range, ray_data);
        ground_ranges2 = [ray_data.ground_range];
        
        ground_ranges(~mask) = Inf;          % exclude non-matching by making them huge
        [m, index] = min(abs(ground_ranges - distance));

        % ground_ranges = [ray_data.ground_range];
        % data111 = [raydata.ground_range; raydata.apogee];
        % x = abs(ground_ranges - distance)
        % [m, index] = min((x)(x));
        height_for_elevation_angle = ray_data(index).virtual_height;
        closest_elevation_angle = ray_data(index).initial_elev;
        fprintf('Freq %.2f: Closest elevation angle to Auburn: %.2f deg (d: %.2f); height: %.2f km \n', ...
            freq, closest_elevation_angle, m, height_for_elevation_angle);
        ray_path_data_here = ray_path_data(index);
        C{end + 1} = ray_path_data_here;
        all_heights_temp(idx) = height_for_elevation_angle;

    end

    all_freqs = frequencies;
    all_heights = all_heights_temp;
    ray_path_datas = C;
end




% lll = cell2mat(jjj);




% figure(2)
% plot(fff, hhh)


function plotionogram(x, y, title_str)
    
    figure;

    ax = axes;
    hold(ax, 'on');
    grid(ax, 'on');
    box(ax, 'on');
    
    scatter(x, y, 25, 'filled', 'o');
    
    xlabel('Frequency (MHz)', ...
        'FontSize', 12, ...
        'FontWeight', 'bold');

    ylabel('Virtual Reflection Height (km)', ...
        'FontSize', 12, ...
        'FontWeight', 'bold');

    title(title_str, ...
        'FontSize', 14, ...
        'FontWeight', 'bold');

    set(ax, ...
        'FontSize', 11, ...
        'LineWidth', 1.2, ...
        'GridLineStyle', '--', ...
        'MinorGridLineStyle', ':', ...
        'XMinorGrid', 'on', ...
        'YMinorGrid', 'on');
    ylim([0 400])
end

function altitude_graph(lat, long, bearing, tim, grid, raypathdatas)
    cfg = config();
    lll = cell2mat(raypathdatas);
    % plot the rays and ionosphere
    figure(1)
    UT_str = [num2str(tim(3)) '/' num2str(tim(2)) '/' num2str(tim(1)) '  ' ...
              num2str(tim(4), '%2.2d') ':' num2str(tim(5), '%2.2d') 'UT'];
    % freq_str = [num2str(freq) 'MHz'];
    R12_str = num2str(cfg.R12);
    lat_str = num2str(lat);
    lon_str = num2str(long);
    bearing_str = num2str(bearing);
    fig_str = [UT_str '  R12 = ' R12_str '   lat = ' lat_str ...
               ', lon = ' lon_str ', bearing = ' bearing_str];
    set(gcf, 'name', fig_str)
    start_range = 0;
    end_range = cfg.max_range;
    end_range_idx = fix((end_range-start_range) ./ cfg.range_inc) + 1;
    start_ht = cfg.start_height;
    start_ht_idx = 1;
    end_ht = 400;
    end_ht_idx = fix(end_ht ./ cfg.height_inc) + 1;
    iono_pf_subgrid = grid(start_ht_idx:end_ht_idx, 1:end_range_idx);
    plot_ray_iono_slice(iono_pf_subgrid, start_range, end_range, cfg.range_inc, ...
        start_ht, end_ht, cfg.height_inc, lll, 'color', [1, 1, 0.99], ...
        'linewidth', .5);

    set(gcf,'units','normal')
    pos = get(gcf,'position');
    pos(2) = 0.55;
    set(gcf,'position', pos)

    % uncomment the following to print figure to hi-res ecapsulated postscript
    % and PNG files
    set(gcf, 'paperorientation', 'portrait')
    set(gcf, 'paperunits', 'cent', 'paperposition', [0 0 61 18])
    set(gcf, 'papertype', 'a4') 
    print -depsc2 -loose -opengl test.ps 
    print -dpng test.png
end

% Chesapeake to Auburn, 6 AM
[frequencies1, heights1, grid1, raypathdatas1] = f1(ches_lat, ches_long, ches_to_au_bearing, ches_to_au_distance, t1);
plotionogram(frequencies1, heights1, 'Virtual Reflection Height vs Frequency (04/15/2002 6:00 AM from Chesapeake, VA to Auburn, AL)')

% Chesapeake to Auburn, 6 PM
[frequencies2, heights2, grid2, raypathdatas2] = f1(ches_lat, ches_long, ches_to_au_bearing, ches_to_au_distance, t2);
plotionogram(frequencies2, heights2, 'Virtual Reflection Height vs Frequency (04/15/2002 6:00 PM from Chesapeake, VA to Auburn, AL)')

% [frequencies3, heights3, grid3, raypathdatas3] = f1(cc_lat, cc_long, cc_to_au_bearing, cc_to_au_distance, t1);
% altitude_graph(cc_lat, cc_long, cc_to_au_bearing, t1, grid3, raypathdatas3);
% plotionogram(frequencies3, heights3, 'Virtual Reflection Height vs Frequency (04/15/2002 6:00 AM from Corpus Christi, TX to Auburn, AL)')

% [frequencies4, heights4, grid4, raypathdatas4] = f1(cc_lat, cc_long, cc_to_au_bearing, cc_to_au_distance, t2);
% altitude_graph(cc_lat, cc_long, cc_to_au_bearing, t2, grid4, raypathdatas4);
% plotionogram(frequencies4, heights4, 'Virtual Reflection Height vs Frequency (04/15/2002 6:00 PM from Corpus Christi, TX to Auburn, AL)')


% figure('Color','w');            % White background (journal-friendly)

% ax = axes;
% hold(ax,'on')
% grid(ax,'on')
% box(ax,'on')
% 
% scatter(frequencies1, heights1, 50, 'filled')
% scatter(frequencies2, heights2, 50, 'filled')
% 
% xlabel('X axis')
% ylabel('Y axis')
% legend({'Dataset A','Dataset B'}, 'Location','best')

% Corpus Christi to Auburn, 6 AM

% Corpus Christi to Auburn, 6 PM



