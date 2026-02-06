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

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% OTHER CONSTANTS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% PARAMETERS
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

UT = [2026 4 15 6 0];        % UT - year, month, day, hour, minute


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




    elevs = [2:.05:60];            % initial ray elevation
    M = length(elevs);
    frequencies = [1:.1:10];

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
        fprintf('Generating %d 2D NRT rays ...', M);
        tic
        [ray_data, ray_path_data] = ...
           raytrace_2d(origin_lat, origin_long, elevs, bearing, freqs, nhops, ...
                       tol, irregs_flag, iono_en_grid, iono_en_grid_5, ...
 	               collision_freq, cfg.start_height, cfg.height_inc, cfg.range_inc, irreg);
        toc;
        
        

        % find frequency that travels closest to Auburn
        % Identify the frequency that travels closest to Auburn
        ground_ranges = [ray_data.ground_range];
        [~, index] = min(abs(ground_ranges - distance));
        height_for_elevation_angle = ray_data(index).apogee;
        closest_elevation_angle = ray_data(index).initial_elev;
        fprintf('Closest elevation angle to Auburn: %.2f deg; height: %.2f km \n', ...
            closest_elevation_angle, height_for_elevation_angle);
        ray_path_data_here = ray_path_data(index);
        C{end + 1} = ray_path_data_here;
        all_heights_temp(idx) = height_for_elevation_angle;

    end

    all_freqs = frequencies;
    all_heights = all_heights_temp;
    ray_path_datas = C;
end


[fff, hhh, iii, jjj] = f1(ches_lat, ches_long, ches_to_au_bearing, ches_to_au_distance, UT);


lll = cell2mat(jjj);

cfg = config();

% plot the rays and ionosphere
figure(1)
UT_str = [num2str(UT(3)) '/' num2str(UT(2)) '/' num2str(UT(1)) '  ' ...
          num2str(UT(4), '%2.2d') ':' num2str(UT(5), '%2.2d') 'UT'];
% freq_str = [num2str(freq) 'MHz'];
R12_str = num2str(cfg.R12);
lat_str = num2str(ches_lat);
lon_str = num2str(ches_long);
bearing_str = num2str(ches_to_au_bearing);
fig_str = [UT_str '  R12 = ' R12_str '   lat = ' lat_str ...
           ', lon = ' lon_str ', bearing = ' bearing_str];
set(gcf, 'name', fig_str)
start_range = 0;
end_range = cfg.max_range;
end_range_idx = fix((end_range-start_range) ./ cfg.range_inc) + 1;
start_ht = start_height;
start_ht_idx = 1;
end_ht = 400;
end_ht_idx = fix(end_ht ./ cfg.height_inc) + 1;
iono_pf_subgrid = iii(start_ht_idx:end_ht_idx, 1:end_range_idx);
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
% print -depsc2 -loose -opengl test.ps 
% print -dpng test.png

figure(2)
plot(fff, hhh)