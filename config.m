function cfg = config()

    cfg.max_range = 2000;      % maximum range for sampling the ionosphere (km)
    cfg.num_range = 201;        % number of ranges (must be < 2000)
    cfg.range_inc = cfg.max_range ./ (cfg.num_range - 1);  % range cell size (km)
    
    cfg.start_height = 0 ;      % start height for ionospheric grid (km)
    cfg.height_inc = 1;         % height increment (km)
    cfg.num_heights = 1000;      % number of  heights (must be < 2000)

    cfg.R12 = 100;                   % R12 index

end