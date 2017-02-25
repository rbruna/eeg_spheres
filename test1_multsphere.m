%%
clc
clear
close all

% Adds the functions folder to the path.
addpath ( sprintf ( '%s/functions', fileparts ( pwd ) ) );

% Adds the functions folder to the path (in development).
addpath ( sprintf ( '%s/functions', pwd ) );

% Adds FT to the path.
ft_path;
ft_defaults

ft_hastoolbox ( 'freesurfer', 1, 1 );
ft_hastoolbox ( 'spm8', 1, 1 );

clear

% Loads the EEG data.
eegdata = load ( '../../data/alc02_restingOA_EEG' );
eegdata.trialdata.elec = ft_convert_units ( eegdata.trialdata.elec, 'm' );
eegdata.trialdata.grad = ft_convert_units ( eegdata.trialdata.grad, 'm' );
eegdata.trialdata.elec.elecpos ( isnan ( eegdata.trialdata.elec.elecpos ) ) = 0;

% Loads the subject segmentation and the generated meshes.
mridata = load ( '../../data/alc02_mri.mat' );
mridata.mesh = ft_convert_units ( mridata.mesh, 'mm' );
mridata.mesh = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.mesh );
mridata.mesh = ft_convert_units ( mridata.mesh, 'm' );
mridata.grid = ft_convert_units ( mridata.grid, 'mm' );
mridata.grid = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.grid );
mridata.grid = ft_convert_units ( mridata.grid, 'm' );

% Fits the three concentric spheres to the meshes.
headmodel    = ft_headmodel_concentricspheres ( mridata.mesh );
headmodelsc  = mymcs_headmodel ( mridata.mesh, eegdata.trialdata.elec, 'singlesphere', 'yes' );
headmodellc  = mymcs_headmodel ( mridata.mesh, eegdata.trialdata.elec );
headmodellc2 = mymcs_headmodel ( mridata.mesh, eegdata.trialdata.elec, 'eqradii', false );

% Selects only the 60 first EEG channels.
sens   = eegdata.trialdata.elec;


%%

% For eachs ensor, draws the head meshes and the spheres.

% If EEG projects the channels onto the outter mesh.
if isfield ( sens, 'elecpos'  )
    for eindex = 1: 60
        [ ~, Pm ] = NFT_dmp ( sens.elecpos ( eindex, : ), mridata.mesh.bnd ( end ).pnt, mridata.mesh.bnd ( end ).tri );
        sens.elecpos ( eindex, : ) = Pm;
    end
    for cindex = 1: 60
        [ ~, Pm ] = NFT_dmp ( sens.chanpos ( cindex, : ), mridata.mesh.bnd ( end ).pnt, mridata.mesh.bnd ( end ).tri );
        sens.chanpos ( cindex, : ) = Pm;
    end
end

for sindex = 1: size ( sens.label, 1 )
    
    label  = sens.label { sindex };
    hindex = find ( strcmp ( headmodellc.label, label ) );
    
    if isempty ( hindex )
        continue
    end
    
    % Creats the figure.
    figure
    hold on
    
    % Plots all the sensors.
    sensors = sens.chanpos;
    plot3 ( sensors ( :, 1 ), sensors ( :, 2 ), sensors ( :, 3 ), 'k.' )
    
    % Plots the sensor.
    sensor = sens.chanpos ( hindex, : );
    plot3 ( sensor (1), sensor (2), sensor (3), 'r*' )
    
    % Plots the three meshes.
    ft_plot_mesh ( mridata.mesh.bnd (1), 'edgecolor', 'none', 'facecolor', 'brain', 'facealpha', 0.1 );
    ft_plot_mesh ( mridata.mesh.bnd (2), 'edgecolor', 'none', 'facecolor', 'white', 'facealpha', 0.1 );
    ft_plot_mesh ( mridata.mesh.bnd (3), 'edgecolor', 'none', 'facecolor', 'skin',  'facealpha', 0.3 );
    
    % Plots the three spheres.
    spheres      = headmodellc;
    spheres.o    = headmodellc.o ( hindex, : );
    spheres.type = 'concentricspheres';
    spheres.r    = headmodellc.r ( hindex, 1 );
    ft_plot_vol ( spheres, 'edgecolor', 'none', 'facecolor', 'brain', 'facealpha', 0.8 )
    spheres.r    = headmodellc.r ( hindex, 2 );
    ft_plot_vol ( spheres, 'edgecolor', 'none',  'facecolor', 'white', 'facealpha', 0.6 )
    spheres.r    = headmodellc.r ( sindex, 3 );
    ft_plot_vol ( spheres, 'edgecolor', 'none',  'facecolor', 'skin',  'facealpha', 0.3 )
    
%     plot3 ( mridata.grid.pos ( [ 719 720 ], 1 ), mridata.grid.pos ( [ 719 720 ], 2 ), mridata.grid.pos ( [ 719 720 ], 3 ), '*b' )
    
%     views = linspace ( 0, 360, 100 );
%     views = views ( 1: end - 1 );
%     for vindex = 1: numel ( views )
%         view ( views ( vindex ), 0 );
%         pause ( 0.04 );
%     end
%     close
    
%     savefig ( sprintf ( '%s.fig', sens.label { sindex } ) );
%     my_savegif ( sprintf ( '%s.gif', sens.label { sindex } ) );
%     close
    
    rotate3d
    uiwait
end

%%

% For each sensor, draws the contribution of each source.

% If EEG projects the channels onto the outter mesh.
if isfield ( sens, 'elecpos'  )
    for eindex = 1: 60
        [ ~, Pm ] = NFT_dmp ( sens.elecpos ( eindex, : ), mridata.mesh.bnd ( end ).pnt, mridata.mesh.bnd ( end ).tri );
        sens.elecpos ( eindex, : ) = Pm;
    end
    for cindex = 1: 60
        [ ~, Pm ] = NFT_dmp ( sens.chanpos ( cindex, : ), mridata.mesh.bnd ( end ).pnt, mridata.mesh.bnd ( end ).tri );
        sens.chanpos ( cindex, : ) = Pm;
    end
end

cfg = [];
cfg.grid = mridata.grid;
cfg.senstype = 'eeg';
cfg.headmodel = headmodel;
cfg.channel = { 'all' '-EEG061' '-EEG062' '-EEG063' '-EEG064' };
cfg.elec = eegdata.trialdata.elec;
cfg.sens = eegdata.trialdata.elec;
cfg.feedback = 'no';

tic
leadfieldX = my_leadfield ( cfg );
leadfieldX = cat ( 2, leadfieldX.leadfield {:} );
leadfieldX = leadfieldX';
toc

cfg2 = cfg;
cfg2.headmodel = headmodellc;

tic
leadfieldY = my_leadfield ( cfg2 );
leadfieldY = cat ( 2, leadfieldY.leadfield {:} );
leadfieldY = leadfieldY';
toc

cfg3 = cfg;
cfg3.headmodel = headmodellc2;

tic
leadfieldZ = my_leadfield ( cfg3 );
leadfieldZ = cat ( 2, leadfieldZ.leadfield {:} );
leadfieldZ = leadfieldZ';
toc

for sindex = 1: size ( sens.label, 1 )
    
    label  = sens.label { sindex };
    hindex = find ( strcmp ( headmodellc.label, label ) );
    
    if isempty ( hindex )
        continue
    end
    
    % Creats the figure.
    figure
    hold on
    
    % Plots all the sensors.
    sensors = sens.chanpos;
    plot3 ( sensors ( :, 1 ), sensors ( :, 2 ), sensors ( :, 3 ), 'k.', 'MarkerSize', 1 )
    
    % Plots the sensor.
    sensor = sens.chanpos ( hindex, : );
    plot3 ( sensor (1), sensor (2), sensor (3), 'r*' )
    
    % Plots the three meshes.
    ft_plot_mesh ( mridata.mesh.bnd (1), 'edgecolor', 'none', 'facecolor', 'brain', 'facealpha', 0.1 );
    ft_plot_mesh ( mridata.mesh.bnd (2), 'edgecolor', 'none', 'facecolor', 'white', 'facealpha', 0.1 );
    ft_plot_mesh ( mridata.mesh.bnd (3), 'edgecolor', 'none', 'facecolor', 'skin',  'facealpha', 0.3 );
    
    % Gets the leadfield for this sensor.
    lf = reshape ( leadfieldY ( :, sindex ), 3, [] );
    lf = sqrt ( sum ( lf .^ 2, 1 ) );
    lf = lf - min ( lf );
    lf = lf / max ( lf );
    
    color = bsxfun ( @times, [ 0 1 1 ], 1 - lf' );
    color = bsxfun ( @plus, [ 1 0 0 ], color );
    scatter3 ( mridata.grid.pos ( :, 1 ), mridata.grid.pos ( :, 2 ), mridata.grid.pos ( :, 3 ), 15, color, 'filled' )
    
    
    views = linspace ( 0, 360, 100 );
    views = views ( 1: end - 1 );
    for vindex = 1: numel ( views )
        view ( views ( vindex ), 0 );
        pause ( 0.04 );
    end
    close
    
%     savefig ( sprintf ( '%s.fig', sens.label { sindex } ) );
%     my_savegif ( sprintf ( '%s.gif', sens.label { sindex } ) );
%     close
    
%     rotate3d
%     uiwait
end
%%

% For each source, draws the effect in each sensor.

% If EEG projects the channels onto the outter mesh.
if isfield ( sens, 'elecpos'  )
    for eindex = 1: 60
        [ ~, Pm ] = NFT_dmp ( sens.elecpos ( eindex, : ), mridata.mesh.bnd ( end ).pnt, mridata.mesh.bnd ( end ).tri );
        sens.elecpos ( eindex, : ) = Pm;
    end
    for cindex = 1: 60
        [ ~, Pm ] = NFT_dmp ( sens.chanpos ( cindex, : ), mridata.mesh.bnd ( end ).pnt, mridata.mesh.bnd ( end ).tri );
        sens.chanpos ( cindex, : ) = Pm;
    end
end

for dindex = 1: 100: size ( mridata.grid.pos, 1 )
    
    label  = sens.label ( 1: 60 );
    hindex = find ( strcmp ( headmodellc.label, label ) );
    
    % Creats the figure.
    figure
    hold on
    
    % Plots the three meshes.
    ft_plot_mesh ( mridata.mesh.bnd (1), 'edgecolor', 'none', 'facecolor', 'brain', 'facealpha', 0.5 );
    ft_plot_mesh ( mridata.mesh.bnd (2), 'edgecolor', 'none', 'facecolor', 'white', 'facealpha', 0.3 );
    ft_plot_mesh ( mridata.mesh.bnd (3), 'edgecolor', 'none', 'facecolor', 'skin',  'facealpha', 0.3 );
    
    % Plots all the soruces.
    plot3 ( mridata.grid.pos ( :, 1 ), mridata.grid.pos ( :, 2 ), mridata.grid.pos ( :, 3 ), '.k', 'MarkerSize', 1 )
    
    % Plots the current source.
    plot3 ( mridata.grid.pos ( dindex, 1 ), mridata.grid.pos ( dindex, 2 ), mridata.grid.pos ( dindex, 3 ), 'b*' )
    
    
    % Gets the leadfield for this source.
    lf = reshape ( leadfieldY', 60, 3, [] );
    lf = sqrt ( sum ( lf ( :, :, dindex ) .^ 2, 2 ) );
    lf = lf - min ( lf );
    lf = lf / max ( lf );
    
    % Plots all the sensors.
    color = bsxfun ( @times, [ 0 1 1 ], 1 - lf );
    color = bsxfun ( @plus, [ 1 0 0 ], color );
    scatter3 ( sens.chanpos ( hindex, 1 ), sens.chanpos ( hindex, 2 ), sens.chanpos ( hindex, 3 ), 15, color, 'filled' )
    
    
    views = linspace ( 0, 360, 100 );
    views = views ( 1: end - 1 );
    for vindex = 1: numel ( views )
        view ( views ( vindex ), 0 );
        pause ( 0.04 );
    end
    close
    
    
%     savefig ( sprintf ( '%s.fig', sens.label { sindex } ) );
%     my_savegif ( sprintf ( '%s.gif', sens.label { sindex } ) );
%     close
    
%     rotate3d
%     uiwait
end
