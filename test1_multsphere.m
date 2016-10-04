
clc
clear
close all

% Adds FT to the path.
ft_path;
ft_defaults

ft_hastoolbox ( 'freesurfer', 1, 1 );
ft_hastoolbox ( 'spm8', 1, 1 );

% Adds the functions folder to the path.
addpath ( sprintf ( '%s/functions', pwd ) );

% Loads the EEG data.
eegdata = load ( 'alc02_restingOA_EEG' );
eegdata.trialdata.elec = ft_convert_units ( eegdata.trialdata.elec, 'm' );
eegdata.trialdata.grad = ft_convert_units ( eegdata.trialdata.grad, 'm' );

% Loads the subject segmentation and the generated meshes.
mridata = load ( 'alc02_mri.mat' );
mridata.mesh = ft_convert_units ( mridata.mesh, 'mm' );
mridata.mesh = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.mesh );
mridata.mesh = ft_convert_units ( mridata.mesh, 'm' );

sens = eegdata.trialdata.elec;
mesh = mridata.mesh;

% headmodel = my_headmodel_eegspheres ( mesh, sens, 'singlesphere', 'yes' );
headmodel = my_headmodel_eegspheres ( mesh, sens );





% Only keeps the sensors with a defined position.
if isfield ( sens, 'coilpos'  )
    nocoil = any ( isnan ( sens.coilpos ), 2 );
    if isfield ( sens, 'coilpos'  ), sens.coilpos  ( nocoil, : ) = []; end
    if isfield ( sens, 'coilori'  ), sens.coilori  ( nocoil, : ) = []; end
    if isfield ( sens, 'tra'      ), sens.tra      ( :, nocoil ) = []; end
end
if isfield ( sens, 'elecpos'  )
    noelec = any ( isnan ( sens.elecpos ), 2 );
    if isfield ( sens, 'elecpos'  ), sens.elecpos  ( noelec, : ) = []; end
    if isfield ( sens, 'tra'      ), sens.tra      ( :, noelec ) = []; end
end
if isfield ( sens, 'chanpos'  )
    nochan = any ( isnan ( sens.chanpos ), 2 );
    if isfield ( sens, 'chanpos'  ), sens.chanpos  ( nochan, : ) = []; end
    if isfield ( sens, 'chantype' ), sens.chantype ( nochan, : ) = []; end
    if isfield ( sens, 'chanunit' ), sens.chanunit ( nochan, : ) = []; end
    if isfield ( sens, 'label'    ), sens.label    ( nochan, : ) = []; end
    if isfield ( sens, 'tra'      ), sens.tra      ( nochan, : ) = []; end
end

% If EEG projects the channels onto the outter mesh.
if isfield ( sens, 'elecpos'  )
    for eindex = 1: size ( sens.elecpos, 1 );
        [ ~, Pm ] = NFT_dmp ( sens.elecpos ( eindex, : ), mesh.bnd ( end ).pnt, mesh.bnd ( end ).tri );
        sens.elecpos ( eindex, : ) = Pm;
    end
    for cindex = 1: size ( sens.chanpos, 1 );
        [ ~, Pm ] = NFT_dmp ( sens.chanpos ( cindex, : ), mesh.bnd ( end ).pnt, mesh.bnd ( end ).tri );
        sens.chanpos ( cindex, : ) = Pm;
    end
end



%%

for sindex = 1: size ( sens.chanpos, 1 )
    
    % Creats the figure.
    figure
    hold on
    
    % Plots all the sensors.
    sensors = sens.chanpos;
    plot3 ( sensors ( :, 1 ), sensors ( :, 2 ), sensors ( :, 3 ), 'k.' )
    
    % Plots the sensor.
    sensor = sens.chanpos ( sindex, : );
    plot3 ( sensor (1), sensor (2), sensor (3), 'r*' )
    
    % Plots the three meshes.
    ft_plot_mesh ( mesh.bnd (1), 'edgecolor', 'none', 'facecolor', 'brain', 'facealpha', 0.3 );
    ft_plot_mesh ( mesh.bnd (2), 'edgecolor', 'none', 'facecolor', 'white', 'facealpha', 0.3 );
    ft_plot_mesh ( mesh.bnd (3), 'edgecolor', 'none', 'facecolor', 'skin',  'facealpha', 0.3 );
    
    % Plots the three spheres.
    spheres      = headmodel;
    spheres.o    = headmodel.o ( sindex, : );
    spheres.type = 'concentricspheres';
    spheres.r    = headmodel.r ( sindex, 1 );
    ft_plot_vol ( spheres, 'edgecolor', 'black', 'facecolor', 'brain', 'facealpha', 0.8 )
    spheres.r    = headmodel.r ( sindex, 2 );
    ft_plot_vol ( spheres, 'edgecolor', 'none',  'facecolor', 'white', 'facealpha', 0.6 )
    spheres.r    = headmodel.r ( sindex, 3 );
    ft_plot_vol ( spheres, 'edgecolor', 'none',  'facecolor', 'skin',  'facealpha', 0.3 )
    
    rotate3d
%     savefig ( sprintf ( '%s.fig', sens.label { sindex } ) );
%     my_savegif ( sprintf ( '%s.gif', sens.label { sindex } ) );
%     close
    uiwait
end
