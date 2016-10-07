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
