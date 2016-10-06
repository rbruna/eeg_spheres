function grid = my_leadfield ( cfg, data )

% Leadfield calculation.
%
% Use as
%   grid = ft_prepare_leadfield ( cfg, data );
%
% where:
%   cfg   Structure as provided to ft_prepare_leadfield.
%   data  FieldTrip data to select channels (optional).
%
% This function mimics the FieldTrip function ft_prepare_leadfield. See
% help on this function for more information.
%
% This function requires FieldTrip 20160222 or newer to work properly.

% Undocumented local options:
% cfg.sel50p      = 'no' (default) or 'yes'
% cfg.lbex        = 'no' (default) or a number that corresponds with the radius
% cfg.mollify     = 'no' (default) or a number that corresponds with the FWHM

% Based on FieldTrip 20160222 functions:
% * ft_prepare_leadfield by Robert Oostenveld
% * ft_compute_leadfield by Robert Oostenveld


% These options are not yet implemented.
if any ( isfield ( cfg, { 'mollify' 'patchsvd' 'sel50p' 'lbex' } ) )
    warning ( 'ft_prepare_leadfield''s options ''mollify'', ''patchsvd'', ''sel50p'' and ''lbex'' are not yet implemented. Ignoring.' );
end

% % Checks that a maximum of one subspace projection is selected.
% if isfield ( cfg, 'sel50p' ) && isfield ( cfg, 'lbex' )
%     error ( 'Only one subspace projection method is allowed.' );
% end


% Checks the data, if provided.
if nargin > 1
    data = ft_checkdata ( data, 'feedback', 'no' );
else
    data = [];
end

% Gets the sensor definition.
if isfield ( cfg, 'sens' )
    
    % Gets the sensor definition from the configuration.
    sens = cfg.sens;
else
    
    error ( 'No sensors'' definition provided.' );
%     
%     % Gets the sensor definition from the data.
%     ft_fetch_sens ( cfg, 'data' );
end

% Gets the grid definition.
if isfield ( cfg, 'grid' )
    
    % Gets the sensor definition from the configuration.
    grid = cfg.grid;
else
    
    error ( 'No grid''s definition provided.' );
    
    % Creates a grid definition.
%     tmpcfg           = keepfields(cfg, {'grid', 'mri', 'headshape', 'symmetry', 'smooth', 'threshold', 'spheremesh', 'inwardshift'});
%     tmpcfg.headmodel = headmodel;
%     tmpcfg.grad      = sens; % either electrodes or gradiometers
%     grid = ft_prepare_sourcemodel(tmpcfg);
end

% Gets the head model.
if isfield ( cfg, 'headmodel' )
    
    % Gets the sensor definition from the configuration.
    headmodel = cfg.headmodel;
else
    
    error ( 'No head model provided.' );
%     
%     % Calculates the headmodel.
end

% Sanitizes the grid.
grid = fixgrid ( grid );


% Determines if the data is EEG or MEG.
ismeg = isfield ( sens, 'coilpos' );
iseeg = isfield ( sens, 'elecpos' );


% Gets the channels present in the sensor definition.
channels = sens.label;

% If a channel restriction has been provided uses it.
if isfield ( cfg, 'channel' )
    channels = ft_channelselection ( cfg.channel, channels );
end

% If data, restricts the channels to that.
if data
    channels = ft_channelselection ( data.label, channels );
end

% If the head model is channel-dependent get only those channels.
if isfield ( headmodel, 'label' )
    channels = ft_channelselection ( headmodel.label, channels );
end

% Removes the unused channels and sensors from the sensor definition.
sens = fixsens ( sens, channels );


% If channel units checks if the unit transformation is possible.
if isfield ( cfg, 'chanunit' )
    
    % If no channel units the transformation is not possible.
    if ~isfield ( sens, 'chanunit' )
        
        % Removes the channel units.
        cfg = rmfield ( cfg, 'chanunit' );
        
    % Transformation onlo possible for SI units.
    elseif ~all ( ismember ( sens.chanunit, { 'V' 'V/m' 'T' 'T/m' } ) )
        warning ( 'Channel units only can be adjusted if sensors are defined in SI units. Ignoring.' );
        
        % Removes the channel units.
        cfg = rmfield ( cfg, 'chanunit' );
    end
end

% If channel or sources units transforms everything to SI units (m).
if any ( isfield ( cfg, { 'chanunit' 'dipoleunit' } ) )
    
    % Makes sure that both the grid and the sensor definition have units.
    grid = ft_convert_units ( grid, 'm' );
    sens = ft_convert_units ( sens, 'm' );
    
    % Checks if the units are SI units.
    if ~strcmp ( grid.unit, 'm' )
        warning ( 'Converting the grid to meters in order to adjust the channel and dipole units.' );
        grid = ft_convert_units ( grid, 'm' );
    end
    if ~strcmp ( sens.unit, 'm' )
        warning ( 'Converting the sensor definition to meters in order to adjust the channel and dipole units.' );
        sens = ft_convert_units ( sens, 'm' );
    end
    
    % Initializes the leadfield scale.
    lfscale = 1;
    
    % Adjust the scale for the sensors' edge.
    if isfield ( cfg, 'chanunit' )
        lfscale = lfscale * cellfun ( @ft_scalingfactor, sens.chanunit (:), cfg.chanunit (:) );
    end
    
    % Adjusts the scale in the dipoles' edge.
    if isfield ( cfg, 'dipoleunit' )
        lfscale = lfscale / ft_scalingfactor ( 'A*m', cfg.dipoleunit );
    end
    
else
    lfscale = 1;
end


% Gets the number of channels and sources.
nchannels = numel ( channels );
nsources  = sum   ( grid.inside, 1 );

% Selects the right function for each head model.
switch headmodel.type
    
    % Local concentric spheres.
    case 'localconcentricspheres'
        
        % Reserves memory for the leadfield matrix in 3D form.
        leadfield = zeros ( nchannels, 3, nsources );
        
        % Iterates along channels.
        for cindex = 1: nchannels
            
            % Gets the sensor label.
            senslabel      = sens.label { cindex };
            
            % Creates a dummy head model containing only this channel.
            sensmodel      = [];
            sensmodel.r    = headmodel.r ( strcmp ( headmodel.label, senslabel ), : );
            sensmodel.o    = headmodel.o ( strcmp ( headmodel.label, senslabel ), : );
            sensmodel.unit = headmodel.unit;
            sensmodel.cond = headmodel.cond;
            sensmodel.type = 'concentricspheres';
            
            % Calculates the leadfield for the current channel.
            leadfield ( cindex, : ) = my_leadfield_eegspheres (  grid.pos ( grid.inside, : ), sens.elecpos ( cindex, : ), sensmodel );
        end
        
    % Concentric spheres.
    case 'concentricspheres'
        
        % Calculates the leadfield.
        leadfield = my_leadfield_eegspheres ( grid.pos ( grid.inside, : ), sens.elecpos, headmodel );
        
        % Rewrites the leadfield as a 3D matrix.
        leadfield = reshape ( leadfield, nchannels, 3, nsources );
        
    % Otherwise relies on FieldTrip.
    otherwise
        
        % Calculates the leadfield using FieldTrip.
        grid      = ft_prepare_leadfield ( cfg, data );
        return
end


% If 'tra' field compose the channel leadfield from the sensors.
if isfield ( sens, 'tra' )
    leadfield = sens.tra * leadfield ( :, : );
    leadfield = reshape ( leadfield, [], 3, nsources );
    
% If no 'tra' field and EEG sensors assumes average reference.
elseif iseeg
    leadfield = bsxfun ( @minus, leadfield, mean ( leadfield, 1 ) );
end


% Determines if apply rank reduction or not.
if ~isfield ( cfg, 'reducerank' ) && ismeg
        cfg.reducerank = 2;
end

% Applies the rank reduction.
if isfield ( cfg, 'reducerank' ) && cfg.reducerank < 3
    
    % Goes through each dipole.
    for sindex = 1: nsources
        
        % Performs a SVD over the data.
        [ u, s, v ] = svd ( leadfield ( :, :, sindex ) );
        
        % Removes the smaller singular values from the matrix V.
        v ( :, cfg.reducerank + 1: end ) = 0;
        
        % Recomposes the leadfield.
        if ~isfield ( cfg, 'backproject' ) || cfg.backproject
            leadfield ( :, :, sindex ) = u * s * v';
            
        % Removes the smaller singular values from the leadfield.
        else
            leadfield ( :, :, sindex ) = leadfield ( :, :, sindex ) * v;
        end
    end
    
    % Deletes the null dimensions of the leadfield matrix.
    if isfield ( cfg, 'backproject' ) && ~cfg.backproject
        leadfield ( :, cfg.reducerank + 1: end, : ) = [];
    end
end


% Normalizes the leadfield, if requested.
if isfield ( cfg, 'normalize' ) 
    switch cfg.normalize
        case 'yes'
            norm = sum ( sum ( leadfield .^ 2, 1 ), 2 ) .^ cfg.normalizeparam;
        case 'column'
            norm = sum ( leadfield .^ 2, 1 ) .^ cfg.normalizeparam;
        otherwise
            norm = 1;
    end
    
    % Normalizes the leadfield.
    leadfield = bsxfun ( @rdivide, leadfield, norm );
end


% Applies a weight to each dipole, if requested.
if isfield ( cfg, 'weight' ) && numel ( cfg.weight ) == nsources
    leadfield = bsxfun ( @times, leadfield, reshape ( cfg.weight, 1, 1, [] ) );
end


% Scales the leadfield, if requested.
leadfield = bsxfun ( @times, leadfield, lfscale );

% Projects the dipoles over the dipole moment, if requested.
% (Not sure about the dimensions).
if isfield ( grid, 'mom' )
    leadfield = bsxfun ( @times, leadfield, grid.mom );
end


% Generates the leadfield cell.
grid.leadfield = cell ( 1, size ( grid.pos, 1 ) );

% Stores the leadfield in the grid.
grid.leadfield ( grid.inside ) = num2cell ( leadfield, [ 1 2 ] );
grid.label     = channels;
grid.leadfielddimord = '{pos}_chan_ori';


% % mollify the leadfields
% if ~strcmp(cfg.mollify, 'no')
%   grid = mollify(cfg, grid);
% end
% 
% % combine leadfields in patches and do an SVD on them
% if ~strcmp(cfg.patchsvd, 'no')
%   grid = patchsvd(cfg, grid);
% end
% 
% % Keeps only the nearest 50% of the channels.
% sens.coilpos=sens.elecpos;
% if isfield ( cfg, 'sel50p' )
%   grid = sel50p ( cfg, grid, sens );
% end

% % compute the local basis function expansion (LBEX) subspace projection
% if ~strcmp(cfg.lbex, 'no')
%   grid = lbex(cfg, grid);
% end


function grid = fixgrid ( grid )

% Fix the absent or numerical 'inside' field.
if ~isfield ( grid, 'inside' ) || islogical ( grid.inside )
    
    % Initializes the 'inside' field to a logical array.
    inside = false ( size ( grid.pos ( :, 1 ) ) );
    
    % If numerical 'inside' field converts it to logical.
    if isfield ( grid, 'inside' )
        inside ( grid.inside ) = true;
        
    % Otherwise all the points are inside.
    else
        inside (:) = true;
    end
    
    % Stores the 'inside' field.
    grid.inside = inside;
end


function sens = fixsens ( sens, channels )

% Modifies the sensors definition, if needed.
if ~all ( ismember ( sens.label, channels ) )
    
    % Gets the indexes of the channels to remove.
    remidx = ~ismember ( sens.label, channels );
    
    % Removes the sensors from the 'label' field.
    sens.label ( remidx ) = [];
    
    % Removes the channels from each channel field.
    if isfield ( sens, 'chanpos' ),  sens.chanpos  ( remidx, : ) = []; end
    if isfield ( sens, 'chanori' ),  sens.chanori  ( remidx, : ) = []; end
    if isfield ( sens, 'chantype' ), sens.chantype ( remidx ) = []; end
    if isfield ( sens, 'chanunit' ), sens.chanunit ( remidx ) = []; end
    
    % If 'tra' removes the channels and the associated sensors.
    if isfield ( sens, 'tra' )
        
        % Removes the channels.
        sens.tra ( remidx, : ) = [];
        
        % Gets the list of useless sensors.
        remidx = ~any ( sens.tra, 1 );
        
        % Removes the sensors from thwe 'tra' field.
        sens.tra ( :, remidx ) = [];
        
        % Removes the sensors from each sensor field.
        if isfield ( sens, 'coilpos' ),  sens.coilpos  ( remidx, : ) = []; end
        if isfield ( sens, 'coilori' ),  sens.coilori  ( remidx, : ) = []; end
        if isfield ( sens, 'elecpos' ),  sens.elecpos  ( remidx, : ) = []; end
    end
    
    % If no 'tra' field assumes that the mapping is 1:1.
    if ~isfield ( sens, 'tra' )
        
        % Removes the sensors from each sensor field.
        if isfield ( sens, 'coilpos' ),  sens.coilpos  ( remidx, : ) = []; end
        if isfield ( sens, 'coilori' ),  sens.coilori  ( remidx, : ) = []; end
        if isfield ( sens, 'elecpos' ),  sens.elecpos  ( remidx, : ) = []; end
    end
end
