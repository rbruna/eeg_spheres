function grid = my_leadfield ( cfg, data )

% Leadfield calculation.
%
% This function requires FieldTrip 20160222 or newer to work properly.

% FT_PREPARE_LEADFIELD computes the forward model for many dipole locations
% on a regular 2D or 3D grid and stores it for efficient inverse modelling
%
% Use as
%   [grid] = ft_prepare_leadfield(cfg, data);
%
% It is neccessary to input the data on which you want to perform the
% inverse computations, since that data generally contain the gradiometer
% information and information about the channels that should be included in
% the forward model computation. The data structure can be either obtained
% from FT_PREPROCESSING, FT_FREQANALYSIS or FT_TIMELOCKANALYSIS. If the data is empty,
% all channels will be included in the forward model.
%
% The configuration should contain
%   cfg.channel            = Nx1 cell-array with selection of channels (default = 'all'),
%                            see FT_CHANNELSELECTION for details
%
% The positions of the sources can be specified as a regular 3-D
% grid that is aligned with the axes of the head coordinate system
%   cfg.grid.xgrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.ygrid      = vector (e.g. -20:1:20) or 'auto' (default = 'auto')
%   cfg.grid.zgrid      = vector (e.g.   0:1:20) or 'auto' (default = 'auto')
%   cfg.grid.resolution = number (e.g. 1 cm) for automatic grid generation
% Alternatively the position of a few sources at locations of interest can
% be specified, for example obtained from an anatomical or functional MRI
%   cfg.grid.pos        = N*3 matrix with position of each source
%   cfg.grid.inside     = N*1 vector with boolean value whether grid point is inside brain (optional)
%   cfg.grid.dim        = [Nx Ny Nz] vector with dimensions in case of 3-D grid (optional)
%
% The volume conduction model of the head should be specified as
%   cfg.headmodel     = structure with volume conduction model, see FT_PREPARE_HEADMODEL
%
% The EEG or MEG sensor positions can be present in the data or can be specified as
%   cfg.elec          = structure with electrode positions, see FT_DATATYPE_SENS
%   cfg.grad          = structure with gradiometer definition, see FT_DATATYPE_SENS
%   cfg.elecfile      = name of file containing the electrode positions, see FT_READ_SENS
%   cfg.gradfile      = name of file containing the gradiometer definition, see FT_READ_SENS
%
% Optionally, you can modify the leadfields by reducing the rank (i.e.
% remove the weakest orientation), or by normalizing each column.
%   cfg.reducerank      = 'no', or number (default = 3 for EEG, 2 for MEG)
%   cfg.normalize       = 'yes' or 'no' (default = 'no')
%   cfg.normalizeparam  = depth normalization parameter (default = 0.5)
%   cfg.backproject     = 'yes' or 'no' (default = 'yes') determines when reducerank is applied
%                         whether the lower rank leadfield is projected back onto the original
%                         linear subspace, or not.
%
% To facilitate data-handling and distributed computing you can use
%   cfg.inputfile   =  ...
% If you specify this option the input data will be read from a *.mat
% file on disk. This mat files should contain only a single variable named 'data',
% corresponding to the input structure.
%
% See also FT_SOURCEANALYSIS, FT_DIPOLEFITTING, FT_PREPARE_HEADMODEL,
% FT_PREPARE_SOURCEMODEL

% Undocumented local options:
% cfg.feedback
% cfg.sel50p      = 'no' (default) or 'yes'
% cfg.lbex        = 'no' (default) or a number that corresponds with the radius
% cfg.mollify     = 'no' (default) or a number that corresponds with the FWHM

% Based on FieldTrip 20160222 functions:
% * ft_prepare_leadfield by Robert Oostenveld
% * ft_compute_leadfield by Robert Oostenveld


% % set the defaults
% cfg.normalize      = ft_getopt(cfg, 'normalize',      'no');
% cfg.normalizeparam = ft_getopt(cfg, 'normalizeparam', 0.5);
% cfg.lbex           = ft_getopt(cfg, 'lbex',           'no');
% cfg.sel50p         = ft_getopt(cfg, 'sel50p',         'no');
% cfg.feedback       = ft_getopt(cfg, 'feedback',       'text');
% cfg.mollify        = ft_getopt(cfg, 'mollify',        'no');
% cfg.patchsvd       = ft_getopt(cfg, 'patchsvd',       'no');
% cfg.backproject    = ft_getopt(cfg, 'backproject',    'yes'); % determines whether after rank reduction the subspace projected leadfield is backprojected onto the original space
% % cfg.reducerank   = ft_getopt(cfg, 'reducerank', 'no');      % the default for this depends on EEG/MEG and is set below
% 
% % put the low-level options pertaining to the dipole grid in their own field
% cfg = ft_checkconfig(cfg, 'renamed', {'tightgrid', 'tight'});  % this is moved to cfg.grid.tight by the subsequent createsubcfg
% cfg = ft_checkconfig(cfg, 'renamed', {'sourceunits', 'unit'}); % this is moved to cfg.grid.unit by the subsequent createsubcfg
% cfg = ft_checkconfig(cfg, 'createsubcfg',  {'grid'});
% 
% % this code expects the inside to be represented as a logical array
% cfg.grid = ft_checkconfig(cfg.grid, 'renamed',  {'pnt' 'pos'});
% cfg = ft_checkconfig(cfg, 'index2logical', 'yes');
% 
% if strcmp(cfg.sel50p, 'yes') && strcmp(cfg.lbex, 'yes')
%   error('subspace projection with either lbex or sel50p is mutually exclusive');
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







% % set the default for reducing the rank of the leadfields
% if ft_senstype(sens, 'eeg')
%   cfg.reducerank = ft_getopt(cfg, 'reducerank', 3);
% else
%   cfg.reducerank = ft_getopt(cfg, 'reducerank', 2);
% end






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
            senslabel = sens.label { cindex };
            
            % Creates a dummy head model containing only this channel.
            sensmodel = headmodel;
            sensmodel.r = headmodel.r ( strcmp ( headmodel.label, senslabel ), : );
            sensmodel.o = headmodel.o ( strcmp ( headmodel.label, senslabel ), : );
            
            % Calculates the leadfield for the current channel.
            leadfield ( cindex, :, : ) = my_leadfield_eegspheres (  grid.pos ( grid.inside, : ), sens.elecpos ( cindex ), sensmodel );
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

% if ~isempty(chanunit) || ~isempty(dipoleunit)
%   assert(strcmp(headmodel.unit,  'm'), 'unit conversion only possible for SI input units');
%   assert(strcmp(sens.unit, 'm'), 'unit conversion only possible for SI input units');
% end
% 
% if ~isempty(chanunit)
%   assert(all(strcmp(sens.chanunit, 'V') | strcmp(sens.chanunit, 'V/m') | strcmp(sens.chanunit, 'T') | strcmp(sens.chanunit, 'T/m')), 'unit conversion only possible for SI input units');
%   % compute conversion factor and multiply each row of the matrix
%   scale = cellfun(@ft_scalingfactor, sens.chanunit(:), chanunit(:));
%   lf = bsxfun(@times, lf, scale(:));
%   % prior to this conversion, the units might be  (T/m)/(A*m) for planar gradients or   (V/m)/(A*m) for bipolar EEG
%   % after this conversion, the units will be     (T/cm)/(A*m)                      or (uV/mm)/(A*m)
% end
% 
% if ~isempty(dipoleunit)
%   scale = ft_scalingfactor('A*m', dipoleunit); % compue the scaling factor from A*m to the desired dipoleunit
%   lf    = lf/scale;                         % the leadfield is expressed in chanunit per dipoleunit, i.e. chanunit/dipoleunit
% end






% if isfield(cfg, 'grid') && isfield(cfg.grid, 'mom')
%     % multiply with the normalized dipole moment to get the leadfield in the desired orientation
%     grid.leadfield{thisindx} = grid.leadfield{thisindx} * grid.mom(:,thisindx);
% end




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
% % compute the 50 percent channel selection subspace projection
% if ~strcmp(cfg.sel50p, 'no')
%   grid = sel50p(cfg, grid, sens);
% end
% 
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
