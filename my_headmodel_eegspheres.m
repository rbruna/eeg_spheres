function headmodel = my_headmodel_eegspheres ( mesh, sens, varargin )

% Creates a concentric spheres volume conductor for EEG.
%
% Use as
%   headmodel = my_headmodel_eegspheres ( mesh, grad, ... )
%
% Optional arguments should be specified in key-value pairs and can include
%   radius    = number, radius of sphere within which headshape points will
%               be included for the fitting algorithm
%   maxradius = number, if for a given sensor the fitted radius exceeds
%               this value, the radius and origin will be replaced with the
%               single sphere fit
%   baseline  = number
%
% This function requires FieldTrip 20160222 or newer to work properly.

% Based on FieldTrip 20160222 functions:
% * ft_headmodel_concentricspheres
% * ft_headmodel_localspheres

% Copyright (C) 2016, Ricardo Bruna

% Defines the default conductivities.
conductivities = [
    0.3300, ... % Brain conductivity.
    1.0000, ... % CSF conductivity.
    0.0042, ... % Skull conductivity.
    0.3300 ];   % Skin conductivity.

% Defines the minimum skull and scalp thicknesses.
minskull = 0.005;
minscalp = 0.007;


% Sanitices the meshes.
mesh = fixmesh ( mesh );

% Sorts the meshes from inner to outter.
mesh = sortmesh ( mesh );

% Makes sure that mesh and sensors have the same geometrical units.
sens = ft_convert_units ( sens, mesh (1).unit );

% Transforms the minimum thicknesses to the meshes units.
minskull = minskull * ft_scalingfactor ( 'm', mesh (1).unit );
minscalp = minscalp * ft_scalingfactor ( 'm', mesh (1).unit );


% Determines if the input is EEG or MEG.
ismeg = isfield ( sens, 'coilpos' );
iseeg = isfield ( sens, 'elecpos' );

% If more than one mesh, the head model is designed for EEG.
if ismeg && numel ( mesh ) > 1
    error ( 'Local concentric spheres only works with EEG' );
end


% Defines the conductivities depending on the tissues.
if numel ( mesh ) == 1
    conductivities = 1;
elseif numel ( mesh ) == 3
    conductivities = conductivities ( [ 1 3 4 ] );
elseif numel ( mesh ) == 4
    conductivities = conductivities ( [ 1 2 3 4 ] );
else
    error ( 'Wrong number of meshes.' );
end


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
if iseeg
    for eindex = 1: size ( sens.elecpos, 1 );
        [ ~, Pm ] = NFT_dmp ( sens.elecpos ( eindex, : ), mesh ( end ).pos, mesh ( end ).tri );
        sens.elecpos ( eindex, : ) = Pm;
    end
    for cindex = 1: size ( sens.chanpos, 1 );
        [ ~, Pm ] = NFT_dmp ( sens.chanpos ( cindex, : ), mesh ( end ).pos, mesh ( end ).tri );
        sens.chanpos ( cindex, : ) = Pm;
    end
end



% Sets the default parameters in the desired untis.
threshold = ft_getopt(varargin, 'radius',    ft_scalingfactor('cm', mesh (1).unit) * 8.5);
maxradius = ft_getopt(varargin, 'maxradius', ft_scalingfactor('cm', mesh (1).unit) * 20);
baseline  = ft_getopt(varargin, 'baseline',  ft_scalingfactor('cm', mesh (1).unit) * 5);

% get the additional inputs and set the defaults
singlesphere  = ft_getopt(varargin, 'singlesphere', 'no');



% Initializes the center and radius arrays.
nchans  = numel ( sens.label );
radii   = zeros ( nchans, numel ( mesh ) );
centers = zeros ( nchans, 3 );


% Suppress the 'nearly singlular matrix' warking.
w = warning ( 'off', 'MATLAB:nearlySingularMatrix' );

% Fits a single sphere to all the points over z = 4 cm.
pos      = unique ( cat ( 1, mesh.pos ), 'rows');
sphereps = pos ( pos ( :, 3 ) > 0.04, : );
[ center, radius ] = my_fitsphere ( sphereps );

% If only one sphere is required uses this as its center.
if strcmp ( singlesphere, 'yes' )
    
    % Reserves memory for the radius.
    radius = zeros ( 1, numel ( mesh ) );
    
    % Gets the radii for each sphere.
    for mindex = 1: numel ( mesh )
        
        % Calculates the radius as the average distance to the center.
        pos      = mesh ( mindex ).pos;
        sphereps = pos ( pos ( :, 3 ) > 0.04, : );
        distance = sqrt ( sum ( bsxfun ( @minus, sphereps, center ) .^ 2, 2 ) );
        radius ( mindex ) = mean ( distance );
    end
    
    % Sets the output.
    headmodel       = [];
    headmodel.type  = 'concentricspheres';
    headmodel.o     = center;
    headmodel.r     = radius;
    headmodel.unit  = mesh (1).unit;
    headmodel.cond  = conductivities;
    
    % Restores the original warnings.
    warning ( w )
    
    return
end

for chan = 1: nchans
    
    % In MEG looks for the coil nearer to the head (?).
    if ismeg
        
        % Gets all the coils of the current channel.
        coilsel = find ( sens.tra ( chan, : ) );
        allpos  = sens.coilpos ( coilsel, : );
        allori  = sens.coilori ( coilsel, : );
        
        % Gets the average position of the coils.
        chanpos = mean ( allpos, 1 );
        
        % Gets the average orientation of the coils.
        [ ~, ~, v ] = svd ( allori );
        chanori = v ( :, 1 )';
        
        % Ensures that the orientation points outwards.
        chanori = chanori * sign ( dot ( chanpos, chanori ) );
        
        % Computes the distance from every coil along that orientation.
        rawdist  = bsxfun ( @minus, allpos, chanpos );
        distance = sum ( bsxfun ( @times, rawdist, chanori ), 2 );
        
        % Checs the distance from the point to all the coils.
        [ m, i ] = min ( distance );
        
        % If too far selects the nearest coil as center of the coils.
        if abs ( m ) > ( baseline / 4 )
            chanpos = allpos ( i, : );
        end
        
    % In EEG gets the position of the electrode.
    else
        chanpos = sens.chanpos ( chan, : );
    end
    
    
    % Find the first mesh's points closer to the sensor.
    distance = sqrt ( sum ( bsxfun ( @minus, mesh (1).pos, chanpos ) .^ 2, 2 ) );
    sphereps = mesh (1).pos ( distance < threshold, : );
    
    % Fits a sphere to those points.
    if size ( sphereps, 1 ) > 10
        [ o, r ] = my_fitsphere ( sphereps );
        
        % Forces the radius to be the same to all the spheres.
        
    else
        fprintf ( 'Not enough surface points for channel %s. Using all points.\n', sens.label { chan } );
        o = center;
        r = radius;
    end
    
    if r > maxradius
        fprintf ( 'Sphere for channel %s too large. Using all points.\n', sens.label { chan } );
        o = center;
        r = radius;
    end
    
    % If several meshes calculares the radii for the other spheres.
    if numel ( mesh ) > 1
        
        % Extends the radii vector.
        r ( numel ( mesh ) ) = 0;
        
        % Goes through each mesh.
        for mindex = 2: numel ( mesh )
            
            % Find the first mesh's points closer to the sensor.
            distance = sqrt ( sum ( bsxfun ( @minus, mesh ( mindex ).pos, chanpos ) .^ 2, 2 ) );
            sphereps = mesh ( mindex ).pos ( distance < threshold, : );
            
            % Calculates the radius as the mean distance to the center.
            distance = sqrt ( sum ( bsxfun ( @minus, sphereps, o ) .^ 2, 2 ) );
            r ( mindex ) = mean ( distance );
        end
    end
    
    % If EEG makes sure that the sensor is on the outter surface.
    if iseeg
        r ( end ) = sqrt ( sum ( ( o - chanpos ) .^ 2 ) );
    end
    
    % Stores the sphere center and radius.
    centers ( chan, : ) = o;
    radii   ( chan, : ) = r;
end

% Gets the thickness of the skull and the scalp.
dradii = diff ( radii, 1, 2 );

% Fixes the thickness if too low.
if any ( dradii ( :, 1 ) < minskull )
    fprintf ( 1, 'The skull is too thin in some points. Trying to fix them.\n' );
    
    % Sets the minimum skull thickness to 5mm.
    radii ( :, 2 ) = max ( radii ( :, 1 ) + minskull, radii ( :, 2 ) );
end
if any ( dradii ( :, 2 ) < minscalp )
    fprintf ( 1, 'The scalp is too thin in some points. Trying to fix them.\n' );
    
    % Sets the minimum scalp thickness to 7mm.
    radii ( :, 2 ) = max ( radii ( :, 1 ) + minscalp, radii ( :, 2 ) );
end


% Sets the output.
headmodel       = [];
headmodel.type  = 'localconcentricspheres';
headmodel.o     = centers;
headmodel.r     = radii;
headmodel.unit  = mesh (1).unit;
headmodel.cond  = conductivities;
headmodel.label = sens.label;

% Restores the original warnings.
warning ( w )



function mesh = fixmesh ( mesh )

% If 'bnd' field takes only this field.
if isfield ( mesh, 'bnd' )
    
    % If no 'unit' field in the 'bnd' definitions uses the global one.
    if ~isfield ( mesh.bnd, 'unit' )
        
        % If mesh units takes those.
        if isfield ( mesh, 'unit' )
            [ mesh.bnd.unit ] = deal ( mesh.unit );
            
        % Otherwise tries to estimates them.
        else
            mesh.bnd = ft_convert_units ( mesh.bnd );
        end
    end
    
    % Takes the 'bnd' field.
    mesh = mesh.bnd;
end

% Replaces the 'pnt' field for 'pos'.
if isfield ( mesh, 'pnt' )
    for mindex = 1: numel ( mesh )
        mesh ( mindex ).pos = mesh ( mindex ).pnt;
    end
    mesh = rmfield ( mesh, 'pnt' );
end

function mesh = sortmesh ( mesh )




