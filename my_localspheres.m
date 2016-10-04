function headmodel = my_localspheres ( mesh, sens, varargin )

% Based on FieldTrip 20160222 functions:
% * ft_headmodel_localspheres

% FT_HEADMODEL_LOCALSPHERES constructs a MEG volume conduction model in
% with a local sphere fitted to the head or brain surface for each separate
% channel
%
% This implements
%   Huang MX, Mosher JC, Leahy RM. "A sensor-weighted overlapping-sphere
%   head model and exhaustive head model comparison for MEG." Phys Med
%   Biol. 1999 Feb;44(2):423-40
%
% Use as
%   headmodel = ft_headmodel_localspheres(mesh, grad, ...)
%
% Optional arguments should be specified in key-value pairs and can include
%   radius    = number, radius of sphere within which headshape points will
%               be included for the fitting algorithm
%   maxradius = number, if for a given sensor the fitted radius exceeds
%               this value, the radius and origin will be replaced with the
%               single sphere fit
%   baseline  = number
%   feedback  = boolean, true or false
%
% See also FT_PREPARE_HEADMODEL, FT_PREPARE_VOL_SENS, FT_COMPUTE_LEADFIELD


% Sanitices the meshes.
mesh = fixmesh ( mesh );

% Makes sure that mesh and sensors have the same geometrical units.
sens = ft_convert_units ( sens, mesh.unit );


% This function accepts only one mesh.
if numel ( mesh ) > 1
    error ( 'This function accepts only one mesh.' );
end


% Sets the default parameters in the desired untis.
threshold = ft_getopt(varargin, 'radius',    ft_scalingfactor('cm', mesh.unit) * 8.5);
maxradius = ft_getopt(varargin, 'maxradius', ft_scalingfactor('cm', mesh.unit) * 20);
baseline  = ft_getopt(varargin, 'baseline',  ft_scalingfactor('cm', mesh.unit) * 5);

% get the additional inputs and set the defaults
singlesphere  = ft_getopt(varargin, 'singlesphere', 'no');



% Initializes the center and radius arrays.
nchans  = numel(sens.label);
radii   = zeros ( nchans, 1 );
centers = zeros ( nchans, 3 );


% Fits a single sphere to all the points.
[ center, radius ] = my_fitsphere ( mesh.pos );

% If only one sphere is required that is all.
if strcmp ( singlesphere, 'yes' )
    
    % Sets the output.
    headmodel       = [];
    headmodel.type  = 'localspheres';
    headmodel.o     = center;
    headmodel.r     = radius;
    headmodel.unit  = mesh.unit;
end


for chan = 1: nchans
    
    % In MEG looks for the coil nearer to the head (?).
    if isfield ( sens, 'coilpos')
        
        % Gets all the coils of the current channel.
        coilsel = find ( sens.tra ( chan, : ) );
        allpos  = sens.coilpos ( coilsel, : );
        allori  = sens.coilori ( coilsel, : );
        
        % Gets the average position of the coils.
        thispos = mean ( allpos, 1 );
        
        % Gets the average orientation of the coils.
        [ ~, ~, v ] = svd ( allori );
        thisori = v ( :, 1 )';
        
        % Ensures that the orientation points outwards.
        thisori = thisori * sign ( dot ( thispos, thisori ) );
        
        % Computes the distance from every coil along that orientation.
        rawdist  = bsxfun ( @minus, allpos, thispos );
        distance = sum ( bsxfun ( @times, rawdist, thisori ), 2 );
        
        % Checs the distance from the point to all the coils.
        [ m, i ] = min ( distance );
        
        % If too far selects the nearest coil as center of the coils.
        if abs ( m ) > ( baseline / 4 )
            thispos = allpos ( i, : );
        end
        
    % In EEG gets the position of the electrode.
    else
        thispos = sens.chanpos ( chan, : );
    end
    
    
    % Find the head shape points closer to the current sensor.
    distance = sqrt ( sum ( bsxfun ( @minus, mesh.pos, thispos ) .^ 2, 2 ) );
    sphereps = mesh.pos ( distance < threshold, : );
    
    % Fits a sphere to those points.
    if size ( sphereps, 1 ) > 10
        [ o, r ] = my_fitsphere ( sphereps );
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
    
    % Stores the sphere center and radius.
    centers ( chan, : ) = o;
    radii   ( chan )    = r;
end


% Sets the output.
headmodel       = [];
headmodel.type  = 'localspheres';
headmodel.o     = centers;
headmodel.r     = radii;
headmodel.unit  = mesh.unit;
headmodel.label = sens.label;



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

