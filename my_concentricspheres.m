function headmodel = my_concentricspheres ( mesh )

% Based on FieldTrip 20160222 functions:
% * ft_headmodel_concentricspheres

% FT_HEADMODEL_CONCENTRICSPHERES creates a volume conduction model
% of the head based on three or four concentric spheres. For a 3-sphere
% model the spheres represent the skin surface, the outside of the
% skull and the inside of the skull For a 4-sphere model, the surfaces
% describe the skin, the outside-skull, the inside-skull and the inside of the
% cerebro-spinal fluid (CSF) boundaries.
%
% The innermost surface is sometimes also referred to as the brain
% surface, i.e. as the outside of the brain volume.
%
% This function takes as input a single headshape described with
% points and fits the spheres to this surface. If you have a set of
% points describing each surface, then this function fits the spheres
% to all individual surfaces.
%
% Use as
%   headmodel = ft_headmodel_concentricspheres(mesh)

% Defines the default conductivities.
conductivities = [
    0.3300, ... % Brain conductivity.
    1.0000, ... % CSF conductivity.
    0.0042, ... % Skull conductivity.
    0.3300 ];   % Skin conductivity.


% Sanitices the meshes.
mesh = fixmesh ( mesh );

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


% Finds the center of all the spheres as the center of all points.
pos    = unique ( cat ( 1, mesh.pos ), 'rows');
center = my_fitsphere ( pos );
radii  = zeros ( size ( mesh ) );

% Fits each sphere radius as the the mean distance to the mesh's points.
for mindex = 1: numel ( mesh )
  dist     = sqrt ( sum ( bsxfun ( @minus, mesh ( mindex ).pos, center ) .^ 2, 2 ) );
  radii ( mindex ) = mean ( dist );
end

% Sorts the spheres from the smallest to the largest radius.
radii = sort ( radii );


% Sets the output.
headmodel      = [];
headmodel.type = 'concentricspheres';
headmodel.o    = center;
headmodel.r    = radii;
headmodel.unit = mesh (1).unit;
headmodel.cond = conductivities;



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

