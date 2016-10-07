function leadfield = my_leadfield_eegspheres ( dips, sens, headmodel )

% Leadfield for a EEG concentric spheres model.
% 
% Use as:
%   leadfield = my_leadfield_eegspheres ( dippos, senspos, headmodel )
%
% Where:
%   dippos      Position of the dipole.
%   senspos     Position of the electrodes.
%   headmodel   FieldTrip concentric spheres definition:
%       headmodel.o     Center of the spheres (optional if origin).
%       headmodel.r     Radius of the spheres.
%       headmodel.cond  Conductivity of each sphere.
%       headmodel.t     Series expansion for Gamma (optional).
%
% Electrodes must be on the surface of the outter sphere. Otherwise it
% position is radially projected over this surface.
%
% This function requires FieldTrip 20160222 or newer to work properly.

% This implementation is adapted from:
%   Cuffin and Cohen. Comparison of the Magnetoencephalogram and the Electroencephalogram. Electroencephalogr Clin Neurophysiol, 47:131-146.
% The method of images is adapted from:
%   Ermer et al 2001. Rapidly recomputable EEG forward models for realistic head shapes. Phys. Med. Biol. 46:1265–1281.
%   Cheng 1989. Field and Wave Electromagnetics 2nd edn. pp 141–143.

% Based on FieldTrip 20160222 functions:
% * eeg_leadfield4 by Robert Oostenveld

% Copyright (C) 2016, Ricardo Bruna

% Sanitizes the volume conductor.
headmodel   = fix_volume ( headmodel );

% Centers the sphere in the origin.
dips        = bsxfun ( @minus, dips,  headmodel.o );
sens        = bsxfun ( @minus, sens,  headmodel.o );
headmodel.o = bsxfun ( @minus, headmodel.o, headmodel.o );

% Gets the radius and the conductivity of the outter sphere.
oradius     = headmodel.r    (4);
ocond       = headmodel.cond (4);

% Gets the radius of the inner sphere.
iradius     = headmodel.r    (1);

% Gets the terms of the expansion.
terms       = headmodel.t;


% Extracts the number of terms of the expansion.
nterms      = numel ( terms );

% Extracts the number of dipoles and sensors.
ndips       = size ( dips, 1 );
nsens       = size ( sens, 1 );

% Calculates the radius of each dipole relative to the outter sphere.
% radii       = sqrt ( sum ( dips .^ 2, 2 ) ) / radius;
radii       = zeros ( ndips, 1, 'single' );
for c = 1: ndips
    radii (c)      = norm ( dips ( c, : ) );
end

% Calculates the correction factor for sources outside the inner sphere.
outcors = iradius ./ radii;
outcors = min ( outcors, 1 );

% if any ( outcors < 1 )
%     fprintf ( 1, 'Some of the dipoles are outside the inner sphere of the model. Fixing it using the method of images.\n' );
% end

% Uses the method of images to moves the dipoles inside the inner sphere.
radii =  radii .* outcors .^ 2;

% Generates the costant factors for each dipole.
% Cuffin & Cohen 1979. Eq A2. Part 1.
const1      = ( 2 * ( 1: nterms ) + 1 ) .^ 4;
const2      = bsxfun ( @power, radii / oradius, ( 1: nterms ) - 1 );
const3      = terms .* 4 * pi * ocond * oradius ^ 2;
const       = bsxfun ( @rdivide, bsxfun ( @times, const1, const2 ), const3 );


% Initializes the leadfield matrix.
leadfield = zeros ( nsens, 3, ndips, 'single' );

% Goes through each dipole.
for dindex = 1: ndips
    
    % Gets the current dipole.
    dippos      = dips ( dindex, : );

    % Rotates the system in order to put the dipole in the positive z-axis.
    % This separates the radial (z) and tangential (xy) components.
    
    % Initializes the rotation matrix.
    rot         = eye ( 3, 'single' );
    
    % If dipole not in the z-axis rotates it.
    if dippos (1) ~= 0 || dippos (2) ~= 0
        val1         = norm ( dippos );
        val2         = norm ( dippos ( 1: 2 ) );
        rot ( 1, 1 ) = dippos (1) * dippos (3) / ( val1 * val2 );
        rot ( 1, 2 ) = dippos (2) * dippos (3) / ( val1 * val2 );
        rot ( 1, 3 ) = -1.0 * val2 / val1;
        rot ( 2, 1 ) = -1.0 * dippos (2) / val2;
        rot ( 2, 2 ) =        dippos (1) / val2;
        rot ( 2, 3 ) =                  0;
        rot ( 3, : ) = dippos / val1;
        
    % If dipole in negative z rotates it around the x-axis.
    elseif dippos (3) < 0
        rot ( 2, 2 ) = -1;
        rot ( 3, 3 ) = -1;
    end
    
    % Creates a rotated version of the electrodes definition.
    dipsens     = sens * rot';
    
    
    % Gets the electrode's position in spherical coordinates.
    [ phi, theta ] = cart2sph ( dipsens ( :, 1 ), dipsens ( :, 2 ), dipsens ( :, 3 ) );
    cos_theta   = cos ( pi / 2 - theta );
    
    P0          = zeros ( nsens, nterms + 1 );
    P1          = zeros ( nsens, nterms );
    for cindex = 1: nsens
        
        % Zeroth order Legendre polynomial.
        P0 ( cindex, : ) = my_plgndr ( nterms, 0, cos_theta ( cindex ), 1 );
        
        % First order Legendre polynomial.
        P1 ( cindex, : ) = my_plgndr ( nterms, 1, cos_theta ( cindex ), 1 );
    end
    
    % Discards P0_0.
    P0 ( :, 1 ) = [];
    
    % Corrects P1.
    % Cuffin & Cohen 1979. Eq A2. Part 2.
    P1  = bsxfun ( @rdivide, P1, 1: nterms );
    
    % Creates the radial and tangential components.
    % Cuffin & Cohen 1979. Eq A2. Part 3.
    s_r         = bsxfun ( @times, const ( dindex, : ), P0 );
    s_t         = bsxfun ( @times, const ( dindex, : ), P1 );
    
    % Creates the leadfield for each direction of the dipole.
    % Cuffin & Cohen 1979. Eq A2. Part 4.
    lf          = zeros ( nsens, 3, 'single' );
    lf ( :, 1 ) = sum ( s_t, 2 ) .* -cos ( phi );
    lf ( :, 2 ) = sum ( s_t, 2 ) .* -sin ( phi );
    lf ( :, 3 ) = sum ( s_r, 2 );
    
    % Rotates back the leadfield and stores it.
    leadfield ( :, :, dindex ) = lf * rot;
end

% Corrects the leadfield amplitud of the dipoles outside the sphere.
leadfield = bsxfun ( @times, leadfield, reshape ( outcors, 1, 1, [] ) );

% % Silences the dipoles outside the inner sphere.%
% if any ( outcors < 1 )
%     fprintf ( 1, 'Some of the dipoles are outside the inner sphere of the model. Setting its leadfield to zero.\n' );
%     leadfield ( :, :, outcors < 1 ) = 0;
% end

% Reshapes the leadfield in matrix form.
leadfield = leadfield ( :, : );



function headmodel = fix_volume ( headmodel )

% Defines the center of the spheres, if needed.
if ~isfield ( headmodel, 'o' )
    headmodel.o    = zeros ( 1, 3 );
end

% Sorts the spheres from the smallest to the largest.
[ ~, indx ]    = sort ( headmodel.r );
headmodel.cond = headmodel.cond ( indx );
headmodel.r    = headmodel.r    ( indx );

% If only one spher creates the othe three.
if numel ( headmodel.r ) == 1
    headmodel.cond = headmodel.cond ( [ 1 1 1 1 ] );
    headmodel.r    = headmodel.r    ( [ 1 1 1 1 ] );
end

% If only three spheres creates the fourth one.
if numel ( headmodel.r ) == 3
    headmodel.cond = headmodel.cond ( [ 1 2 3 3 ] );
    headmodel.r    = headmodel.r    ( [ 1 2 3 3 ] );
end

% Computes the constant factors for the sphere, if needed.
if ~isfield ( headmodel, 't' )
    headmodel.t    = my_leadfield_eeggamma ( headmodel );
end
