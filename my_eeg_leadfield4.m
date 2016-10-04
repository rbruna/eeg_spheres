function [lf, vol] = my_eeg_leadfield4(R, elc, vol)

% EEG_LEADFIELD4 electric leadfield for a dipole in 4 concentric spheres
% 
% [lf] = eeg_leadfield4(R, elc, vol)
%
% with input arguments
%   R          position of the dipole
%   elc        position of the electrodes
% and vol being a structure with the elements
%   vol.r      radius of the 4 spheres 
%   vol.cond   conductivity of the 4 spheres
%   vol.t      constant factors for series expansion (optional)
%
% The center of the spheres should be at the origin.
%
% This implementation is adapted from
%   Lutkenhoner, Habilschrift 1992.
% The original reference is
%  Cuffin BN, Cohen D. Comparison of the magnetoencephalogram and electroencephalogram. Electroencephalogr Clin Neurophysiol. 1979 Aug;47(2):132-46. 
%
% See also EEG_LEADFIELD4_PREPARE for precomputing the constant factors,
% which can save time when multiple leadfield computations are done.

% Copyright (C) 2002, Robert Oostenveld
% Copyright (C) 2016, Ricardo Bruna
%
% This file is part of FieldTrip, see http://www.fieldtriptoolbox.org
% for the documentation and details.
%
%    FieldTrip is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    FieldTrip is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with FieldTrip. If not, see <http://www.gnu.org/licenses/>.
%
% $Id$


% Sorts the spheres from the smallest to the largest.
[ ~, indx ] = sort ( vol.r );
vol.cond    = vol.cond ( indx );
vol.r       = vol.r    ( indx );

% If the sphere is not centered in the origin corrects it.
if isfield ( vol, 'o' ) && sum ( vol.o .^ 2 ) > 1e-10
    R     = bsxfun ( @minus, R,     vol.o );
    elc   = bsxfun ( @minus, elc,   vol.o );
    vol.o = bsxfun ( @minus, vol.o, vol.o );
end


% Extracts the number of dipoles.
ndipoles = size ( R, 1 );

if ndipoles > 1
    error ( 'This function works with only one dipole.' );
end

% % Checks that the dipole is inside the inner sphere.
% if any ( sqrt ( sum ( R .^ 2, 2 ) ) > vol.r (1) )
%     error ( 'There are dipoles is outside the brain.' );
% end


% Extracts the number of sensors.
nchans = size ( elc, 1 );

% Projects the electrodes to the surface of the sphere.
dist = sqrt ( sum ( elc .^ 2, 2 ) );
elc  = bsxfun ( @times, elc, vol.r (4) ) ./ dist;


% Computes the constant factors for the sphere.
if ~isfield ( vol, 't' )
    vol.t = ft_eeg_leadfield4_prepare ( vol );
end

% Extracts the number of terms of the expansion.
nterms  = numel ( vol.t );



% use more convenient names for the radii and conductivities
r4 = vol.r(4); c4 = vol.cond(4);
n      = 1:nterms;
f      = norm(R)/r4;

% Generates the constant factors.
% The constant factors are immune to rotations.
const = (2*n+1).^4.*f.^(n-1) ./ (vol.t.*4*pi*c4*r4^2);






if R (1) ~= 0 || R (2) ~= 0
  % compute the rotation matrix
  % the inverse rotation matrix is the transposed of this one
  val1 = norm(R);
  val2 = norm(R(1:2));
  rot(1,1) = R(1) * R(3) / (val1 * val2); 
  rot(1,2) = R(2) * R(3) / (val1 * val2);
  rot(1,3) = -1.0 * val2 / val1;
  rot(2,1) = -1.0 * R(2) / val2;
  rot(2,2) =        R(1) / val2;
  rot(2,3) =                  0; 
  rot(3,:) = R ./ val1;
  % rotate the electrodes
  elc = elc*rot';
elseif R (3) < 0
  % dipole on negative z-axis, rotation is very simple: around x-axis
  elc(2,:) = -elc(2,:);
  elc(3,:) = -elc(3,:);
else
  % dipole is on positive z-axis, nothing has to be done
end


% Gets the electrode's position in spherical coordinates.
[ phi, the ] = cart2sph ( elc ( :, 1 ), elc ( :, 2 ), elc ( :, 3 ) );
cos_the = cos ( pi / 2 - the );

P0 = zeros ( nchans, nterms + 1 );
P1 = zeros ( nchans, nterms );
for cindex = 1: nchans
    P0 ( cindex, : ) = my_plgndr ( nterms, 0, cos_the ( cindex ), 1 ); % Zero'th order Legendre.
    P1 ( cindex, : ) = my_plgndr ( nterms, 1, cos_the ( cindex ), 1 ) ./ ( 1: nterms ); % First order Legendre.
end

s_x = bsxfun ( @times, const, P1 ); % s_y is identical.
s_z = bsxfun ( @times, const, P0 ( :, 2: end ) );

lf          = zeros ( nchans, 3 );
lf ( :, 1 ) = sum ( s_x, 2 ) .* -cos ( phi );
lf ( :, 2 ) = sum ( s_x, 2 ) .* -sin ( phi ); % s_y is identical to s_x
lf ( :, 3 ) = sum ( s_z, 2 );

% apply the inverse rotation to the leadfield matrix
if R(1)~=0 || R(2)~=0
  lf = lf*rot;
elseif R(1)==0 && R(2)==0 && R(3)<0
  lf(2,:) = -lf(2,:);
  lf(3,:) = -lf(3,:);
end

