clc
clear
close all

% Adds the functions folder to the path (in development).
addpath ( sprintf ( '%s/functions', pwd ) );
addpath ( sprintf ( '%s/xtra', pwd ) );

% Adds FT to the path.
ft_path;
ft_defaults

ft_hastoolbox ( 'freesurfer', 1, 1 );
ft_hastoolbox ( 'spm8', 1, 1 );

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
headmodel = ft_headmodel_concentricspheres ( mridata.mesh );

% Selects only the 60 first EEG channels.
sens   = eegdata.trialdata.elec;
senpos = sens.elecpos;
senpos = senpos ( 1: 60, : );
senpos = bsxfun ( @minus, senpos, headmodel.o );

% Selects the dipoles.
dippos = mridata.grid.pos;
dippos = bsxfun ( @minus, dippos, headmodel.o );

headmodel.o = headmodel.o - headmodel.o;
headmodel.r (4) = headmodel.r (3);
headmodel.cond (4) = headmodel.cond (3);

% Places (projects) the sensors on the surface of the sphere.
dist = sqrt ( sum ( senpos .^ 2, 2 ) );
senposp = bsxfun ( @times, senpos, headmodel.r (4) ./ dist );

tic
leadfield0 = zeros ( size ( senpos, 1 ), 3 * size ( dippos, 1 ) );
for i = 1: size ( dippos, 1 )
    for c = 1: size ( senpos, 1 )
        leadfield0 ( c, ( i - 1 ) * 3 + ( 1: 3 ) ) = ft_eeg_leadfield4 ( dippos ( i, : ), senpos ( c, : ), headmodel );
    end
end

% Applies the average reference.
leadfield0 = bsxfun ( @minus, leadfield0, mean ( leadfield0, 1 ) );
toc

tic
leadfield1 = zeros ( size ( senpos, 1 ), 3 * size ( dippos, 1 ) );
for i = 1: size ( dippos, 1 )
    for c = 1: size ( senpos, 1 )
        leadfield1 ( c, ( i - 1 ) * 3 + ( 1: 3 ) ) = my_eeg_leadfield4 ( dippos ( i, : ), senpos ( c, : ), headmodel );
    end
end

% Applies the average reference.
leadfield1 = bsxfun ( @minus, leadfield1, mean ( leadfield1, 1 ) );
toc

tic
leadfield2 = mymsc_leadfield ( dippos, senpos, headmodel, false );

% Applies the average reference.
leadfield2 = bsxfun ( @minus, leadfield2, mean ( leadfield2, 1 ) );
toc

figure
hold on
plot ( leadfield0 (:) );
plot ( leadfield1 (:) );
plot ( leadfield2 (:) );
