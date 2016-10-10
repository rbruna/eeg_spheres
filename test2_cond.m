clc
clear
close all

% Adds the functions folder to the path (in development).
addpath ( sprintf ( '%s/functions', pwd ) );

% Adds FT to the path.
ft_path;
ft_defaults

ft_hastoolbox ( 'freesurfer', 1, 1 );
ft_hastoolbox ( 'spm8', 1, 1 );

clear

% Loads the EEG data.
eegdata = load ( 'alc02_restingOA_EEG' );
eegdata.trialdata.elec = ft_convert_units ( eegdata.trialdata.elec, 'm' );
eegdata.trialdata.grad = ft_convert_units ( eegdata.trialdata.grad, 'm' );
eegdata.trialdata.elec.elecpos ( isnan ( eegdata.trialdata.elec.elecpos ) ) = 0;

% Loads the subject segmentation and the generated meshes.
mridata = load ( 'alc02_mri.mat' );
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


% Loads the EEG data.
eegdata = load ( 'alc02_restingOA_EEG' );
eegdata.trialdata.elec = ft_convert_units ( eegdata.trialdata.elec, 'm' );
eegdata.trialdata.grad = ft_convert_units ( eegdata.trialdata.grad, 'm' );
eegdata.trialdata.elec.elecpos ( isnan ( eegdata.trialdata.elec.elecpos ) ) = 0;

% Loads the subject segmentation and the generated meshes.
mridata = load ( 'alc02_mri.mat' );
mridata.mesh = ft_convert_units ( mridata.mesh, 'mm' );
mridata.mesh = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.mesh );
mridata.mesh = ft_convert_units ( mridata.mesh, 'm' );
mridata.grid = ft_convert_units ( mridata.grid, 'mm' );
mridata.grid = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.grid );
mridata.grid = ft_convert_units ( mridata.grid, 'm' );

% Fits the three concentric spheres to the meshes.
headmodel = ft_headmodel_concentricspheres ( mridata.mesh );

% Selects only the 60 first EEG channels.
sens = eegdata.trialdata.elec;
senpos = sens.elecpos;
senpos = senpos ( 1: 60, : );

% Selects the dipoles.
dippos = mridata.grid.pos;


tic
leadfield2 = my_leadfield_eegspheres ( dippos, senpos, headmodel );

% Applies the average reference.
leadfield2 = bsxfun ( @minus, leadfield2, mean ( leadfield2, 1 ) );
toc


tic
leadfield3 = zeros ( size ( senpos, 1 ), 3 * size ( dippos, 1 ), 'single' );
for c = 1: 60
    leadfield3 ( c, : ) = my_leadfield_eegspheres ( dippos, senpos ( c, : ), headmodel );
end

% Applies the average reference.
leadfield3 = bsxfun ( @minus, leadfield3, mean ( leadfield3, 1 ) );
toc
% return



% Loads the EEG data.
eegdata = load ( 'alc02_restingOA_EEG' );
eegdata.trialdata.elec = ft_convert_units ( eegdata.trialdata.elec, 'm' );
eegdata.trialdata.grad = ft_convert_units ( eegdata.trialdata.grad, 'm' );
eegdata.trialdata.elec.elecpos ( isnan ( eegdata.trialdata.elec.elecpos ) ) = 0;

% Loads the subject segmentation and the generated meshes.
mridata = load ( 'alc02_mri.mat' );
mridata.mesh = ft_convert_units ( mridata.mesh, 'mm' );
mridata.mesh = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.mesh );
mridata.mesh = ft_convert_units ( mridata.mesh, 'm' );
mridata.grid = ft_convert_units ( mridata.grid, 'mm' );
mridata.grid = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.grid );
mridata.grid = ft_convert_units ( mridata.grid, 'm' );

% Fits the three concentric spheres to the meshes.
headmodel = ft_headmodel_concentricspheres ( mridata.mesh );


cfg = [];
cfg.grid = mridata.grid;
cfg.senstype = 'eeg';
cfg.headmodel = headmodel;
cfg.channel = { 'all' '-EEG061' '-EEG062' '-EEG063' '-EEG064' };
cfg.elec = eegdata.trialdata.elec;
cfg.sens = eegdata.trialdata.elec;
cfg.feedback = 'no';

leadfield = ft_prepare_leadfield ( cfg );
leadfield = cat ( 2, leadfield.leadfield {:} );





% Loads the EEG data.
eegdata = load ( 'alc02_restingOA_EEG' );
eegdata.trialdata.elec = ft_convert_units ( eegdata.trialdata.elec, 'm' );
eegdata.trialdata.grad = ft_convert_units ( eegdata.trialdata.grad, 'm' );
eegdata.trialdata.elec.elecpos ( isnan ( eegdata.trialdata.elec.elecpos ) ) = 0;

% Loads the subject segmentation and the generated meshes.
mridata = load ( 'alc02_mri.mat' );
mridata.mesh = ft_convert_units ( mridata.mesh, 'mm' );
mridata.mesh = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.mesh );
mridata.mesh = ft_convert_units ( mridata.mesh, 'm' );
mridata.grid = ft_convert_units ( mridata.grid, 'mm' );
mridata.grid = ft_transform_geometry ( eegdata.mriinfo.transform, mridata.grid );
mridata.grid = ft_convert_units ( mridata.grid, 'm' );
% mridata.grid.pos ( 1: end-10, : ) = [];
% mridata.grid.inside ( 1: end-10 ) = [];


% Fits the three concentric spheres to the meshes.
% headmodel = ft_headmodel_concentricspheres ( mridata.mesh );
headmodel = my_headmodel_eegspheres ( mridata.mesh, eegdata.trialdata.elec, 'singlesphere', 'yes' );
headmodellc = my_headmodel_eegspheres ( mridata.mesh, eegdata.trialdata.elec );
headmodellc2 = my_headmodel_eegspheresB ( mridata.mesh, eegdata.trialdata.elec );
% idx = 1: 4;
% headmodellc.o ( idx, : ) = headmodel.o ( ones ( size ( idx ) ), : );
% headmodellc.r ( idx, : ) = headmodel.r ( ones ( size ( idx ) ), : );

cfg = [];
cfg.grid = mridata.grid;
cfg.senstype = 'eeg';
cfg.headmodel = headmodel;
cfg.channel = { 'all' '-EEG061' '-EEG062' '-EEG063' '-EEG064' };
cfg.elec = eegdata.trialdata.elec;
cfg.sens = eegdata.trialdata.elec;
cfg.feedback = 'no';

tic
leadfieldX = my_leadfield ( cfg );
leadfieldX = cat ( 2, leadfieldX.leadfield {:} );
leadfieldX = leadfieldX';
toc

cfg2 = cfg;
cfg2.headmodel = headmodellc;

tic
leadfieldY = my_leadfield ( cfg2 );
leadfieldY = cat ( 2, leadfieldY.leadfield {:} );
leadfieldY = leadfieldY';
toc

cfg3 = cfg;
cfg3.headmodel = headmodellc2;

tic
leadfieldZ = my_leadfield ( cfg3 );
leadfieldZ = cat ( 2, leadfieldZ.leadfield {:} );
leadfieldZ = leadfieldZ';
toc

leaddata = load('alc02_restingOA_EEG_leadfield.mat');
leadfieldOM = cat ( 3, leaddata.grid.leadfield {:} );
leadfieldOM = cat ( 1, zeros ( 1, 3, 2459 ), leadfieldOM );

lfX  = reshape ( leadfieldX', 60, 3, [] );
lfY  = reshape ( leadfieldY', 60, 3, [] );
lfZ  = reshape ( leadfieldZ', 60, 3, [] );
lfOM = bsxfun ( @minus, leadfieldOM, mean ( leadfieldOM, 1 ) );
lfX  = sqrt ( sum ( lfX .^ 2, 2 ) );
lfY  = sqrt ( sum ( lfY .^ 2, 2 ) );
lfZ  = sqrt ( sum ( lfZ .^ 2, 2 ) );
lfOM = sqrt ( sum ( lfOM .^ 2, 2 ) );
figure
plot ( [ lfZ(:) lfOM(:) ] )
legend ( { 'Multiple spheres' 'OM' } )
figure
plot ( [ lfX(:) lfZ(:) lfOM(:) ] )
legend ( { 'Sphere' 'Multiple spheres' 'OM' } )
