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
headmodel    = ft_headmodel_concentricspheres ( mridata.mesh );
headmodelsc  = mymcs_headmodel ( mridata.mesh, eegdata.trialdata.elec, 'singlesphere', 'yes' );
headmodellc  = mymcs_headmodel ( mridata.mesh, eegdata.trialdata.elec );
headmodellc2 = mymcs_headmodel ( mridata.mesh, eegdata.trialdata.elec, 'eqradii', false );


cfg = [];
cfg.grid = mridata.grid;
cfg.senstype = 'eeg';
cfg.headmodel = headmodel;
cfg.channel = { 'all' '-EEG061' '-EEG062' '-EEG063' '-EEG064' };
cfg.elec = eegdata.trialdata.elec;
cfg.sens = eegdata.trialdata.elec;
cfg.feedback = 'no';


cfgA0 = cfg;
cfgA0.headmodel = headmodel;

tic
leadfieldA0 = ft_prepare_leadfield ( cfg );
leadfieldA0 = cat ( 2, leadfieldA0.leadfield {:} );
leadfieldA0 = leadfieldA0';
toc


cfgB0 = cfg;
cfgB0.headmodel = headmodel;
cfgB0.moveinside = false;

tic
leadfieldB0 = my_leadfield ( cfgB0 );
leadfieldB0 = cat ( 2, leadfieldB0.leadfield {:} );
leadfieldB0 = leadfieldB0';
toc

cfgB1 = cfg;
cfgB1.headmodel = headmodel;
cfgB1.moveinside = true;

tic
leadfieldB1 = my_leadfield ( cfgB1 );
leadfieldB1 = cat ( 2, leadfieldB1.leadfield {:} );
leadfieldB1 = leadfieldB1';
toc


cfgC0 = cfg;
cfgC0.headmodel = headmodelsc;
cfgC0.moveinside = false;

tic
leadfieldC0 = my_leadfield ( cfgC0 );
leadfieldC0 = cat ( 2, leadfieldC0.leadfield {:} );
leadfieldC0 = leadfieldC0';
toc

cfgC1 = cfg;
cfgC1.headmodel = headmodelsc;
cfgC1.moveinside = true;

tic
leadfieldC1 = my_leadfield ( cfgC1 );
leadfieldC1 = cat ( 2, leadfieldC1.leadfield {:} );
leadfieldC1 = leadfieldC1';
toc


cfgD0 = cfg;
cfgD0.headmodel = headmodellc;
cfgD0.moveinside = false;

tic
leadfieldD0 = my_leadfield ( cfgD0 );
leadfieldD0 = cat ( 2, leadfieldD0.leadfield {:} );
leadfieldD0 = leadfieldD0';
toc

cfgD1 = cfg;
cfgD1.headmodel = headmodellc;
cfgD1.moveinside = true;

tic
leadfieldD1 = my_leadfield ( cfgD1 );
leadfieldD1 = cat ( 2, leadfieldD1.leadfield {:} );
leadfieldD1 = leadfieldD1';
toc


cfgE0 = cfg;
cfgE0.headmodel = headmodellc2;
cfgE0.moveinside = false;

tic
leadfieldE0 = my_leadfield ( cfgE0 );
leadfieldE0 = cat ( 2, leadfieldE0.leadfield {:} );
leadfieldE0 = leadfieldE0';
toc

cfgE1 = cfg;
cfgE1.headmodel = headmodellc2;
cfgE1.moveinside = true;

tic
leadfieldE1 = my_leadfield ( cfgE1 );
leadfieldE1 = cat ( 2, leadfieldE1.leadfield {:} );
leadfieldE1 = leadfieldE1';
toc


leadfieldOM = load ( '../../data/alc02_restingOA_EEG_leadfield.mat' );
leadfieldOM = cat ( 2, leadfieldOM.grid.leadfield {:} );
leadfieldOM = cat ( 1, zeros ( 1, 3 * 2459 ), leadfieldOM )';
leadfieldOM = bsxfun ( @minus, leadfieldOM, mean ( leadfieldOM, 1 ) );


lfA0 = reshape ( leadfieldA0', 60, 3, [] );
lfB0 = reshape ( leadfieldB0', 60, 3, [] );
lfB1 = reshape ( leadfieldB1', 60, 3, [] );
lfC0 = reshape ( leadfieldC0', 60, 3, [] );
lfC1 = reshape ( leadfieldC1', 60, 3, [] );
lfD0 = reshape ( leadfieldD0', 60, 3, [] );
lfD1 = reshape ( leadfieldD1', 60, 3, [] );
lfE0 = reshape ( leadfieldE0', 60, 3, [] );
lfE1 = reshape ( leadfieldE1', 60, 3, [] );
lfOM = reshape ( leadfieldE1', 60, 3, [] );

lfA0 = sqrt ( sum ( lfA0 .^ 2, 2 ) );
lfB0 = sqrt ( sum ( lfB0 .^ 2, 2 ) );
lfB1 = sqrt ( sum ( lfB1 .^ 2, 2 ) );
lfC0 = sqrt ( sum ( lfC0 .^ 2, 2 ) );
lfC1 = sqrt ( sum ( lfC1 .^ 2, 2 ) );
lfD0 = sqrt ( sum ( lfD0 .^ 2, 2 ) );
lfD1 = sqrt ( sum ( lfD1 .^ 2, 2 ) );
lfE0 = sqrt ( sum ( lfE0 .^ 2, 2 ) );
lfE1 = sqrt ( sum ( lfE1 .^ 2, 2 ) );
lfOM = sqrt ( sum ( lfOM .^ 2, 2 ) );

figure
hold on
plot ( lfA0 (:) )
plot ( lfB0 (:) )
plot ( lfB1 (:) )
legend ( { 'Concentric spheres FT' 'Concentric spheres' 'Concentric spheres (ourside corrected)' } )

figure
hold on
plot ( lfC0 (:) )
plot ( lfD0 (:) )
plot ( lfOM (:) )
legend ( { 'Concentric spheres' 'Multiple spheres equal' 'OM' } )

figure
hold on
plot ( lfC1 (:) )
plot ( lfD1 (:) )
plot ( lfE1 (:) )
plot ( lfOM (:) )
legend ( { 'Concentric spheres (ourside corrected)' 'Multiple spheres equal (ourside corrected)' 'Multiple spheres no equal (ourside corrected)' 'OM' } )
