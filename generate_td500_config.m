function lsnconfig = generate_td500_config
% generate_td500_config.m
% Generates Config file for TEMPORAL AND POPULATION DYNAMICS analysis of rsvp500 data
% written by AHB, May 2014
% updated Feb 16, 2024

global lsnconfig %#ok<*REDEFGI>

% ROOT DIRECTORIES
lsnconfig.ephysdata_rootdir = '/Volumes/WD_BLACK8TB mac/_CompleteArchive_June2023_ORIGof3/__Current_Projects/ephysData';
lsnconfig.in_rootdir = '/Volumes/WD_BLACK8TB mac/_CompleteArchive_June2023_ORIGof3/_LONG_TERM_STORAGE/_@_PLEXONDATA_NIH';
lsnconfig.out_rootdir = '~/OneDrive-King''sCollegeLondon/ephysProjects/T500_TemporalDynamics_Study';

% RAW DATA DIRECTORIES
lsnconfig.spikedir  = [in_rootdir,'COMPILED_output',filesep]; % output dir for spike matrices
lsnconfig.neurondir = [in_rootdir,'UNIT_output',filesep]; % output dir for MUA
lsnconfig.LFPdir    = [in_rootdir,'LFPS',filesep]; % output dir for uncorrected LFP channels
lsnconfig.LFPdir_corr=[in_rootdir,'LFPscorr',filesep]; % output dir for LFP channels (after reverse filtering)
lsnconfig.startend  = [in_rootdir,'StartEnd',filesep]; % output dir for start & end text files
lsnconfig.locationdir=[in_rootdir,'LocationFiles',filesep]; % output dir for location text files
lsnconfig.wavetemp   =[in_rootdir,'WAVEFORM_templates',filesep]; % output dir for waveform templates
lsnconfig.wave_raw   =[in_rootdir,'WAVEFORM_raw',filesep]; % output dir for raw waveforms
lsnconfig.unitdir   = [in_rootdir,'UnitFiles',filesep]; % output dir for spiketrains/trial files
lsnconfig.rsvp500spks=[in_rootdir,'rsvp500spks',filesep]; % output for spiketrains for dms400 neurons
lsnconfig.rsvp500lfps=[in_rootdir,'rsvp500lfps',filesep]; % output for lfps for rsvp500 files

% ANALYSIS/OUTPUT DIRECTORIES
lsnconfig.figurepath    = [out_rootdir,'figure_source_images',filesep]; % output dir for all data structures
lsnconfig.datadir = [out_rootdir,'matlab_data',filesep]; % output dir for all figures



% DEFAULT FILE LISTS
lsnconfig.excelfile = [out_rootdir,'td500_Neurons.xlsx']; % excel spreadsheet

% ANALYSIS DEFAULTS
lsnconfig.gausskernel = 10; % gaussian spike density function kernel (10ms)
lsnconfig.gausskernelsml = 5; % smaller gaussian spike density function kernel (10ms)
lsnconfig.printer     = 0;
lsnconfig.LFP_kernel  = [in_rootdir,'PRA2kernelEmpPRA2HST20Gelec_2kHz.mat'];
lsnconfig.LFP_kernel_HST1X = [in_rootdir,'PRA2kernel_2kHz.mat'];
lsnconfig.xscale      = -100:500; % default time window % this is for RASTER CREATION (Note: can't be any lower than 100 given short ITI)
lsnconfig.xScaleRange = [-100 500]; % this is for figures
lsnconfig.fontsize_sml=10; lsnconfig.fontsize_med=12; lsnconfig.fontsize_lrg=14;

% LFP ANALYSIS 
lsnconfig.chronux_params=struct('tapers',[3 5],'Fs',1000,'pad',-1,'err',0,'trialave',1,'fpass',[0 120]);

% PATTERN CLASSIFIER ANALYSIS
lsnconfig.pClass_time_skip = 10; % slides TIME_SKIP ms each step
lsnconfig.pClass_time_win = 20; % plus/minus TIME_WIN ms
lsnconfig.pClass_time_range = -100:lsnconfig.pClass_time_skip:500;

% BEHAV_MATRIX
lsnconfig.COL_parnum    = 3; % column for paradigm code
lsnconfig.COL_condno    = 4; % column for condition number
lsnconfig.COL_outcome   = 30; % column for trial outcome
lsnconfig.COL_CTOA      = 31; % column for CTOA
lsnconfig.COL_SRT       = 32; % column for SRT
lsnconfig.COL_remark    = 39; % column for trial remark
lsnconfig.COL_duration  = 33; % column for trial duration
lsnconfig.COL_cueduration = 34; % column for cue duration
lsnconfig.COL_starttrial= 40; % column for trial start time
lsnconfig.COL_endtrial  = 41; % column for trial end time
lsnconfig.COL_FPon      = 42; % column for FP onset time
lsnconfig.COL_FPoff     = 43; % column for FP offset time
lsnconfig.COL_CUEon     = 44; % column for cue onset time
lsnconfig.COL_CUEoff    = 45; % column for cue offset time
lsnconfig.COL_TARG1on   = 46; % column for target 1 onset time
lsnconfig.COL_TARG2on   = 48; % column for target 2 onset time
lsnconfig.COL_eyewind   = 50; % column for eye left fixation window time
lsnconfig.COL_reward    = 51; % column for reward given time
lsnconfig.COL_blocknum  = 53; % column for cortex assigned block number
lsnconfig.COL_trialnum  = 54; % column for cortex assigned trialnumber

% EVENT CODES - Lists the definitions for the event codes
% GENERAL CODES
lsnconfig.start_trial     = 9; % start of the trial
lsnconfig.eye_buffer_on   = 10; % eye buffer on
lsnconfig.end_iti         = 11; % end of the ITI - duration of ITI = end_iti - start_trial
lsnconfig.end_trial       = 101; % end of trial
lsnconfig.eye_left_fp     = 32; % eye left fp window

% STIMULUS CODES
lsnconfig.FP_on           = 23; % FP on
lsnconfig.FP_off          = 24; % FP off
lsnconfig.cue_on          = 25; % cue on
lsnconfig.cue_off         = 26; % cue off
lsnconfig.choice_on       = 27; % choice stimuli presented
lsnconfig.go_signal       = 30; % go signal (choices on, FP off)

% CORRECT CODES
lsnconfig.eye_on_fp       = 12; % eye on fixation point
lsnconfig.eye_on_targ     = 28; % eye on target
lsnconfig.reward          = 34; % reward given, correct trial (should be 3!!!)

% ERROR CODES
lsnconfig.never_on_fp     = 14; % broke fixation during first fixation period
lsnconfig.broke_1stfix    = 15; % broke fixation during first fixation period
lsnconfig.broke_cue       = 17; % broke fixation during cue presentation
lsnconfig.broke_2ndfix    = 18; % broke fixation during second fixation period
lsnconfig.eye_on_incorrect= 33; % eye on incorrect target
lsnconfig.eye_on_easyincorrect = 28; % selected easy but incorrect target
lsnconfig.no_movement     = 31; % monkey never left FP after choices appeared


% EVENT CODES - Lists the definitions for the event codes
% rsvp500
lsnconfig.faces500  = 1:20; % condition numbers corresponding to face stimuli
lsnconfig.fruit500  = 21:40;
lsnconfig.places500 = 41:60;
lsnconfig.bodyp500  = 61:80;
lsnconfig.objct500  = 81:100;
lsnconfig.mkfac500  = 1:20; % condition numbers corresponding to face stimuli
lsnconfig.hmfac500  = 21:40;
lsnconfig.places500 = 41:60;
lsnconfig.bodyp500  = 61:80;
lsnconfig.objct500  = 81:100;


% SAVE CONFIG FILE
configname = [out_rootdir,'td500_configplex.mat'];
save(configname,'lsnconfig');
