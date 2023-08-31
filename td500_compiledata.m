function [UnitData,monkeyname,grids]=td500_compiledata(monkeyname,sheetname)
% by AHB, May, 2014
global lsnconfig
%%%  LOAD FILE LIST
tic;
fprintf(['\nLoading file list from Excel: ',sheetname,'\n']);
[~,UnitData.plxname]=xlsread(lsnconfig.excelfile,sheetname,'B5:B1000'); % Filename
[~,UnitData.unitname]=xlsread(lsnconfig.excelfile,sheetname,'C5:C1000'); % Cell Number
[~,UnitData.gridloc]=xlsread(lsnconfig.excelfile,sheetname,'E5:E1000'); % Grid Location
UnitData.depth=xlsread(lsnconfig.excelfile,sheetname,'F5:F1000'); % Depth (um)
[~,UnitData.apindex]=xlsread(lsnconfig.excelfile,sheetname,'G5:G1000'); % AP Index
[~,UnitData.est_location]=xlsread(lsnconfig.excelfile,sheetname,'H5:H1000'); % Estimated Location
[~,UnitData.autoSensory]=xlsread(lsnconfig.excelfile,sheetname,'I5:I1000'); % Automated Sensory Classification
[~,UnitData.confSensory]=xlsread(lsnconfig.excelfile,sheetname,'J5:J1000'); % Confirmed Sensory Classification
[~,UnitData.autoPrefCat]=xlsread(lsnconfig.excelfile,sheetname,'K5:K1000'); % Automated Category Preference
[~,UnitData.confPrefCat]=xlsread(lsnconfig.excelfile,sheetname,'L5:L1000'); % Confirmed Category Preference
[~,UnitData.autoCatSelect]=xlsread(lsnconfig.excelfile,sheetname,'M5:M1000'); % Automated Category Selectivity
[~,UnitData.confCatSelect]=xlsread(lsnconfig.excelfile,sheetname,'N5:N1000'); % Confirmed Category Selectivity
[~,UnitData.autoRespDir]=xlsread(lsnconfig.excelfile,sheetname,'O5:O1000'); % Automated Response Direction (Excite/Inhibit)
[~,UnitData.confRespDir]=xlsread(lsnconfig.excelfile,sheetname,'P5:P1000'); % Confirmed Response Direction (Excite/Inhibit)
UnitData.quality=xlsread(lsnconfig.excelfile,sheetname,'Q5:Q1000'); % Quality
UnitData.inclWaveform=xlsread(lsnconfig.excelfile,sheetname,'AA1:AA1000'); % include in Waveform analysis?
toc

% Convert Grid Location to AP coordinates
for uu=1:length(UnitData.gridloc),
    UnitData.APcoords(uu,:)=plx_convertgrid2ap(UnitData.gridloc(uu));
end

%%% Load information from RESPSTRUCTSINGLE files
h = waitbar(0,'Loading data from individual neuron structures...'); tic
for un=1:size(UnitData.plxname,1),
    
    % Load individual file
    NewUnitName=char(UnitData.plxname(un)); newunit=char(UnitData.unitname(un));
    % disp(['...',NewUnitName(1:12),'-',newunit,'...'])
    load([lsnconfig.rsvp500spks,NewUnitName(1:12),'-',newunit,'-500_NeuronData.mat']); % load unit data
    waitbar(un/size(UnitData.plxname,1),h,['Loading data from ',NewUnitName(1:12),'-',newunit])
    
    % Labels and General Descriptives
    UnitData.UnitName{un}               = NewUnitName;
    UnitData.monkeyname{un}             = monkeyname;
    UnitData.label{un}                  = respstructsingle.label;
    UnitData.datemodified{un}           = respstructsingle.datemodified;
    UnitData.quality(un)                = respstructsingle.quality;
    UnitData.APcoords(un,:)             = respstructsingle.APcoords;
    UnitData.APIndex(un)                = respstructsingle.APIndex;
    % Calculate AP in number
    temp=char(UnitData.APIndex(un));
    UnitData.APindexNumber(un)          = str2num(temp(2:end));
    UnitData.depth(un)                  = respstructsingle.depth;
    UnitData.gridloc(un)                = respstructsingle.gridlocation;
    UnitData.confCatSelect(un)          = respstructsingle.conf_selective;
    UnitData.confPrefCat(un)            = respstructsingle.conf_preferred_cat;
    UnitData.confRespDir(un)            = respstructsingle.conf_excite;
    UnitData.confSensory(un)            = respstructsingle.conf_neurtype; % Sensory Classification (Visual vs. Non-responsive)
    UnitData.selective_nofruit{un}      = respstructsingle.selective_nofruit;
    UnitData.preferred_sensory{un}      = respstructsingle.preferred_sensory;
    UnitData.excite_inhibit(un,:)       = respstructsingle.excite_inhibit;
    UnitData.excitetype{un}             = respstructsingle.excitetype;
    UnitData.excitetype_nofruit{un}     = respstructsingle.excitetype_nofruit;
    % UnitData.prefcat_nofruit{un}        = respstructsingle.prefcat_nofruit;
    UnitData.preferred_category{un}     = respstructsingle.preferred_category;
    UnitData.pref_excite{un}            = respstructsingle.pref_excite;
    UnitData.prefcat_excite_nofruit{un} = respstructsingle.prefcat_excite_nofruit;
    UnitData.pref_inhibit{un}           = respstructsingle.pref_inhibit;
    UnitData.prefcat_inhibit_nofruit{un}= respstructsingle.prefcat_inhibit_nofruit;
    UnitData.face_trad(un)              = respstructsingle.face_trad;
    UnitData.face_trad_nofruit(un)      = respstructsingle.face_trad_nofruit;
    UnitData.validrsp(un,:)             = respstructsingle.validrsp(:,2); % is average condition response > 2*mean baseline?
    UnitData.wf_autoinclude(un)         = respstructsingle.wf_autoinclude;
    UnitData.wf_include(un)             = respstructsingle.wf_include;
    UnitData.wf_params(un,:)            = respstructsingle.wf_params;
    UnitData.wf_type(un)                = respstructsingle.wf_type;
    
    % Spiking Activity
    UnitData.m_baseline(un,:)           = respstructsingle.m_baseline;
    UnitData.m_epoch1(un,:)             = respstructsingle.m_epoch1; % 50-300ms
    UnitData.m_epoch1_nobase(un,:)      = respstructsingle.m_epoch1_nobase;
    UnitData.m_epoch2(un,:)             = respstructsingle.m_epoch2; % 100-400ms
    UnitData.m_epoch2_nobase(un,:)      = respstructsingle.m_epoch2_nobase;
    UnitData.m_epoch3(un,:)             = respstructsingle.m_epoch3; % 50-200ms
    UnitData.m_epoch3_nobase(un,:)      = respstructsingle.m_epoch3_nobase;
    
    UnitData.latency(un,:)              = respstructsingle.latency'; % face fruit place body object
    UnitData.cat_latency(un,:)          = reshape(respstructsingle.cat_latency,10,1); % face fruit place body object
    
    % Selectivity Indices
    UnitData.cat_si(un,:)               = respstructsingle.cat_si(:,2)'; % face fruit place body object (based on mean firing rate of epoch 1)
    UnitData.cat_si_nobase(un,:)        = reshape(respstructsingle.cat_si_nobase,12,1);
    UnitData.excite_rawsi(un)           = respstructsingle.excite_rawsi';
    UnitData.excite_rawsi_nofruit(un)   = respstructsingle.excite_rawsi_nofruit';
    UnitData.inhibit_rawsi(un)          = respstructsingle.inhibit_rawsi';
    UnitData.inhibit_rawsi_nofruit(un)  = respstructsingle.inhibit_rawsi_nofruit';
    UnitData.pure_si(un,:)                = reshape(respstructsingle.pure_si,25,1); % all five vs. all five
    
    % Statistics
    UnitData.anova_baseline(un,:)       = respstructsingle.anova_baseline;
    UnitData.anova_epoch(un,:)          = respstructsingle.anova_epoch;
    UnitData.anova_epoch_nobase(un,:)   = respstructsingle.anova_epoch_nobase;
    UnitData.anova_latency(un)          = respstructsingle.anova_latency';
    UnitData.anova_within_group(un,:)   = respstructsingle.anova_within_group(:,:,2)';
    UnitData.catanova_nofruit(un)       = respstructsingle.catanova_nofruit';
    UnitData.stats_prefexcite_v_others_nofruit(un)  = respstructsingle.stats_prefexcite_v_others_nofruit';
    UnitData.stats_prefinhibit_v_others_nofruit(un) = respstructsingle.stats_prefinhibit_v_others_nofruit';
    UnitData.stats_rsp_matrix_avg(un,:)             = reshape(respstructsingle.stats_rsp_matrix_avg,25,1);
    UnitData.stats_rsp_matrix_trial(un,:)           = reshape(respstructsingle.stats_rsp_matrix_trial,25,1);
    temproc = reshape(respstructsingle.roc_analysis,1,25);
    for rr=1:25, if temproc(rr)<0.5, temproc(rr)=1-temproc(rr); end; end % correct ROC values
    UnitData.roc_analysis(un,:)         = temproc;
    
    UnitData.cat_avg(un,:)              = respstructsingle.cat_avg(:,2)';
    UnitData.cat_sem(un,:)              = respstructsingle.cat_sem(:,2)';
    UnitData.cat_avg_nobase(un,:)       = respstructsingle.cat_avg_nobase(:,2)';
    UnitData.cat_sem_nobase(un,:)       = respstructsingle.cat_sem_nobase(:,2)';
    UnitData.norm_cat_avg(un,:)         = respstructsingle.norm_cat_avg';
    UnitData.cat_bst(un,:)              = reshape(respstructsingle.cat_bst,10,1);
    UnitData.cat_bst_nobase(un,:)       = reshape(respstructsingle.cat_bst_nobase,10,1);
    UnitData.cat_id(un,:)               = respstructsingle.cat_id;
    UnitData.cat_sensory(un,:)          = reshape(respstructsingle.cat_sensory,10,1);
    
    % LFP
    unitdata.LFP_trial_epoch=respstructsingle.LFP_trial_epoch;
    unitdata.LFP_trial_epoch_rect=respstructsingle.LFP_trial_epoch_rect;
    try % LFP fields (need catch loop in case LFP processing failed
        UnitData.LFP_lfp_average_epoch_rect(un,:)   = respstructsingle.LFP_lfp_average_epoch_rect;
        UnitData.LFP_cat_avg_epoch_tr(un,:)         = respstructsingle.LFP_cat_avg_epoch_tr;
        UnitData.LFP_cat_avg_rect_epoch(un,:)       = respstructsingle.LFP_cat_avg_rect_epoch;
        unitdata.LFP_cat_avg_rect_epoch_tr          = respstructsingle.LFP_cat_avg_rect_epoch_tr;
        unitdata.LFP_anova_stim                     = respstructsingle.LFP_anova_stim;
        UnitData.LFP_bestlabel(un)                  = respstructsingle.LFP_bestlabel;
        UnitData.LFP_cat_anova_rect(un,:)           = respstructsingle.LFP_cat_anova_rect;
        UnitData.LFP_evoked_cat_si(un,:)            = respstructsingle.LFP_evoked_cat_si;
        unitdata.LFP_evokedpure_cat_si(un,:)        = respstructsingle.LFP_evokedpure_cat_si;
        UnitData.LFP_freq_epoch_cat(un,:)           = reshape(respstructsingle.LFP_freq_epoch_cat,1,24);
        unitdata.LFP_freq_within_anova              = respstructsingle.LFP_freq_within_anova;
        UnitData.LFP_freq_bestlabel(un,:)           = respstructsingle.LFP_freq_bestlabel;
        UnitData.LFP_freq_across_anova(un,:)        = respstructsingle.LFP_freq_across_anova;
        UnitData.LFP_freq_cat_si(un,:)              = reshape(respstructsingle.LFP_freq_cat_si,1,24);
        unitdata.LFP_tr_min_max                     = respstructsingle.LFP_tr_min_max;
        unitdata.LFP_cond_min_max                   = respstructsingle.LFP_cond_min_max;
        unitdata.LFP_cat_min_max                    = respstructsingle.LFP_cat_min_max;
        UnitData.LFP_stats_pref_v_others_evoked(un) =respstructsingle.LFP_stats_pref_v_others_evoked';
        UnitData.LFP_stats_pref_v_others_0_120Hz(un)=respstructsingle.LFP_stats_pref_v_others_0_120Hz';
        UnitData.LFP_stats_pref_v_others_0_20Hz(un) =respstructsingle.LFP_stats_pref_v_others_0_20Hz';
    catch
        UnitData.LFP_lfp_average_epoch_rect(un,:)   = respstructsingle.LFP_lfp_average_epoch_rect;
        UnitData.LFP_cat_avg_epoch_tr(un,:)         = respstructsingle.LFP_cat_avg_epoch_tr;
        UnitData.LFP_cat_avg_rect_epoch(un,:)       = respstructsingle.LFP_cat_avg_rect_epoch;
        unitdata.LFP_cat_avg_rect_epoch_tr          = respstructsingle.LFP_cat_avg_rect_epoch_tr;
        unitdata.LFP_anova_stim                     = respstructsingle.LFP_anova_stim;
        UnitData.LFP_bestlabel(un)                  = respstructsingle.LFP_bestlabel;
        UnitData.LFP_cat_anova_rect(un,:)           = respstructsingle.LFP_cat_anova_rect;
        UnitData.LFP_evoked_cat_si(un,:)            = respstructsingle.LFP_evoked_cat_si;
        unitdata.LFP_evokedpure_cat_si(un,:)        = respstructsingle.LFP_evokedpure_cat_si;
        UnitData.LFP_freq_epoch_cat(un,:)           = reshape(respstructsingle.LFP_freq_epoch_cat,1,24);
        unitdata.LFP_freq_within_anova              = respstructsingle.LFP_freq_within_anova;
        UnitData.LFP_freq_bestlabel(un,:)           = respstructsingle.LFP_freq_bestlabel;
        UnitData.LFP_freq_across_anova(un,:)        = respstructsingle.LFP_freq_across_anova;
        UnitData.LFP_freq_cat_si(un,:)              = reshape(respstructsingle.LFP_freq_cat_si,1,24);
        unitdata.LFP_tr_min_max                     = respstructsingle.LFP_tr_min_max;
        unitdata.LFP_cond_min_max                   = respstructsingle.LFP_cond_min_max;
        unitdata.LFP_cat_min_max                    = respstructsingle.LFP_cat_min_max;
        UnitData.LFP_stats_pref_v_others_evoked(un) = respstructsingle.LFP_stats_pref_v_others_evoked';
        UnitData.LFP_stats_pref_v_others_0_120Hz(un)= respstructsingle.LFP_stats_pref_v_others_0_120Hz';
        UnitData.LFP_stats_pref_v_others_0_20Hz(un) = respstructsingle.LFP_stats_pref_v_others_0_20Hz';
    end
    
    % FMRI
    UnitData.AFNIcoords(un,:)=respstructsingle.AFNIcoords;
    %UnitData.AFNItimeseries(un,:)=respstructsingle.AFNItimeseries;
    UnitData.FMRIdat_coords(un,:)=respstructsingle.FMRIdat_coords;
    UnitData.fmri_rsp(un,:)=respstructsingle.fmri_rsp;
    UnitData.fmri_catsi(un,:)=respstructsingle.fmri_catsi;
    UnitData.fmri_excite_rawsi(un)=respstructsingle.fmri_excite_rawsi';
    UnitData.fmri_inhibit_rawsi(un)=respstructsingle.fmri_inhibit_rawsi';
    UnitData.fmri_norm_rsp(un,:) = respstructsingle.fmri_norm_rsp;
    
    % INSERT NEW FIELDS AS NEEDED...
    
    if strcmp(monkeyname,'Stewie')==1.
        
        % Grid location groups for comparison
        grids.grp(1).grids={'A7L2','A7L1','A6L3'}; % BodyPart Selective
        grids.grp(2).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Face Selective
        grids.grp(3).grids={'A4L2','A4L1','A4R1'}; % No Category Selectivity
        grids.grp(4).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4', 'P4L2','P4L4','P5L3','P6L3'}; % Object Selective
        grids.grp(5).grids={'P6L2','P6L1','P7L2'}; % Face Selective
        
        grids.sampledlocations={'A7R1','A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'};
        % Grid location groups for comparison
        grids.grpf(1).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Anterior Face Selective
        grids.grpf(2).grids={'P6L2','P6L1','P7L2'}; % Posterior Face Selective
        grids.grpf(3).grids={'A6L2','A6L0','A5L2','A5L1','A5L0','P6L2','P6L1','P7L2'}; % Inside All
        grids.grpf(4).grids={'A7L2','A7L1','A7R1','A4L2','A4L1','A4R1','A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P5L3','P6L3'}; % Outside All Face Selective
        grids.grpbp(1).grids={'A7L2','A7L1'}; % Inside Bodypart Selective
        grids.grpbp(2).grids={'A7R1','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Outside All Bodypart Selective
        grids.grpob(1).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L2'}; % Inside Object Selective
        grids.grpob(2).grids={'A7R1','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','P5L0','P6L1','P6L3','P7L2'}; % Outside All Object Selective
        grids.grppl(1).grids={'A4L2','A4L1','A4R1','A7R1'}; % Non Category Selective
        grids.grppl(2).grids={'A7L2','A7L1','A6L2','A6L0','A5L2','A5L1','A5L0','A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L3','P6L2','P6L1','P7L2'}; % Category Selective
        
        grids.grpfnf(1).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Anterior Face Selective
        grids.grpfnf(2).grids={'A7L2','A7L1','A7R1','A4L2','A4L1','A4R1'}; % Near Face Selective (Anterior)
        grids.grpfnf(3).grids={'P6L2','P6L1','P7L2'}; % Posterior Face Selective
        grids.grpfnf(4).grids={'P4L2','P4L4','P5L2','P5L3','P6L3','P5L0'}; % Near Face Selective (Posterior)
        grids.grpfnf(5).grids={'A7R1','A7L1','A7L2','A6L3','A4R1','A4L1','A4L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Outside Face Selective (All)
        grids.grpbpnf(1).grids={'A7L2','A7L1','A6L3'}; % Inside Bodypart Selective
        grids.grpbpnf(2).grids={'A7R1','A6L2','A6L0','A5L2','A5L1','A5L0','A4L2','A4L1','A4R1',}; % Near Bodypart Selective
        grids.grpbpnf(3).grids={'A7R1','A6L0','A6L2','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Outside BodySelective
        grids.grpobnf(1).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L2'}; % Inside Object Selective
        grids.grpobnf(2).grids={'A4L3','A4L2','A4L1','A4R1','P6L3','P6L1','P7L2'}; % Near Object Selective
        grids.grpobnf(3).grids={'A7R1','A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','P5L0','P6L1','P6L3','P7L2'}; % Outside All Object Selective
        grids.grpplnf(1).grids={'A4L2','A4L1','A4R1','A7R1'}; % Non Category Selective
        grids.grpplnf(2).grids={'A7L2','A7L1','A6L2','A6L0','A5L2','A5L1','A5L0','A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L3','P6L2','P6L1','P7L2'}; % Category Selective
        grids.grpplnf(3).grids={'A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Category Selective
        
        % In Near Far
        grids.grpfnf(1).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Anterior Face Selective
        grids.grpfnf(2).grids={'A7L2','A7L1','A7R1','A4L2','A4L1','A4R1'}; % Near Face Selective (Anterior)
        grids.grpfnf(3).grids={'P6L2','P6L1','P7L2'}; % Posterior Face Selective
        grids.grpfnf(4).grids={'P4L2','P4L4','P5L2','P5L3','P6L3','P5L0'}; % Near Face Selective (Posterior)
        grids.grpfnf(5).grids={'A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5'}; % Outside Face Selective (All)
        grids.grpbpnf(1).grids={'A7L2','A7L1','A6L3'}; % Inside Bodypart Selective
        grids.grpbpnf(2).grids={'A7R1','A6L2','A6L0','A5L2','A5L1','A5L0','A4L2','A4L1','A4R1',}; % Near Bodypart Selective
        grids.grpbpnf(3).grids={'A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Outside BodySelective
        grids.grpobnf(1).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L2'}; % Inside Object Selective
        grids.grpobnf(2).grids={'A4L3','A4L2','A4L1','A4R1','P6L3','P6L1','P7L2'}; % Near Object Selective
        grids.grpobnf(3).grids={'A7R1','A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','P5L0'}; % Outside All Object Selective
        grids.grpplnf(1).grids={'A4L2','A4L1','A4R1','A7R1'}; % Non Category Selective
        grids.grpplnf(2).grids={'A7L2','A7L1','A6L2','A6L0','A5L2','A5L1','A5L0','A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L3','P6L2','P6L1','P7L2'}; % Category Selective
        grids.grpplnf(3).grids={'A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Category Selective
    
    elseif strcmp(monkeyname,'Wiggum')==1,
        % Grid location groups for comparison OLD
        %grids.grp(1).grids={'A6R2','A5R0','A4R3'}; % Bodypart Selective
        %grids.grp(2).grids={'AR0','A2R1','A2R3','A2R5'}; % Face Selective
        %grids.grp(3).grids={'P1R0','P1R3'}; % Bodypart Selective
        %grids.grp(4).grids={'P3R0','P3R2','P5R0'}; % Place Selective
        %grids.grp(5).grids={'P3R0','P3R2','P5R0'}; % Place Selective
        % Grid location groups for comparison NEW
        
        % Current Locations:
        % 'A0R0','A1R0','A2R1','A2R3','A2R5','A3R0','A3R2','A4R3','A5R0','A6R2'
        % 'A7L1','P1R0','P1R3','P3R0','P3R2','P5R0','P6R0'
        %grids.grp(1).grids={'A7L2','A7L1','A7L0','A7R0','A7R1','A6L2','A6L1','A6L0','A6R0','A6R1'}; % Face Selective
        %grids.grp(2).grids={'A7R2','A6R2','A5R3','A5R4'}; % Object Selective
        %grids.grp(3).grids={'A1L1','A1L0','A1R0','A1R1','A1R2','A0L1','A0L0','A0R0','A0R1','A0R2','P1R0','P1L0','P1R1','P1R2','P1R3'}; % Bodypart Selective
        %grids.grp(4).grids={'P3L1','P3L0','P3R0','P3R1','P3R2','P4L0','P4R0','P4R1','P4R2','P4R3'}; % Face Selective
        %grids.grp(5).grids={'P6L1','P6L0','P6R0','P6R1','P6R2','P7L1','P7L0','P7R0','P7R1','P7R2'}; % Object Selective
        
        % updated May 10, 2010
        grids.grp(1).grids={'A7R0','A7L1','A6L1'}; %  Places
        grids.grp(2).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Faces
        grids.grp(3).grids={'P1R0','P1R3'}; % Bodyparts
        grids.grp(4).grids={'P3R0','P3R2'}; % Objects
        grids.grp(5).grids={'P7R0','P7R2'}; %  Face Selective
        
        grids.sampledlocations={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'};
        % Grid location groups for comparison
        grids.grpf(1).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Anterior Face Selective
        grids.grpf(2).grids={'P7R0','P7R2'}; % Posterior Face Selective
        grids.grpf(3).grids={'A0R0','A0R1','A1R0','A2R1','A3R0','P7R0','P7R2'}; % Inside All
        grids.grpf(4).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R2','A2R3','A2R5','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0'}; % Outside Face Selective
        grids.grpbp(1).grids={'P1R0','P1R3'}; % Inside Bodypart Selective
        grids.grpbp(2).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Bodypart Selective
        grids.grpob(1).grids={'P3R0','P3R2'}; % Inside Object Selective
        grids.grpob(2).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Object Selective
        grids.grppl(1).grids={'A7R0','A7L1','A6L1'}; % Inside Place Category Selective
        grids.grppl(2).grids={'A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Place Selective
        
        grids.grpfnf(1).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Anterior Face Selective
        grids.grpfnf(2).grids={'A5R0','A4R3','A3R2','A2R3','A2R5','P1R0','P1R3'}; % Near Face Selective (Anterior)
        grids.grpfnf(3).grids={'P7R0','P7R2'}; % Posterior Face Selective
        grids.grpfnf(4).grids={'P4R1','P5R0','P6R0','P7R0','P7R2'}; % Near Face Selective (Posterior)
        grids.grpfnf(5).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R2','A2R3','A2R5','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0'}; % Outside all Face Selective
        grids.grpbpnf(1).grids={'P1R0','P1R3'}; % Inside Bodypart Selective
        grids.grpbpnf(2).grids={'A2R1','A1R0','A0R0','A0R1','P3R0','P3R2'}; % Near Bodypart Selective
        grids.grpbpnf(3).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Bodypart Selective
        grids.grpobnf(1).grids={'P3R0','P3R2'}; % Inside Object Selective
        grids.grpobnf(2).grids={'A0R0','A0R1','P1R0','P1R3','P4R1','P5R0'}; % Near Object Selective
        grids.grpobnf(3).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside Object Selective
        grids.grpplnf(1).grids={'A7R0','A7L1','A6L1'}; % Inside Place Category Selective
        grids.grpplnf(2).grids={'A6R2','A5R0','A4R3'}; % Near Place Category Selective
        grids.grpplnf(3).grids={'A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Category Selective
        
        % In Near Far
        grids.grpfnf(1).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Anterior Face Selective
        grids.grpfnf(2).grids={'A5R0','A4R3','A3R2','A2R3','A2R5','P1R0','P1R3'}; % Near Face Selective (Anterior)
        grids.grpfnf(3).grids={'P7R0','P7R2'}; % Posterior Face Selective
        grids.grpfnf(4).grids={'P4R1','P5R0','P6R0'}; % Near Face Selective (Posterior)
        grids.grpfnf(5).grids={'A7L1','A7R0','A6L1','A6R2','P3R0','P3R2'}; % Outside all Face Selective
        grids.grpbpnf(1).grids={'P1R0','P1R3'}; % Inside Bodypart Selective
        grids.grpbpnf(2).grids={'A2R1','A1R0','A0R0','A0R1','P3R0','P3R2'}; % Near Bodypart Selective
        grids.grpbpnf(3).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R3','A2R5','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Bodypart Selective
        grids.grpobnf(1).grids={'P3R0','P3R2'}; % Inside Object Selective
        grids.grpobnf(2).grids={'A0R0','A0R1','P1R0','P1R3','P4R1','P5R0'}; % Near Object Selective
        grids.grpobnf(3).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','P6R0','P7R0','P7R2'}; % Outside Object Selective
        grids.grpplnf(1).grids={'A7R0','A7L1','A6L1'}; % Inside Place Category Selective
        grids.grpplnf(2).grids={'A6R2','A5R0','A4R3'}; % Near Place Category Selective
        grids.grpplnf(3).grids={'A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Category Selective
    end
end
toc; close (h)
fprintf(['\n  (Loaded data from ',num2str(size(UnitData.plxname,1)),' neurons)\n']);

UnitData=orderfields(UnitData); grids=orderfields(grids);

save([lsnconfig.datadir,'td500data_',monkeyname,'.mat'],'UnitData','grids','monkeyname');
return