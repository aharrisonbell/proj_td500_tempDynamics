function START_C2_estimateRDMs
% >>> needs debugging

% Version 2: linear-discriminant contrast

% This function extracts the single-trial spike density functions for each
% subject, neuron, and stimulus. The extracted spd functions serve as a
% basis for computing a stimulus x stimulus RDM at each timepoint. The RDMs
% are used to investigate the representational dynamics of object responses
% in primate IT.

% There are different versions, each of which computes a (slightly)
% different distance measure. All versions use visually-responsive neurons
% for which we have information about the neuron's location along the AP
% axis (about 500 neurons in each subject) or a subset of these neurons.
% Furthermore, all versions use a sliding-window approach to temporally
% smooth the spd functions before computing pattern distances.

% Version 1: correlation distance. Stimulus response patterns are based on
%               trial-average spd functions. All visually-responsive
%               AP-index neurons participate. Inference on the presence of
%               category information is based on stimulus-label
%               randomisation tests and FDR correction for multiple
%               comparisons (time points).
% Version 2: linear-discriminant contrast. Data is split into training and
%               test sets before computing response patterns. Only neurons
%               that have > 2 trials for each stimulus participate.
%               Inference same as version 1.
% Version 3: linear-discriminant contrast. Data is split into training and
%               test sets before computing response patterns. All
%               visually-responsive AP-index neurons participate. RDMs are
%               computed for each neuron separately, and combined using
%               weights for each stimulus pair and neuron. Weights
%               reflect how much data each neuron contributes to a
%               particular response-pattern distance (i.e. they are the
%               number of folds for that particular stimulus comparison).
%               Inference on the presence of exemplar and category
%               information is based on bootstrap-resampling of the neurons
%               and FDR correction for multiple comparisons (time points).

% Note: the linear-discriminant contrast is the cross-validated Mahalanobis
% distance between two stimulus response patterns. We assume that the
% neurons are independent, so the covariance matrix reduces to the diagonal
% (variances only).

% Details version 2:
% (1) compute the LDC using only neurons that have >2 trials for each of
%     the 100 stimuli (needed for noise normalisation and cross-validation)
% (2) use the same number of folds (i.e. 2) across neurons
% (3) for each fold and stimulus, use one trial as test data and the
%     average of the remaining trials as training data. This implies that
%     the amount of training data might differ across neurons.
% (4) compute each neuron's response variance before applying temporal
%     smoothing

% Author: Marieke Mur; last edit: 06-10-2015


%% preparation
clear; close all;

resultsPath='/imaging/mm07/mITReprDynamics/analysis/results';
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));


%% control variables
subjStr={'Stew' 'Wigg'};
nStimuli=100;
monitor=0;

RSA_timepoints=-100:700; % ms
RSA_tempSm_window=21; % ms (should be odd)
RSA_nDataSplits=2; % not a real variable (needs to be set at 2)


%% estimate RDMs
for subjectI=1:numel(subjStr)
    % load extracted patterns and trial information
    load(fullfile(resultsPath,['trialInfo_',subjStr{subjectI}]),'trialInfo_neuron__trial_stimCat');
    load(fullfile(resultsPath,['SPDsingleTrial_',subjStr{subjectI}]));
    load(fullfile(resultsPath,['SPDexemplar_',subjStr{subjectI}]));
        
    % how many neurons contain 3 or more trials for each stimulus? (these neurons can be used for LDC estimation)
    min3ValidTrials_LOG__neuron_stim=false(size(nValidTrials__neuron_stim));
    min3ValidTrials_LOG__neuron_stim(nValidTrials__neuron_stim>2)=1;
    min3ValidTrials_LOG_summedAcrossStim__neuron=sum(min3ValidTrials_LOG__neuron_stim,2);
    cvAmenableNeuronIs=find(min3ValidTrials_LOG_summedAcrossStim__neuron==100);
    
    % inspect distribution of valid trials over neurons (for each stimulus)
    if monitor
        % show number of valid trials over neurons and stimuli
        figI=199; clf;
        figure(figI);
        subplot(3,1,1); imagesc(nValidTrials__neuron_stim); colormap(gca, RDMcolormap); colorbar;
        xlabel('stimuli'); ylabel('neurons'); title('number of valid trials (no threshold)')
        subplot(3,1,2); imagesc(min3ValidTrials_LOG__neuron_stim); colormap(gca, RDMcolormap); colorbar;
        xlabel('stimuli'); ylabel('neurons'); title({'number of valid trials (thresholded)','\fontsize{8}[0 indicates < 3, 1 indicates >= 3]'});
        subplot(3,1,3); plot(histc(min3ValidTrials_LOG_summedAcrossStim__neuron,1:nStimuli));
        xlabel('number of stimuli for which there are >= 3 trials'); ylabel('number of neurons'); title('do most neurons have >= 3 trials for all 100 stimuli?');
        description(1)={'\fontsize{12}number of valid trials'};
        addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),figI);
        % estimate distributions
        for stimulusI=1:nStimuli
            nValidTrials__neuron_cStim=nValidTrials__neuron_stim(:,stimulusI);
            table_cStim=tabulate(nValidTrials__neuron_cStim); % col 1 = n valid trials; col 2 = how many neurons; col 3 = percentage of neurons
            tables_stim{stimulusI}=table_cStim;
            minnValidTrials_stim(stimulusI)=min(table_cStim(:,1));
            maxnValidTrials_stim(stimulusI)=max(table_cStim(:,1));
            meannValidTrials_stim(stimulusI)=mean(nValidTrials__neuron_cStim);
        end
        maxnValidTrials=max(maxnValidTrials_stim);
        % plot min, mean, and max number of valid trials
        figI=200; clf;
        figure(figI); plot(minnValidTrials_stim); hold on;
        plot(meannValidTrials_stim,'k'); hold on;
        plot(maxnValidTrials_stim,'r');
        legend('min','mean','max');
        xlabel('stimulus number');
        ylabel('number of valid trials');
        description(1)={'\fontsize{12}number of valid trials (across neurons)'};
        addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),figI);
        % plot histograms
        figI=201; clf;
        anchorCols=[.9 .9 .9; 0 0 1];
        cols=colorScale(anchorCols,100,0);
        for stimulusI=1:nStimuli
            table_cStim=tables_stim{stimulusI};
            histogram_cStim=zeros(maxnValidTrials+1,1);
            histogram_cStim(table_cStim(:,1)+1)=table_cStim(:,2);
            figure(figI); plot(histogram_cStim,'Color',cols(stimulusI,:)); hold on;
        end
        set(gca,'XTick',1:10:maxnValidTrials+1);
        set(gca,'XTickLabel',{'0','10','20','30','40','50','60','70'});
        xlabel('number of valid trials');
        ylabel('number of neurons');
        text(0.62*maxnValidTrials,20,'each line is a stimulus');
        description(1)={'\fontsize{12}distribution of number of valid trials over neurons'};
        description(2)={['\fontsize{10}visually-responsive neurons (n=',num2str(size(nValidTrials__neuron_stim,1)),')']};
        description(3)={'\fontsize{10}all stimuli'};
        addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),figI);
        % save
        save(fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),'min3ValidTrials_LOG__neuron_stim','min3ValidTrials_LOG_summedAcrossStim__neuron','cvAmenableNeuronIs',...
            'tables_stim','minnValidTrials_stim','maxnValidTrials_stim','meannValidTrials_stim');
        clear tables_stim minnValidTrials_stim maxnValidTrials_stim meannValidTrials_stim description
    end
    
    % select visually-responsive neurons for which we have trial information (stimulus and category) and AP index
    spd_visRespNeuron__trial_timepoints=spd_neuron__trial_timepoints(neuronSelect_LOG);
    trialInfo_visRespNeuron__trial_stimCat=trialInfo_neuron__trial_stimCat(neuronSelect_LOG);
    % furthermore, select neurons that have >= 3 trials for each stimulus, and exclude trials with zero mean and/or zero variance
    for neuronIsI=1:numel(cvAmenableNeuronIs)
        neuronI=cvAmenableNeuronIs(neuronIsI);
        spd_cNeuron__trial_timepoints=spd_visRespNeuron__trial_timepoints{neuronI};
        trialInfo_cNeuron__trial_stimCat=trialInfo_visRespNeuron__trial_stimCat{neuronI};
        spd_cNeuron__trial_timepointAvg=mean(spd_cNeuron__trial_timepoints,2);
        spd_cNeuron__trial_timepointVar=var(spd_cNeuron__trial_timepoints,0,2);
        avg_LOG=false(size(spd_cNeuron__trial_timepoints,1),1); avg_LOG(spd_cNeuron__trial_timepointAvg>0)=true;
        var_LOG=false(size(spd_cNeuron__trial_timepoints,1),1); var_LOG(spd_cNeuron__trial_timepointVar>0)=true;
        trialSelect_LOG=logical(floor((avg_LOG+var_LOG)/2));
        trialSelect_LOG__neuron{neuronIsI}=trialSelect_LOG;
        spd_neuron_4LDC__trial_timepoints{neuronIsI}=spd_cNeuron__trial_timepoints(trialSelect_LOG,:);
        trialInfo_neuron_4LDC__trial_stimCat{neuronIsI}=trialInfo_cNeuron__trial_stimCat(trialSelect_LOG,:);
    end
    save(fullfile(resultsPath,['SPDsingleTrial_4LDC_',subjStr{subjectI}]),'neuronSelect_LOG','cvAmenableNeuronIs','trialSelect_LOG__neuron',...
        'spd_neuron_4LDC__trial_timepoints','trialInfo_neuron_4LDC__trial_stimCat','-v7.3');
    
    % split data in test and training sets (for each split: 1 trial for testing, the rest for training)
    nNeurons=numel(spd_neuron_4LDC__trial_timepoints);
    for neuronI=1:nNeurons
        trialInfo_cNeuron__trial_stimCat=trialInfo_neuron_4LDC__trial_stimCat{neuronI};
        spd_cNeuron__trial_timepoints=spd_neuron_4LDC__trial_timepoints{neuronI};
        for stimulusI=1:nStimuli
            trialIs=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimulusI);
            trialIs_rand=trialIs(randperm(numel(trialIs)));
            for splitI=1:RSA_nDataSplits
                % assign trials to test and training sets
                trialI_test=trialIs_rand(splitI);
                trialIs_train=setxor(trialIs_rand,trialI_test);
                % extract test and training data
                spd_test__neuron_stim_split_timepoints(neuronI,stimulusI,splitI,:)=spd_cNeuron__trial_timepoints(trialI_test,:);
                spd_train_cNeuron_cStim_cSplit__trial_timepoints=spd_cNeuron__trial_timepoints(trialIs_train,:);
                spd_train_cNeuron_cStim_cSplit_trialAvg__timepoints=mean(spd_train_cNeuron_cStim_cSplit__trial_timepoints,1);
                spd_train__neuron_stim_split_timepoints(neuronI,stimulusI,splitI,:)=spd_train_cNeuron_cStim_cSplit_trialAvg__timepoints;
                % compute residuals for (co)variance matrix
                spdResid_train_cNeuron_cStim_cSplit__trial_timepoints=spd_train_cNeuron_cStim_cSplit__trial_timepoints-repmat(spd_train_cNeuron_cStim_cSplit_trialAvg__timepoints,[numel(trialIs_train) 1]);
                % store to pass on
                trialIs_test__trial_split(splitI)=trialI_test;
                trialIs_train__trial_split(:,splitI)=trialIs_train;
                spdResid_train_cNeuron__stim_split__trial_timepoints{stimulusI,splitI}=spdResid_train_cNeuron_cStim_cSplit__trial_timepoints;
            end % splitI
            trialIs_test__neuron_stim{neuronI,stimulusI}=trialIs_test__trial_split;
            trialIs_train__neuron_stim{neuronI,stimulusI}=trialIs_train__trial_split;
            clear trialIs_test__trial_split trialIs_train__trial_split
        end % stimulusI
        % compute response variance for current neuron (for each split)
        for splitI=1:RSA_nDataSplits
            residuals_cNeuron=cell2mat(spdResid_train_cNeuron__stim_split__trial_timepoints(:,splitI));
            residualsSSQ_cNeuron=sum(residuals_cNeuron.^2);
            df=size(residuals_cNeuron,1)-nStimuli;
            var_cNeuron=residualsSSQ_cNeuron./df;
            var__neuron_split_timepoints(neuronI,splitI,:)=var_cNeuron;
        end
    end % neuronI
    save(fullfile(resultsPath,['SPDtrainTestVar_4LDC_',subjStr{subjectI}]),'trialIs_test__neuron_stim','trialIs_train__neuron_stim',...
        'spd_test__neuron_stim_split_timepoints','spd_train__neuron_stim_split_timepoints','var__neuron_split_timepoints');    
    
    % estimate RDM at each timepoint
    RSA_timepoints_recod=RSA_timepoints+1000; % spd data ranges from 0 to 5000 ms; 1000 = stimulus onset
    windowWidthOneSide=floor(RSA_tempSm_window/2);
    RDMs_timepoints=zeros(nStimuli,nStimuli,numel(RSA_timepoints));
    for timepointI=1:numel(RSA_timepoints)
        % average single-trial spike-density functions within 20-ms windows
        windowStart=RSA_timepoints_recod(timepointI)-windowWidthOneSide; windowEnd=RSA_timepoints_recod(timepointI)+windowWidthOneSide;
        spd_test_tempSm__neuron_stim_split_timepoints(:,:,:,timepointI)=mean(spd_test__neuron_stim_split_timepoints(:,:,:,windowStart:windowEnd),4);
        spd_train_tempSm__neuron_stim_split_timepoints(:,:,:,timepointI)=mean(spd_train__neuron_stim_split_timepoints(:,:,:,windowStart:windowEnd),4);
        % for each stimulus pair, estimate the linear-discriminant distance (LDC)
        RDM=zeros(nStimuli,nStimuli);
        RDM_ltvMask=logical(tril(ones(nStimuli,nStimuli),-1));
        dissimIs=find(RDM_ltvMask);
        for stimPairI=1:numel(dissimIs)
            [stimIa,stimIb]=ind2sub(size(RDM),dissimIs(stimPairI));
            for foldI=1:RSA_nDataSplits
                pattern_a_test=spd_test_tempSm__neuron_stim_split_timepoints(:,stimIa,foldI,timepointI);
                pattern_a_train=spd_train_tempSm__neuron_stim_split_timepoints(:,stimIa,foldI,timepointI);
                pattern_b_test=spd_test_tempSm__neuron_stim_split_timepoints(:,stimIb,foldI,timepointI);
                pattern_b_train=spd_train_tempSm__neuron_stim_split_timepoints(:,stimIb,foldI,timepointI);
                var_train=var__neuron_split_timepoints(:,foldI,RSA_timepoints_recod(timepointI));
                LDC_folds(foldI)=[pattern_b_train-pattern_a_train]'./var_train'*(pattern_b_test-pattern_a_test);
            end % foldI
            RDM(stimIa,stimIb)=mean(LDC_folds);
            clear LDC_folds
        end % stimPairI
        RDM_ltv=RDM(RDM_ltvMask);
        RDM=squareform(RDM_ltv);
        RDMs_timepoints(:,:,timepointI)=RDM;
        % check: correlation-distance RDMs
        patterns_train1=[spd_train_tempSm__neuron_stim_split_timepoints(:,:,1,timepointI)]';
        patterns_train2=[spd_train_tempSm__neuron_stim_split_timepoints(:,:,2,timepointI)]';
        patterns_test1=[spd_test_tempSm__neuron_stim_split_timepoints(:,:,1,timepointI)]';
        patterns_test2=[spd_test_tempSm__neuron_stim_split_timepoints(:,:,2,timepointI)]';
        RDM_ltv_train1=pdist(patterns_train1,'correlation'); RDM_train1=squareform(RDM_ltv_train1);
        RDMs_train1__timepoints(:,:,timepointI)=RDM_train1;
        RDM_ltv_train2=pdist(patterns_train2,'correlation'); RDM_train2=squareform(RDM_ltv_train2);
        RDMs_train2__timepoints(:,:,timepointI)=RDM_train2;
        RDM_ltv_test1=pdist(patterns_test1,'correlation'); RDM_test1=squareform(RDM_ltv_test1);
        RDMs_test1__timepoints(:,:,timepointI)=RDM_test1;
        RDM_ltv_test2=pdist(patterns_test2,'correlation'); RDM_test2=squareform(RDM_ltv_test2);
        RDMs_test2__timepoints(:,:,timepointI)=RDM_test2;
    end % timepointI
    
    % save
    save(fullfile(resultsPath,['SPDtrainTestVar_4LDC_tempSm_',subjStr{subjectI}]),...
        'spd_test_tempSm__neuron_stim_split_timepoints','spd_train_tempSm__neuron_stim_split_timepoints');
    save(fullfile(resultsPath,['RDMs_ldc_',subjStr{subjectI}]),'RDMs_timepoints','RSA_timepoints');
    save(fullfile(resultsPath,['RDMs_cor_',subjStr{subjectI}]),'RDMs_train1__timepoints','RDMs_train2__timepoints','RDMs_test1__timepoints','RDMs_test2__timepoints','RSA_timepoints');
    
    
    clear trialSelect_LOG__neuron spd_neuron_4LDC__trial_timepoints trialInfo_neuron_4LDC__trial_stimCat
    clear spd_test__neuron_stim_split_timepoints spd_train__neuron_stim_split_timepoints spdResid_train_cNeuron__stim_split__trial_timepoints trialIs_test__neuron_stim trialIs_train__neuron_stim var__neuron_split_timepoints
    clear spd_test_tempSm__neuron_stim_split_timepoints spd_train_tempSm__neuron_stim_split_timepoints
    clear RDMs_timepoints
    clear RDMs_train1__timepoints RDMs_train2__timepoints RDMs_test1__timepoints RDMs_test2__timepoints
end % subjectI



