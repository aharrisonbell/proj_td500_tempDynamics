function START_C3_estimateRDMs
% >>> needs debugging

% Version 3: single-neuron linear-discriminant contrast

% Author: Marieke Mur; last edit 06-10-2015


%% preparation
clear; close all;

resultsPath='/imaging/mm07/mITReprDynamics/analysis/results';
resultsPath_vSpec=fullfile(resultsPath,'RDMs_singleNeuron');
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));


%% control variables
subjStr={'Stew' 'Wigg'};
nStimuli=100;
monitor=1;

RSA_timepoints=-100:700; % ms
RSA_tempSm_window=21; % ms (should be odd)


%% estimate RDMs
for subjectI=1:numel(subjStr)
    % load extracted patterns and trial information
    load(fullfile(resultsPath,['trialInfo_',subjStr{subjectI}]),'trialInfo_neuron__trial_stimCat');
    load(fullfile(resultsPath,['SPDsingleTrial_',subjStr{subjectI}]));
    load(fullfile(resultsPath,['SPDexemplar_',subjStr{subjectI}]));
    
    % estimate RDM at each timepoint (per neuron!)
    nNeurons=numel(spd_neuron_4LDC__trial_timepoints);
    RSA_timepoints_recod=RSA_timepoints+1000; % spd data ranges from 0 to 5000 ms; 1000 = stimulus onset
    windowWidthOneSide=floor(RSA_tempSm_window/2);
    RDMweights_neuron=zeros(nStimuli,nStimuli,nNeurons);
    for neuronI=1:nNeurons
        tic
        trialInfo_cNeuron__trial_stimCat=trialInfo_neuron_4LDC__trial_stimCat{neuronI};
        spd_cNeuron__trial_timepoints=spd_neuron_4LDC__trial_timepoints{neuronI};
        % average single-trial spike-density functions within 20-ms windows
        for timepointI=1:numel(RSA_timepoints)
            windowStart=RSA_timepoints_recod(timepointI)-windowWidthOneSide; windowEnd=RSA_timepoints_recod(timepointI)+windowWidthOneSide;
            spd_cNeuron_tempSm__trial_timepoints(:,timepointI)=mean(spd_cNeuron__trial_timepoints(:,windowStart:windowEnd),2);
        end
        % compute residuals for each stimulus (single-trial spd MINUS trial-avg spd)
        for stimulusI=1:nStimuli
            trialIs_cStim=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimulusI);
            spd_cNeuron_tempSm_cStim__trial_timepoints=spd_cNeuron_tempSm__trial_timepoints(trialIs_cStim,:);
            spd_cNeuron_tempSm_cStim_trialAvg__timepoints=mean(spd_cNeuron_tempSm_cStim__trial_timepoints,1);
            spdResid_cNeuron_tempSm__stim__trial_timepoints{stimulusI}=spd_cNeuron_tempSm_cStim__trial_timepoints-repmat(spd_cNeuron_tempSm_cStim_trialAvg__timepoints,[numel(trialIs_cStim) 1]);
        end
        % for each stimulus pair, estimate the linear-discriminant distance (LDC)
        RDMs_timepoints=zeros(nStimuli,nStimuli,numel(RSA_timepoints));
        RDM_ltvMask=logical(tril(ones(nStimuli,nStimuli),-1));
        dissimIs=find(RDM_ltvMask);
        for stimPairI=1:numel(dissimIs)
            [stimIa,stimIb]=ind2sub([nStimuli nStimuli],dissimIs(stimPairI));
            trialIsA=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimIa); nTrialsA=numel(trialIsA);
            trialIsB=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimIb); nTrialsB=numel(trialIsB);
            nFolds=nTrialsA*nTrialsB;
            if numel(trialIsA)>2 && numel(trialIsB)>2
                [trialIsI_test_stimA,trialIsI_test_stimB]=ind2sub([nTrialsA nTrialsB],1:nFolds);
                for foldI=1:nFolds
                    % extract test and training sets (for each split: 1 trial for testing, the rest for training)
                    trialI_test_stimA=trialIsA(trialIsI_test_stimA(foldI));
                    trialI_test_stimB=trialIsB(trialIsI_test_stimB(foldI));
                    resp_a_test__timepoints=spd_cNeuron_tempSm__trial_timepoints(trialI_test_stimA,:);
                    resp_a_train__timepoints=mean(spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsA,trialI_test_stimA),:));
                    resp_b_test__timepoints=spd_cNeuron_tempSm__trial_timepoints(trialI_test_stimB,:);
                    resp_b_train__timepoints=mean(spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsB,trialI_test_stimB),:));
                    % compute response variance for current neuron and fold
                    residuals_cNeuron_stimA__trial_timepoints=spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsA,trialI_test_stimA),:)-repmat(resp_a_train__timepoints,[nTrialsA-1 1]);
                    residuals_cNeuron_stimB__trial_timepoints=spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsB,trialI_test_stimB),:)-repmat(resp_b_train__timepoints,[nTrialsB-1 1]);
                    residuals_cNeuron_otherStim__trial_timepoints=cell2mat([spdResid_cNeuron_tempSm__stim__trial_timepoints(setxor(1:nStimuli,[stimIa stimIb]))]');
                    residuals_cNeuron__trial_timepoints=cat(1,residuals_cNeuron_stimA__trial_timepoints,residuals_cNeuron_stimB__trial_timepoints,residuals_cNeuron_otherStim__trial_timepoints);
                    residualsSSQ_cNeuron__timepoints=sum(residuals_cNeuron__trial_timepoints.^2);
                    df=size(residuals_cNeuron__trial_timepoints,1)-nStimuli;
                    var_cNeuron__timepoints=residualsSSQ_cNeuron__timepoints./df;
                    var__neuron__fold_timepoints{neuronI}(foldI,:)=var_cNeuron__timepoints;
                    % compute LDC
                    LDC_folds_timepoints(foldI,:)=(resp_b_train__timepoints-resp_a_train__timepoints)./var_cNeuron__timepoints.*(resp_b_test__timepoints-resp_a_test__timepoints);
                end % foldI
                RDMs_timepoints(stimIa,stimIb,:)=mean(LDC_folds_timepoints,1);
                RDMweights_neuron(stimIa,stimIb,neuronI)=nFolds;
                clear LDC_folds_timepoints
            else
                RDMs_timepoints(stimIa,stimIb,:)=nan(1,numel(RSA_timepoints));
                RDMweights_neuron(stimIa,stimIb,neuronI)=0;
            end
        end % stimPairI
        for timepointI=1:size(RDMs_timepoints,3)
            cRDM=RDMs_timepoints(:,:,timepointI);
            cRDM_ltv=cRDM(RDM_ltvMask);
            cRDM=squareform(cRDM_ltv);
            RDMs_neuron_timepoints(:,:,neuronI,timepointI)=cRDM;
        end
        clear spd_cNeuron_tempSm__trial_timepoints spdResid_cNeuron_tempSm__stim__trial_timepoints
        toc
    end % neuronI
    RDMs_timepoints=(RDMs_neuron_timepoints.*repmat(RDMweights_neuron,[1 1 1 size(RDMs_timepoints,3)]))./repmat(sum(RDMweights_neuron,3),[1 1 nNeurons size(RDMs_timepoints,3)]);
    save(fullfile(resultsPath_vSpec,['SPDsingleTrial_4LDC_tempSm_',subjStr{subjectI}]),'spd_cNeuron_tempSm__trial_timepoints','var__neuron__fold_timepoints','-v7.3');
    save(fullfile(resultsPath_vSpec,['RDMs_',subjStr{subjectI}]),'RDMweights_neuron','RDMs_neuron_timepoints','RDMs_timepoints','-v7.3');
    
    clear trialSelect_LOG__neuron spd_neuron_4LDC__trial_timepoints trialInfo_neuron_4LDC__trial_stimCat
    clear spd_test__neuron_stim_split_timepoints spd_train__neuron_stim_split_timepoints spdResid_train_cNeuron__stim_split__trial_timepoints trialIs_test__neuron_stim trialIs_train__neuron_stim var__neuron_split_timepoints
    clear spd_test_tempSm__neuron_stim_split_timepoints spd_train_tempSm__neuron_stim_split_timepoints
    clear RDMs_timepoints movieFrames
    clear RDMs_train1__timepoints RDMs_train2__timepoints RDMs_test1__timepoints RDMs_test2__timepoints
end % subjectI





