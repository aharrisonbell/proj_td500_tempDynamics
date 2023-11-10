function START_C4a_estimateRDMs(subjectIs,neuronIs)

% Version 4: single-neuron spike rate distance (SRD) = the cross-validated 
%            difference in spike rate between two stimuli. Subsequent
%            analyses average the SRDs across neurons (weighted by SNR) 
%            (see START_D4).
%         a: each neuron can be ran separately, but user can also give in a
%            list of neurons 

% Author: Marieke Mur; last edit 26-04-2016; SRD developed in collaboration
% with Nikolaus Kriegeskorte


%% preparation
close all;

resultsPath='/imaging/mm07/mITReprDynamics/analysis/results';
resultsPath_singleNeuron=fullfile(resultsPath,'singleNeuron');
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));
try mkdir(resultsPath_singleNeuron); end


%% control variables
subjStr={'Stew' 'Wigg'};
nStimuli=100;

RSA_timepoints=-100:700; % ms
RSA_tempSm_window=21; % ms (should be odd)
nTimepoints=numel(RSA_timepoints);

save(fullfile(resultsPath,'START_C4a_variables'));


%% estimate RDMs
for subjectI=subjectIs
    % load extracted patterns and trial information    
    load(fullfile(resultsPath,['SPDsingleTrial_4RSA_',subjStr{subjectI}]),'spd_4RSA__neuron__trial_timepoints');
    load(fullfile(resultsPath,['trialInfo_4RSA_',subjStr{subjectI}]),'trialInfo_4RSA__neuron__trial_stimCat');
    
    % estimate RDM at each timepoint (per neuron!)
    RSA_timepoints_recod=RSA_timepoints+1000; % spd data ranges from 0 to 5000 ms; 1000 = stimulus onset
    windowWidthOneSide=floor(RSA_tempSm_window/2);            
    for neuronI=neuronIs                
        spd_cNeuron__trial_timepoints=spd_4RSA__neuron__trial_timepoints{neuronI};        
        if numel(neuronIs)==1, clear spd_4RSA__neuron__trial_timepoints; end % to free up memory
        trialInfo_cNeuron__trial_stimCat=trialInfo_4RSA__neuron__trial_stimCat{neuronI};
        if numel(neuronIs)==1, clear trialInfo_4RSA__neuron__trial_stimCat; end % to free up memory
        
        % average single-trial spike-density functions within a temporal window
        spd_cNeuron_tempSm__trial_timepoints=nan(size(spd_cNeuron__trial_timepoints,1),nTimepoints);
        for timepointI=1:nTimepoints
            windowStart=RSA_timepoints_recod(timepointI)-windowWidthOneSide; windowEnd=RSA_timepoints_recod(timepointI)+windowWidthOneSide;
            spd_cNeuron_tempSm__trial_timepoints(:,timepointI)=mean(spd_cNeuron__trial_timepoints(:,windowStart:windowEnd),2);
        end
        clear spd_cNeuron__trial_timepoints % to free up memory

        % for each stimulus pair, estimate the spike rate distance (srd)
        RDMweights=zeros(nStimuli,nStimuli);
        RDMs_srd=zeros(nStimuli,nStimuli,nTimepoints);        
        RDM_ltvMask=logical(tril(ones(nStimuli,nStimuli),-1));
        dissimIs=find(RDM_ltvMask);
        for stimPairI=1:numel(dissimIs)
            [stimIa,stimIb]=ind2sub([nStimuli nStimuli],dissimIs(stimPairI));
            trialIsA=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimIa); nTrialsA=numel(trialIsA);
            trialIsB=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimIb); nTrialsB=numel(trialIsB);
            nFolds=nTrialsA*nTrialsB;
            if numel(trialIsA)>1 && numel(trialIsB)>1
                srd__folds_timepoints=nan(nFolds,nTimepoints);
                [trialIsI_test_stimA,trialIsI_test_stimB]=ind2sub([nTrialsA nTrialsB],1:nFolds);                
                for foldI=1:nFolds                      
                    % extract test and training sets (for each split: 1 trial for testing, the rest for training)
                    trialI_test_stimA=trialIsA(trialIsI_test_stimA(foldI));
                    trialI_test_stimB=trialIsB(trialIsI_test_stimB(foldI));
                    resp_a_test__timepoints=spd_cNeuron_tempSm__trial_timepoints(trialI_test_stimA,:);
                    resp_a_train__timepoints=mean(spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsA,trialI_test_stimA),:),1);
                    resp_b_test__timepoints=spd_cNeuron_tempSm__trial_timepoints(trialI_test_stimB,:);
                    resp_b_train__timepoints=mean(spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsB,trialI_test_stimB),:),1);
                    % compute the spike rate distance
                    train_diff=resp_b_train__timepoints-resp_a_train__timepoints;
                    test_diff=resp_b_test__timepoints-resp_a_test__timepoints;
                    srd(1,:)=sign(train_diff).*test_diff;
                    srd(2,:)=sign(test_diff).*train_diff;
                    srd__folds_timepoints(foldI,:)=mean(srd,1);
                    clear srd
                end % foldI
                RDMs_srd(stimIa,stimIb,:)=mean(srd__folds_timepoints,1);
                RDMweights(stimIa,stimIb)=nFolds;
            else
                RDMs_srd(stimIa,stimIb,:)=nan(1,nTimepoints);
                RDMweights(stimIa,stimIb)=0;
            end
        end % stimPairI      
        
        % mirror ltv to utv
        for timepointI=1:size(RDMs_srd,3)
            cRDM=RDMs_srd(:,:,timepointI);
            cRDM_ltv=cRDM(RDM_ltvMask);
            cRDM=squareform(cRDM_ltv);
            RDMs_srd(:,:,timepointI)=cRDM;
        end        
        RDMweights_ltv=RDMweights(RDM_ltvMask);
        RDMweights=squareform(RDMweights_ltv);
        
        % save and clear
        save(fullfile(resultsPath_singleNeuron,['SPDsingleTrial_4RSA_tempSm_',subjStr{subjectI},'_neuron',num2str(neuronI)]),'spd_cNeuron_tempSm__trial_timepoints');
        save(fullfile(resultsPath_singleNeuron,['RDMs_srd_',subjStr{subjectI},'_neuron',num2str(neuronI)]),'RDMs_srd','RDMweights');                
        clear spd_cNeuron_tempSm__trial_timepoints RDMs_srd RDMweights         
    end % neuronI              
end % subjectI



