function START_C3a_estimateRDMs(subjectIs,neuronIs)

% Version 3: single-neuron linear-discriminant contrast
%         a: each neuron can be ran separately, but user can also give in a
%            list of neurons 

% Author: Marieke Mur; last edit 07-10-2015


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

save(fullfile(resultsPath,'START_C3a_variables'));


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
        
        % compute residuals for each stimulus (single-trial spd MINUS trial-avg spd)
        spdResid_cNeuron_tempSm__stim__trial_timepoints=cell(nStimuli,1);
        stim_LOG=true(nStimuli,1);
        for stimulusI=1:nStimuli
            trialIs_cStim=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimulusI);
            if isempty(trialIs_cStim), stim_LOG(stimulusI)=false; end
            spd_cNeuron_tempSm_cStim__trial_timepoints=spd_cNeuron_tempSm__trial_timepoints(trialIs_cStim,:);
            spd_cNeuron_tempSm_cStim_trialAvg__timepoints=mean(spd_cNeuron_tempSm_cStim__trial_timepoints,1);
            spdResid_cNeuron_tempSm__stim__trial_timepoints{stimulusI}=spd_cNeuron_tempSm_cStim__trial_timepoints-repmat(spd_cNeuron_tempSm_cStim_trialAvg__timepoints,[numel(trialIs_cStim) 1]);
        end
        
        % for each stimulus pair, estimate the linear-discriminant contrast (ldc)
        RDMweights=zeros(nStimuli,nStimuli);
        RDMs_ldc=zeros(nStimuli,nStimuli,nTimepoints);        
        RDM_ltvMask=logical(tril(ones(nStimuli,nStimuli),-1));
        dissimIs=find(RDM_ltvMask);
        validStimPairI=0; % >>>
        for stimPairI=1:numel(dissimIs)
            [stimIa,stimIb]=ind2sub([nStimuli nStimuli],dissimIs(stimPairI));
            trialIsA=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimIa); nTrialsA=numel(trialIsA);
            trialIsB=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimIb); nTrialsB=numel(trialIsB);
            nFolds=nTrialsA*nTrialsB;
            if numel(trialIsA)>1 && numel(trialIsB)>1
                validStimPairI=validStimPairI+1; % >>>
                ldc__folds_timepoints=nan(nFolds,nTimepoints);
                [trialIsI_test_stimA,trialIsI_test_stimB]=ind2sub([nTrialsA nTrialsB],1:nFolds);                
                for foldI=1:nFolds                      
                    % extract test and training sets (for each split: 1 trial for testing, the rest for training)
                    trialI_test_stimA=trialIsA(trialIsI_test_stimA(foldI));
                    trialI_test_stimB=trialIsB(trialIsI_test_stimB(foldI));
                    resp_a_test__timepoints=spd_cNeuron_tempSm__trial_timepoints(trialI_test_stimA,:);
                    resp_a_train__timepoints=mean(spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsA,trialI_test_stimA),:),1);
                    resp_b_test__timepoints=spd_cNeuron_tempSm__trial_timepoints(trialI_test_stimB,:);
                    resp_b_train__timepoints=mean(spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsB,trialI_test_stimB),:),1);
                    % compute response variance for current neuron and fold
                    residuals_cNeuron_stimA__trial_timepoints=spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsA,trialI_test_stimA),:)-repmat(resp_a_train__timepoints,[nTrialsA-1 1]);
                    residuals_cNeuron_stimB__trial_timepoints=spd_cNeuron_tempSm__trial_timepoints(setxor(trialIsB,trialI_test_stimB),:)-repmat(resp_b_train__timepoints,[nTrialsB-1 1]);
                    residuals_cNeuron_otherStim__trial_timepoints=cell2mat(spdResid_cNeuron_tempSm__stim__trial_timepoints(setxor(1:nStimuli,[stimIa stimIb])));
                    residuals_cNeuron__trial_timepoints=cat(1,residuals_cNeuron_stimA__trial_timepoints,residuals_cNeuron_stimB__trial_timepoints,residuals_cNeuron_otherStim__trial_timepoints);
                    residualsSSQ_cNeuron__timepoints=sum(residuals_cNeuron__trial_timepoints.^2,1);
                    df=size(residuals_cNeuron__trial_timepoints,1)-sum(stim_LOG);
                    var_cNeuron__timepoints=residualsSSQ_cNeuron__timepoints./df;
                    var__stimPair__fold_timepoints{validStimPairI,1}(foldI,:)=var_cNeuron__timepoints; % >>>
                    % compute ldc
                    ldc__folds_timepoints(foldI,:)=(resp_b_train__timepoints-resp_a_train__timepoints)./var_cNeuron__timepoints.*(resp_b_test__timepoints-resp_a_test__timepoints);
                end % foldI
                RDMs_ldc(stimIa,stimIb,:)=mean(ldc__folds_timepoints,1);
                RDMweights(stimIa,stimIb)=nFolds;
            else
                RDMs_ldc(stimIa,stimIb,:)=nan(1,nTimepoints);
                RDMweights(stimIa,stimIb)=0;
            end
        end % stimPairI      
        
        % mirror ltv to utv
        for timepointI=1:size(RDMs_ldc,3)
            cRDM=RDMs_ldc(:,:,timepointI);
            cRDM_ltv=cRDM(RDM_ltvMask);
            cRDM=squareform(cRDM_ltv);
            RDMs_ldc(:,:,timepointI)=cRDM;
        end        
        RDMweights_ltv=RDMweights(RDM_ltvMask);
        RDMweights=squareform(RDMweights_ltv);
        
        % save and clear
        save(fullfile(resultsPath_singleNeuron,['SPDsingleTrial_4RSA_tempSm_',subjStr{subjectI},'_neuron',num2str(neuronI)]),'spd_cNeuron_tempSm__trial_timepoints');
        save(fullfile(resultsPath_singleNeuron,['RDMs_ldc_',subjStr{subjectI},'_neuron',num2str(neuronI)]),'RDMs_ldc','RDMweights');                
        clear spd_cNeuron_tempSm__trial_timepoints RDMs_ldc RDMweights
    
        % >>>
        save(fullfile(resultsPath_singleNeuron,['var_ldc_',subjStr{subjectI},'_neuron',num2str(neuronI)]),'var__stimPair__fold_timepoints','spdResid_cNeuron_tempSm__stim__trial_timepoints');
        clear var_cNeuron__stimPair__fold_timepoints spdResid_cNeuron_tempSm__stim__trial_timepoints                             
        % >>>
        
    end % neuronI              
end % subjectI



