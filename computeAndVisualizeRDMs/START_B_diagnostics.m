function START_B_diagnostics

% DESCRIPTION   
% This script computes diagnostics for visually-responsive neurons. Best to
% run the code for each subject separately (we are clearing variables 
% throughout the script but might have missed some).

% AUTHOR 
% Marieke Mur; last edit: 05-10-2015 


%% preparation
clear; close all; 

resultsPath='/imaging/mm07/mITReprDynamics/analysis/results';
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));


%% control variables
subjStr={'Stew' 'Wigg'};
nStimuli=100;

plotSingleNeuronResp=0;


%% diagnostics
for subjectI=1:numel(subjStr)
    % load extracted patterns and trial information
    load(fullfile(resultsPath,['dataSelectionInfo_4RSA_',subjStr{subjectI}]),'nTrials__neuron_stim','nValidTrials__neuron_stim');
    load(fullfile(resultsPath,['SPDexemplar_4RSA_',subjStr{subjectI}]),'spdTrialAvg_4RSA__neuron_stim_timepoints');
    load(fullfile(resultsPath,['trialInfo_4RSA_',subjStr{subjectI}]),'trialInfo_4RSA__neuron__trial_stimCat');
    load(fullfile(resultsPath,['SPDsingleTrial_4RSA_',subjStr{subjectI}]),'spd_4RSA__neuron__trial_timepoints');
    
    % (1) inspect distribution of valid trials over neurons (for each stimulus)
    % how many neurons have 2 or more trials for each stimulus?
    min2ValidTrials_LOG__neuron_stim=false(size(nValidTrials__neuron_stim));
    min2ValidTrials_LOG__neuron_stim(nValidTrials__neuron_stim>1)=1;
    min2ValidTrials_LOG_summedAcrossStim__neuron=sum(min2ValidTrials_LOG__neuron_stim,2);
    cvAmenableNeuronIs=find(min2ValidTrials_LOG_summedAcrossStim__neuron==100);
    % show number of valid trials over neurons and stimuli
    figI=199; clf;
    figure(figI);
    subplot(3,1,1); imagesc(nValidTrials__neuron_stim); colormap(gca, RDMcolormap); colorbar;
    xlabel('stimuli'); ylabel('neurons'); title('number of valid trials (no threshold)')
    subplot(3,1,2); imagesc(min2ValidTrials_LOG__neuron_stim); colormap(gca, RDMcolormap); colorbar;
    xlabel('stimuli'); ylabel('neurons'); title({'number of valid trials (thresholded)','\fontsize{8}[0 indicates < 2, 1 indicates >= 2]'});
    subplot(3,1,3); plot(histc(min2ValidTrials_LOG_summedAcrossStim__neuron,1:nStimuli));
    xlabel('number of stimuli for which there are >= 2 trials'); ylabel('number of neurons'); title('do most neurons have >= 2 trials for all 100 stimuli?');
    description(1)={'\fontsize{12}number of valid trials'};
    description(2)={['\fontsize{10}visually-responsive neurons (n=',num2str(size(nValidTrials__neuron_stim,1)),')']};
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
    description(3)={'\fontsize{10}all stimuli'};
    addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),figI);
    % save
    save(fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),'min2ValidTrials_LOG__neuron_stim','min2ValidTrials_LOG_summedAcrossStim__neuron','cvAmenableNeuronIs',...
        'tables_stim','minnValidTrials_stim','maxnValidTrials_stim','meannValidTrials_stim');
    clear tables_stim minnValidTrials_stim maxnValidTrials_stim meannValidTrials_stim description
    close all
    
    % (2) plot single-neuron responses (for each stimulus)
    if plotSingleNeuronResp
        nNeurons=size(spdTrialAvg_4RSA__neuron_stim_timepoints,1);
        anchorCols=[1 1 1; .7 .7 .7; 0 0 0];
        zeroNanColormap=colorScale(anchorCols,100,0);
        description(1)={'\fontsize{12}trial-average responses'};
        
        for neuronI=1:nNeurons
            description(2)={['\fontsize{10}visually-responsive neuron ',num2str(neuronI)]};
            spdTrialAvg_cNeuron__stim_timepoints=squeeze(spdTrialAvg_4RSA__neuron_stim_timepoints(neuronI,:,:));
            minimum=min(spdTrialAvg_cNeuron__stim_timepoints(:)); % minimum is zero for all neurons and subjects
             
            % show single-neuron responses (trial-average for each stimulus)
            figI=301; figure(figI); clf;
            imagesc(spdTrialAvg_cNeuron__stim_timepoints); colormap(gca,RDMcolormap); colorbar;
            xlabel({'timepoints (ms)','\fontsize{8}[1000 = stim onset]'}); ylabel('stimuli');
            title({'trial-average spike-density functions',['\fontsize{8}minimum = ',num2str(minimum)]});
            axis square
            addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_exemplarResp_',subjStr{subjectI}]),figI);
            
            % find zeros and nans            
            zeroIs=find(spdTrialAvg_cNeuron__stim_timepoints==0);
            nan_LOG=isnan(spdTrialAvg_cNeuron__stim_timepoints); nanIs=find(nan_LOG);
            zerosAndNans=0.5*ones(size(spdTrialAvg_cNeuron__stim_timepoints));
            zerosAndNans(zeroIs)=0;
            zerosAndNans(nanIs)=1;
            % show zeros and nans
            figI=302; figure(figI); clf;
            imagesc(zerosAndNans,[0 1]); colormap(gca,zeroNanColormap); colorbar;
            xlabel({'timepoints (ms)','\fontsize{8}[1000 = stim onset]'}); ylabel('stimuli');
            title('gray = nonzero data, white = zeros, black = NaNs (no data)');                        
            axis square
            addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_exemplarResp_',subjStr{subjectI}]),figI);                            
        end % neuronI            
        clear description; close all;
    end
    
    % (3) plot single-neuron response variance
    nNeurons=numel(spd_4RSA__neuron__trial_timepoints);
    for neuronI=1:nNeurons
        trialInfo_cNeuron__trial_stimCat=trialInfo_4RSA__neuron__trial_stimCat{neuronI};
        spd_cNeuron__trial_timepoints=spd_4RSA__neuron__trial_timepoints{neuronI};
        stim_LOG=true(nStimuli,1);        
        stim_singleTrial_LOG=false(nStimuli,1);        
        
        for stimulusI=1:nStimuli
            trialIs=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimulusI);
            if isempty(trialIs), stim_LOG(stimulusI)=false; 
            elseif numel(trialIs)==1, stim_singleTrial_LOG(stimulusI)=true;
            end
            % extract single-trial and trial-avg spd functions for the current stimulus
            spd_cNeuron_cStim__trial_timepoints=spd_cNeuron__trial_timepoints(trialIs,:);
            spd_cNeuron_cStim_trialAvg__timepoints=mean(spd_cNeuron_cStim__trial_timepoints,1);
            % compute residuals
            spdResid_cNeuron_cStim__trial_timepoints=spd_cNeuron_cStim__trial_timepoints-repmat(spd_cNeuron_cStim_trialAvg__timepoints,[numel(trialIs) 1]);
            % store to pass on
            spdResid_cNeuron__stim__trial_timepoints{stimulusI,1}=spdResid_cNeuron_cStim__trial_timepoints;
        end % stimulusI
        
        % compute response variance for current neuron
        residuals_cNeuron=cell2mat(spdResid_cNeuron__stim__trial_timepoints);
        residualsSSQ_cNeuron=sum(residuals_cNeuron.^2,1);
        df=size(residuals_cNeuron,1)-sum(stim_LOG);
        var_cNeuron=residualsSSQ_cNeuron./df;
        
        % assemble for saving
        var__neuron_timepoints(neuronI,:)=var_cNeuron; % variance should be nonnegative (checked and is indeed the case for both subjects)
        nStims__neuron(neuronI,1)=sum(stim_LOG);
        nSingleTrialStims__neuron(neuronI,1)=sum(stim_singleTrial_LOG);
        
        clear spdResid_cNeuron__stim__trial_timepoints        
    end % neuronI
    save(fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),'var__neuron_timepoints','nStims__neuron','nSingleTrialStims__neuron','-append');    
    
    % plot single-neuron "stimulus coverage"
    figI=401; figure(figI); clf; 
    plot(nStims__neuron); hold on;
    plot(nSingleTrialStims__neuron,'r');
    xlabel('neurons'); ylabel('stimuli');    
    legend({'nr of stimuli that neuron has data for','nr of stimuli that neuron has 1 trial for'});
    description(1)={'\fontsize{12}single-neuron "stimulus coverage"'};
    description(2)={['\fontsize{10}visually-responsive neurons (n=',num2str(nNeurons),')']};
    addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),figI);
    
    % show single-neuron response variance
    figI=501; figure(figI); clf;
    imagesc(var__neuron_timepoints); colormap(gca, RDMcolormap); colorbar;
    xlabel({'timepoints (ms)','\fontsize{8}[1000 = stim onset]'}); ylabel('neurons');
    title({'response variance',...
        '\fontsize{8}[residuals = single-trial MINUS trial-avg spd functions]',...
        '\fontsize{8}[df = nr of trials MINUS nr of stimuli]'});
    axis square
    description(1)={'\fontsize{12}single-neuron response variance'};
    description(2)={['\fontsize{10}visually-responsive neurons (n=',num2str(nNeurons),')']};
    addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),figI);    
    % find zeros (no NaNs present)
    zeroIs=find(var__neuron_timepoints==0);
    nan_LOG=isnan(var__neuron_timepoints); nanIs=find(nan_LOG); % >>>
    var_zeros=ones(size(var__neuron_timepoints));
    var_zeros(zeroIs)=0;    
    % show zeros 
    figI=502; figure(figI); clf;
    anchorCols=[1 1 1; .7 .7 .7];
    zeroColormap=colorScale(anchorCols,100,0);
    imagesc(var_zeros,[0 1]); colormap(gca,zeroColormap); colorbar;
    xlabel({'timepoints (ms)','\fontsize{8}[1000 = stim onset]'}); ylabel('neurons');
    title('gray = nonzero variance, white = zero variance');
    axis square
    addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),figI);
    
    clear var__neuron_timepoints nStims__neuron nSingleTrialStims__neuron description trialIs    
    close all
    
    % (4) how many neurons have identical trial information? (and were thus likely measured simultaneously with the same electrode)
    nNeurons=numel(trialInfo_4RSA__neuron__trial_stimCat);
    identicalTrialInfo_neuron2neuronMat=zeros(nNeurons,nNeurons);
    neuron2neuronMat_mask=logical(triu(ones(nNeurons,nNeurons)));
    for neuronI=1:nNeurons, nTrials(neuronI)=size(trialInfo_4RSA__neuron__trial_stimCat{neuronI},1); end
    nTrialsUnique=unique(nTrials);
    for uniqueI=1:numel(nTrialsUnique)
        neuronIs_sameNTrials=find(nTrials==nTrialsUnique(uniqueI));
        if numel(neuronIs_sameNTrials) > 1
            neuron2neuronMat=zeros(nNeurons,nNeurons);
            neuron2neuronMat(neuronIs_sameNTrials,neuronIs_sameNTrials)=1;
            neuron2neuronMat(neuron2neuronMat_mask)=0;
            pairIs=find(neuron2neuronMat==1);
            [nIsa,nIsb]=ind2sub(size(neuron2neuronMat),pairIs);
            for pairIsI=1:numel(pairIs)
                trialIs(:,1)=trialInfo_4RSA__neuron__trial_stimCat{nIsa(pairIsI)}(:,1);
                trialIs(:,2)=trialInfo_4RSA__neuron__trial_stimCat{nIsb(pairIsI)}(:,1);
                overlap=size(find(trialIs(:,1)==trialIs(:,2)),1);
                if overlap==size(trialIs,1)
                    identicalTrialInfo_neuron2neuronMat(nIsa(pairIsI),nIsb(pairIsI))=uniqueI;
                end
                clear trialIs
            end % pairIsI
        end
    end % uniqueI
    identicalTrialInfo_neuron2neuronMat_LOG=false(nNeurons,nNeurons);
    identicalTrialInfo_neuron2neuronMat_LOG(identicalTrialInfo_neuron2neuronMat>0)=true;
    save(fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),'identicalTrialInfo_neuron2neuronMat','identicalTrialInfo_neuron2neuronMat_LOG','-append');
    
    % show
    anchorCols=[.4 .4 .4; 1 1 1];
    grayWhiteColormap=colorScale(anchorCols,100,0);    
    figI=601; figure(figI); clf; imagesc(identicalTrialInfo_neuron2neuronMat_LOG); colormap(gca,grayWhiteColormap); colorbar;
    xlabel('neuron'); ylabel('neuron'); title({'\fontsize{8}[1 = identical trial sequence]'});
    axis square;
    description(1)={['\fontsize{12}neurons with identical trial sequence (',num2str(sum(identicalTrialInfo_neuron2neuronMat_LOG(:))),' pairs)']};
    description(2)={['\fontsize{10}visually-responsive neurons (n=',num2str(nNeurons),')']};
    addHeadingAndPrint(description,fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),figI);
    clear description; close all;
    
end % subjectI






