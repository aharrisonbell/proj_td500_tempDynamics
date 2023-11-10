function START_D3_visualiseRDMs

% Dissimilarity measure = single-neuron linear-discriminant contrast (START_C3a)

% Author: Marieke Mur; last edit 10-10-2015


%% preparation
clear; close all;

resultsPath='/imaging/mm07/mITReprDynamics/analysis/results';
resultsPath_singleNeuron=fullfile(resultsPath,'singleNeuron');
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));


%% control variables
subjStr={'Stew' 'Wigg'};
nStimuli=100;

monitor=1;


%% visualise RDMs
for subjectI=1:numel(subjStr)
    load(fullfile(resultsPath,['respProp_4RSA_',subjStr{subjectI}]),'responseProp_4RSA__neuron_prop','responseProp_str');
    load(fullfile(resultsPath,['dataSelectionInfo_4RSA_',subjStr{subjectI}]),'neuronSelect_LOG','nTrials__neuron_stim','nValidTrials__neuron_stim');
    load(fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),'nStims__neuron','nSingleTrialStims__neuron');
    load(fullfile(resultsPath,'START_C3a_variables'),'RSA_timepoints');
    nNeurons=sum(neuronSelect_LOG);
    nTimepoints=numel(RSA_timepoints);
    
    % pre-allocate
    RDMs_neuron_timepoints=zeros(nStimuli,nStimuli,nNeurons,nTimepoints);
    RDMweights_neuron=zeros(nStimuli,nStimuli,nNeurons);
    if monitor, minVariance__neuron=nan(nNeurons,1); maxVariance__neuron=minVariance__neuron; end
    
    for neuronI=1:nNeurons
        % load RDMs
        load(fullfile(resultsPath_singleNeuron,['RDMs_ldc_',subjStr{subjectI},'_neuron',num2str(neuronI)]));
        
        % perform sanity check on LDCs and variance
        if monitor
            description(1)={['\fontsize{12}neuron ',num2str(neuronI)]};
            PSfilespec=fullfile(resultsPath,['RDMdiagnostics_',subjStr{subjectI},'.ps']);
            
            % store minimum and maximum variance for current neuron
            load(fullfile(resultsPath_singleNeuron,['var_ldc_',subjStr{subjectI},'_neuron',num2str(neuronI)]),'var__stimPair__fold_timepoints');
            var_cNeuron__stimPair__fold_timepoints_mat=cell2mat(var__stimPair__fold_timepoints);
            minVariance__neuron(neuronI)=min(var_cNeuron__stimPair__fold_timepoints_mat(:));
            maxVariance__neuron(neuronI)=max(var_cNeuron__stimPair__fold_timepoints_mat(:));
            clear var__stimPair__fold_timepoints var_cNeuron__stimPair__fold_timepoints_mat % to free up memory            
            
            % LDC histogram for current neuron
            nan_LOG=isnan(RDMs_ldc); % dissimilarities are NaN if either both or one of the contributing stimuli have <2 trials (= not enough data)
            nonnanIs=find(nan_LOG==0);
            RDMs_ldc_nonnans=RDMs_ldc(nonnanIs);
            
            figI=100; pageFigure(figI); clf;
            subplot(2,1,1); hist(RDMs_ldc_nonnans(:));
            xlabel('ldc values'); ylabel('frequency');
            title({'Histogram of LDCs across stimulus pairs and timepoints','\fontsize{8}[NaNs excluded]'});
            
            minLDC__neuron(neuronI)=min(RDMs_ldc_nonnans);
            maxLDC__neuron(neuronI)=max(RDMs_ldc_nonnans);
            
            % check that NaNs in RDMs indeed only occur if there was not enough data (indicated by a zero weight in RDMweights)
            diagMask=logical(eye(nStimuli,nStimuli));
            RDMweights(diagMask)=1;
            nanIs_inferredFromWeights=find(RDMweights==0);
            for timepointI=1:size(RDMs_ldc,3)
                nan_LOG_cTimepoint=isnan(RDMs_ldc(:,:,timepointI));
                nanIs=find(nan_LOG_cTimepoint);
                nNans(neuronI,timepointI)=numel(nanIs);
                if numel(nanIs_inferredFromWeights)~=numel(nanIs)
                    disp(['neuron ',num2str(neuronI),', timepoint ',num2str(timepointI),': nans present in nonzero weight entries']);
                else
                    nanOverlap(neuronI,timepointI)=numel(find(nanIs_inferredFromWeights(:)==nanIs(:)));
                    if nanOverlap~=numel(nanIs)
                        disp(['neuron ',num2str(neuronI),', timepoint ',num2str(timepointI),': nan locations do not match nonzero weight locations']);
                    end
                    nanMatch(neuronI,timepointI)=100;
                end
            end % timepointI
            figure(figI);
            subplot(2,1,2);
            plot(nNans(neuronI,:)); hold on;
            plot(nanOverlap(neuronI,:),'r'); hold on;
            legend({'nr NaNs in RDM','nr zeros in RDM weight matrix'});
            xlabel('timepoints'); ylabel('nr of missing dissimilarities (NaNs or zeros)');
            title({'Number of missing dissimilarities at each timepoint','\fontsize{8}[dissimilarities are missing if one or both of the contributing stimuli have <2 trials]'});
            
            % print figure
            addHeadingAndPrint(description,PSfilespec,figI);
        end % monitor
        
        % assemble RDMs across neurons
        RDMs_neuron_timepoints(:,:,neuronI,:)=RDMs_ldc;
        RDMweights_neuron(:,:,neuronI)=RDMweights;
        clear RDMs_ldc RDMweights
    end % neuronI
    save(fullfile(resultsPath,['RDMs_ldc_',subjStr{subjectI}]),'RDMweights_neuron','RDMs_neuron_timepoints','-v7.3');
    if monitor, save(fullfile(resultsPath,['RDMdiagnostics_',subjStr{subjectI}]),'minVariance__neuron','maxVariance__neuron','minLDC__neuron','maxLDC__neuron'); end
%     load(fullfile(resultsPath,['RDMs_ldc_',subjStr{subjectI}]),'RDMweights_neuron','RDMs_neuron_timepoints');
%     load(fullfile(resultsPath,['RDMdiagnostics_',subjStr{subjectI}]),'minVariance__neuron','maxVariance__neuron','minLDC__neuron','maxLDC__neuron');
    
    if monitor
        % plot single-neuron variance (min and max)
        figI=250; pageFigure(figI); clf;
        plot(minVariance__neuron); hold on;
        plot(maxVariance__neuron,'Color',[.8 .8 1]); hold on;
        varSmaller1_neuronIs=find(minVariance__neuron<1);
        for lowVarNeuronI=1:numel(varSmaller1_neuronIs)
            plot(varSmaller1_neuronIs(lowVarNeuronI),0,'o','MarkerEdgeColor','r','MarkerFaceColor','r','MarkerSize',4);
        end
        xlabel('neurons'); ylabel('variance');
        legend({'min variance','max variance','variance < 1'});
        title('Minimum and maximum stimulus-response variance');
        description(1)={'\fontsize{11}single-neuron response variance [min max]'};
        addHeadingAndPrint(description,PSfilespec,figI);
        
        % plot single-neuron LDCs (min and max)
        figI=255; pageFigure(figI); clf;
        plot(minLDC__neuron,'g'); hold on;
        plot(maxLDC__neuron,'Color',[.8 1 .8]); hold on;
        for lowVarNeuronI=1:numel(varSmaller1_neuronIs)
            plot(varSmaller1_neuronIs(lowVarNeuronI),0,'o','MarkerEdgeColor',[.8 .8 .8],'MarkerFaceColor',[.8 .8 .8],'MarkerSize',4);
        end
        xlabel('neurons'); ylabel('LDC value');
        title('Minimum and maximum LDC');
        legend({'min LDC','max LDC','variance < 1'});
        description(1)={'\fontsize{11}single-neuron LDCs [min max]'};
        addHeadingAndPrint(description,PSfilespec,figI);
        
        figI=260; pageFigure(figI); clf;
        LDCsmaller1000_neuronIs=find(minLDC__neuron<-1e3);
        LDCbigger1000_neuronIs=find(maxLDC__neuron>1e3);
        removeIs=unique([LDCsmaller1000_neuronIs,LDCbigger1000_neuronIs]);
        varSmaller1_neuronIs_woExtrLDCs=setxor(varSmaller1_neuronIs,removeIs);
        plot(minLDC__neuron(setxor(1:nNeurons,removeIs)),'g'); hold on;
        plot(maxLDC__neuron(setxor(1:nNeurons,removeIs)),'Color',[.8 1 .8]); hold on;
        for lowVarNeuronI=1:numel(varSmaller1_neuronIs_woExtrLDCs)
            plot(varSmaller1_neuronIs_woExtrLDCs(lowVarNeuronI),0,'o','MarkerEdgeColor',[.8 .8 .8],'MarkerFaceColor',[.8 .8 .8],'MarkerSize',4);
        end
        xlabel('neurons'); ylabel('LDC');
        title({'Minimum and maximum LDC (extreme neurons removed)','\fontsize{8}[extreme = LDC < -1000 or LDC > 1000]',['[',num2str(numel(removeIs)),' neurons removed]']});
        legend({'min LDC','max LDC','variance < 1'});
        description(1)={'\fontsize{11}single-neuron LDCs [min max]'};
        addHeadingAndPrint(description,PSfilespec,figI);
    end % monitor
    
    %     excInh=responseProp_4RSA__neuron_prop(:,4);
    %     inhNeuronIs=find(excInh==2);
    %     LDCextremeInh_LOG=ismember(removeIs,inhNeuronIs);
    %
    %     lessThan9foldIs=find(RDMweights_neuron<9);
    %     for timepointI=1:nTimepoints
    %         RDMs_neuron_cTimepoint=RDMs_neuron_timepoints(:,:,:,timepointI);
    %         RDMs_neuron_cTimepoint(lessThan9foldIs)=0;
    %         RDMs_neuron_timepoints(:,:,:,timepointI)=RDMs_neuron_cTimepoint;
    %     end
    
    % compute weighted-average over neurons
    nan_LOG=isnan(RDMs_neuron_timepoints); nanIs=find(nan_LOG); nNans=numel(nanIs);
    RDMs_neuron_timepoints_weighted=RDMs_neuron_timepoints.*repmat(RDMweights_neuron,[1 1 1 size(RDMs_neuron_timepoints,4)]);
    RDMs_neuron_timepoints_weighted_summed=squeeze(nansum(RDMs_neuron_timepoints_weighted,3));
    RDMs_timepoints=RDMs_neuron_timepoints_weighted_summed./repmat(sum(RDMweights_neuron,3),[1 1 size(RDMs_neuron_timepoints,4)]);
    save(fullfile(resultsPath,['RDMs_ldc_weighted_',subjStr{subjectI}]),'RDMs_timepoints');
    
    % create an RDM movie
    cols=RDMcolormap;
    for timepointI=1:numel(RSA_timepoints)
        cRDM=RDMs_timepoints(:,:,timepointI);
        cRDM_rtf=scale01(rankTransform_equalsStayEqual(cRDM,1));
        cRDM_rtf_cmap=(cRDM_rtf*255)+1;
        movieFrames(timepointI)=im2frame(cRDM_rtf_cmap,cols);
    end
    figure(201); clf; h=figure(201);
    movie(h,movieFrames,1,25);
    fps=25;
    save(fullfile(resultsPath,['RDMmovieFrames_ldc_',subjStr{subjectI}]),'movieFrames','RSA_timepoints');
    movie2avi(movieFrames,fullfile(resultsPath,['RDMmovie_ldc_',subjStr{subjectI},'_',num2str(RSA_timepoints(1)),'--',num2str(RSA_timepoints(end)),'ms_',num2str(fps),'fps.avi']),'colormap',cols,'compression','None','fps',fps);
    
end % subjectI