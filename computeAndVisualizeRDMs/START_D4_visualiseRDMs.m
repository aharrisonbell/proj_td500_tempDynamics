function START_D4_visualiseRDMs

% Dissimilarity measure = single-neuron spike rate distance (START_C4a)

% Author: Marieke Mur; last edit 04-12-2017


%% preparation
clear; close all;

resultsPath='/imaging/mm07/mITReprDynamics/analysis/results';
resultsPath_singleNeuron=fullfile(resultsPath,'singleNeuron');
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));


%% control variables
subjStr={'Stew' 'Wigg'};
nStimuli=100;

monitor=0;


%% visualise RDMs
for subjectI=1:numel(subjStr)
    load(fullfile(resultsPath,['respProp_4RSA_',subjStr{subjectI}]),'responseProp_4RSA__neuron_prop','responseProp_str');
    load(fullfile(resultsPath,['dataSelectionInfo_4RSA_',subjStr{subjectI}]),'neuronSelect_LOG','nTrials__neuron_stim','nValidTrials__neuron_stim');
    load(fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),'nStims__neuron','nSingleTrialStims__neuron');
    load(fullfile(resultsPath,'START_C4a_variables'),'RSA_timepoints');
    nNeurons=sum(neuronSelect_LOG);
    nTimepoints=numel(RSA_timepoints);
    
    % pre-allocate
    RDMs_neuron_timepoints=zeros(nStimuli,nStimuli,nNeurons,nTimepoints);
    RDMweights_neuron=zeros(nStimuli,nStimuli,nNeurons);   
    if monitor, minSRD__neuron=nan(nNeurons,1); maxSRD__neuron=minSRD__neuron; end

    for neuronI=1:nNeurons
        % load RDMs
        load(fullfile(resultsPath_singleNeuron,['RDMs_srd_',subjStr{subjectI},'_neuron',num2str(neuronI)]));
        
        % perform sanity check on SRDs and variance
        if monitor
            description(1)={['\fontsize{12}neuron ',num2str(neuronI)]};
            PSfilespec=fullfile(resultsPath,['RDMdiagnostics_',subjStr{subjectI},'.ps']);
                        
            % SRD histogram for current neuron
            nan_LOG=isnan(RDMs_srd); % dissimilarities are NaN if either both or one of the contributing stimuli have <2 trials (= not enough data)
            nonnanIs=find(nan_LOG==0);
            RDMs_srd_nonnans=RDMs_srd(nonnanIs);
            
            figI=100; pageFigure(figI); clf;
            subplot(2,1,1); hist(RDMs_srd_nonnans(:));
            xlabel('srd values'); ylabel('frequency');
            title({'Histogram of SRDs across stimulus pairs and timepoints','\fontsize{8}[NaNs excluded]'});
            
            minSRD__neuron(neuronI)=min(RDMs_srd_nonnans);
            maxSRD__neuron(neuronI)=max(RDMs_srd_nonnans);
            
            % check that NaNs in RDMs indeed only occur if there was not enough data (indicated by a zero weight in RDMweights)
            diagMask=logical(eye(nStimuli,nStimuli));
            RDMweights(diagMask)=1;
            nanIs_inferredFromWeights=find(RDMweights==0);
            for timepointI=1:size(RDMs_srd,3)
                nan_LOG_cTimepoint=isnan(RDMs_srd(:,:,timepointI));
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
            for timepointI=1:size(RDMs_srd,3)
                nan_LOG_cTimepoint=isnan(RDMs_srd(:,:,timepointI));
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
            % print figure
            addHeadingAndPrint(description,PSfilespec,figI);
        end % monitor
        
        % assemble RDMs across neurons
        RDMs_neuron_timepoints(:,:,neuronI,:)=RDMs_srd;
        RDMweights_neuron(:,:,neuronI)=RDMweights;
        clear RDMs_srd RDMweights
    end % neuronI
    save(fullfile(resultsPath,['RDMs_srd_',subjStr{subjectI}]),'RDMweights_neuron','RDMs_neuron_timepoints','-v7.3');
    if monitor, save(fullfile(resultsPath,['RDMdiagnostics_',subjStr{subjectI}]),'minSRD__neuron','maxSRD__neuron'); end
    load(fullfile(resultsPath,['RDMs_srd_',subjStr{subjectI}]),'RDMweights_neuron','RDMs_neuron_timepoints');
    load(fullfile(resultsPath,['RDMdiagnostics_',subjStr{subjectI}]),'minVariance__neuron','maxVariance__neuron','minSRD__neuron','maxSRD__neuron');
    
    if monitor      
        % plot single-neuron SRDs (min and max)
        figI=255; pageFigure(figI); clf;
        plot(minSRD__neuron,'g'); hold on;
        plot(maxSRD__neuron,'Color',[.8 1 .8]); hold on;
        xlabel('neurons'); ylabel('SRD value');
        title('Minimum and maximum SRD');
        legend({'min SRD','max SRD'});
        description(1)={'\fontsize{11}single-neuron SRDs [min max]'};
        addHeadingAndPrint(description,PSfilespec,figI);
        
        figI=260; pageFigure(figI); clf;
        SRDsmaller1000_neuronIs=find(minSRD__neuron<-1e3);
        SRDbigger1000_neuronIs=find(maxSRD__neuron>1e3);
        removeIs=unique([SRDsmaller1000_neuronIs,SRDbigger1000_neuronIs]);        
        plot(minSRD__neuron(setxor(1:nNeurons,removeIs)),'g'); hold on;
        plot(maxSRD__neuron(setxor(1:nNeurons,removeIs)),'Color',[.8 1 .8]); hold on;
        xlabel('neurons'); ylabel('SRD');
        title({'Minimum and maximum SRD (extreme neurons removed)','\fontsize{8}[extreme = SRD < -1000 or SRD > 1000]',['[',num2str(numel(removeIs)),' neurons removed]']});
        legend({'min SRD','max SRD'});
        description(1)={'\fontsize{11}single-neuron SRDs [min max]'};
        addHeadingAndPrint(description,PSfilespec,figI);
    end % monitor   
    
    % compute weighted-average spike rate distance over neurons
    nan_LOG=isnan(RDMs_neuron_timepoints); nanIs=find(nan_LOG); nNans=numel(nanIs);
    RDMs_neuron_timepoints_weighted=RDMs_neuron_timepoints.*repmat(RDMweights_neuron,[1 1 1 size(RDMs_neuron_timepoints,4)]);
    RDMs_neuron_timepoints_weighted_summed=squeeze(nansum(RDMs_neuron_timepoints_weighted,3));
    RDMs_timepoints=RDMs_neuron_timepoints_weighted_summed./repmat(sum(RDMweights_neuron,3),[1 1 size(RDMs_neuron_timepoints,4)]);
    save(fullfile(resultsPath,['RDMs_srd_weighted_',subjStr{subjectI}]),'RDMs_timepoints');
    
    figI=200; clf; pageFigure(figI);
    hist(RDMs_timepoints(:)); 
    xlabel({'spike rate distance','\fontsize{8}[weighted-average, cross-validated, squared(?) Eucl distance in repr space spanned by the neurons]'});
    ylabel('frequency');
    title('Distribution of spike rate distances in final population RDM');
    description(1)={'\fontsize{11}population SRDs'};
    description(2)={['\fontsize{10}',num2str(nNeurons),' neurons']};
    addHeadingAndPrint(description,PSfilespec,figI);
    
    % create an RDM movie
    cols=RDMcolormap;
    
    % ranktransform dissims and scale into [0 1]
    for timepointI=1:numel(RSA_timepoints)
        cRDM=RDMs_timepoints(:,:,timepointI);
        cRDM_rtf=scale01(rankTransform_equalsStayEqual(cRDM,1));
        cRDM_rtf_cmap=(cRDM_rtf*255)+1;
        movieFrames(timepointI)=im2frame(cRDM_rtf_cmap,cols);
    end
    figure(201); clf; h=figure(201);
    movie(h,movieFrames,1,25);
    fps=25;
    save(fullfile(resultsPath,['RDMmovieFrames_srd_',subjStr{subjectI}]),'movieFrames','RSA_timepoints');
    movie2avi(movieFrames,fullfile(resultsPath,['RDMmovie_srd_',subjStr{subjectI},'_',num2str(RSA_timepoints(1)),'--',num2str(RSA_timepoints(end)),'ms_',num2str(fps),'fps.avi']),'colormap',cols,'compression','None','fps',fps);
    
    % show original spike rate distances    
    minimum=min(RDMs_timepoints(:));
    maximum=max(RDMs_timepoints(:));
    extreme=max([abs(minimum) abs(maximum)]);
    range=2*extreme;
    colBinSize=range/size(cols,1);
    for timepointI=1:numel(RSA_timepoints)
        cRDM=RDMs_timepoints(:,:,timepointI);
        cRDM_sc_cmap=(ceil((cRDM+extreme)./colBinSize)); % rescale: neg extreme corresponds now to zero, then determine for each SRD which colour bin it is in        
        cRDM_sc_cmap(cRDM_sc_cmap==0)=1; % zeros are in the first colour bin
        showRDMs(cRDM,500,0,[-1*abs(extreme) abs(extreme)]);
        showRDMs(cRDM_sc_cmap,501,0,[1 256]);
        movieFrames(timepointI)=im2frame(cRDM_sc_cmap,cols);        
    end
    figure(202); clf; h=figure(202);
    movie(h,movieFrames,1,25);
    fps=25;
    save(fullfile(resultsPath,['RDMmovieFrames_srd_origScale_',subjStr{subjectI}]),'movieFrames','RSA_timepoints','extreme');
    movie2avi(movieFrames,fullfile(resultsPath,['RDMmovie_srd_origScale_',subjStr{subjectI},'_',num2str(RSA_timepoints(1)),'--',num2str(RSA_timepoints(end)),'ms_',num2str(fps),'fps.avi']),'colormap',cols,'compression','None','fps',fps);
    
end % subjectI