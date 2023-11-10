function START_D4_visualiseRDMs_selectNeurons(responseProp,valueIs)

% FUNCTION
% Select subset of neurons and then compute the RDMs 
%
% INPUT
% responseProp      according to which response property do we want to select neurons?
%                   options: 'visResp' 'catSel' 'stimSel' 'excInh' 'quality'
% valueIs           which value(s) on the response property do we want to select?
%                   can be a scalar or vector of indices
%                   visResp: 0 = no, 1 = yes
%                            loaded data already only include visually-responsive neurons
%                   catSel:  1 = faces, 2 = fruit, 3 = places, 4 = bodyparts, 5 = objects    
%                   stimSel: 0 = no, 1 = yes
%                   excInh:  1 = exc, 2 = inh, 3 = both
%                   quality: 0-3; 3 = good?
%
% AUTHOR 
% Marieke Mur; last edit: 15-12-2016 


%% preparation
close all;

resultsPath='/imaging/mm07/mITReprDynamics/analysis/results';
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));

responseVals(1).strings={'visResp' 'notVisResp'}; 
responseVals(2).strings={'F' 'Fr' 'P' 'B' 'O'};
responseVals(3).strings={'stimSel' 'notStimSel'};
responseVals(4).strings={'E' 'I' 'EI'};
responseVals(5).strings={'Q0' 'Q1' 'Q2' 'Q3'};    


%% control variables
subjStr={'Stew' 'Wigg'};
monitor=1;


%% load single-neuron RDMs and select neurons
for subjectI=1:numel(subjStr)
    PSfilespec=fullfile(resultsPath,['RDMdiagnostics_',subjStr{subjectI},'.ps']);
    load(fullfile(resultsPath,['RDMs_srd_',subjStr{subjectI}]),'RDMweights_neuron','RDMs_neuron_timepoints');
    load(fullfile(resultsPath,['respProp_4RSA_',subjStr{subjectI}]),'responseProp_4RSA__neuron_prop','responseProp_str');
    load(fullfile(resultsPath,'START_C4a_variables'),'RSA_timepoints');
    
    % select response property
    propI=find(strcmp(responseProp_str,responseProp));
    selectedProp=responseProp_4RSA__neuron_prop(:,propI);
    
    % select neurons
    neuronIs=[];
    for i=1:numel(valueIs)
        valueI=valueIs(i);
        neuronIs=[neuronIs;find(selectedProp==valueI)];
    end
    neuronIs=sort(neuronIs);
    RDMs_neuron_timepoints=RDMs_neuron_timepoints(:,:,neuronIs,:);
    RDMweights_neuron=RDMweights_neuron(:,:,neuronIs);
    
    % compute weighted-average spike rate distance over the selected neurons
    RDMs_neuron_timepoints_weighted=RDMs_neuron_timepoints.*repmat(RDMweights_neuron,[1 1 1 size(RDMs_neuron_timepoints,4)]);
    RDMs_neuron_timepoints_weighted_summed=squeeze(nansum(RDMs_neuron_timepoints_weighted,3));
    RDMs_timepoints=RDMs_neuron_timepoints_weighted_summed./repmat(sum(RDMweights_neuron,3),[1 1 size(RDMs_neuron_timepoints,4)]);
    save(fullfile(resultsPath,['RDMs_srd_weighted_',responseVals(propI).strings{valueIs},'_',subjStr{subjectI}]),'RDMs_timepoints');
    
    if monitor
        figI=200; figure(figI); clf;
        hist(RDMs_timepoints(:));
        xlabel({'spike rate distance','\fontsize{8}[weighted-average, cross-validated, squared(?) Eucl distance in repr space spanned by the neurons]'});
        ylabel('frequency');
        title('Distribution of spike rate distances in final population RDM');
        description(1)={'\fontsize{11}population SRDs'};
        description(2)={['\fontsize{10}',responseVals(propI).strings{valueIs}]};
        description(3)={['\fontsize{9}',num2str(numel(neuronIs)),' neurons']};
        addHeadingAndPrint(description,PSfilespec,figI);
    end
    
    % create an RDM movie
    cols=RDMcolormap;
    
    % ranktransform dissims and scale into [0 1]
    for timepointI=1:numel(RSA_timepoints)
        cRDM=RDMs_timepoints(:,:,timepointI);
        cRDM_rtf=scale01(rankTransform_equalsStayEqual(cRDM,1));
        cRDM_rtf_cmap=(cRDM_rtf*255)+1;
        movieFrames(timepointI)=im2frame(cRDM_rtf_cmap,cols);
    end
    if monitor 
        figure(201); clf; h=figure(201);
        movie(h,movieFrames,1,25);
    end
    fps=25;
    save(fullfile(resultsPath,['RDMmovieFrames_srd_',responseVals(propI).strings{valueIs},'_',subjStr{subjectI}]),'movieFrames','RSA_timepoints');
    movie2avi(movieFrames,fullfile(resultsPath,['RDMmovie_srd_',responseVals(propI).strings{valueIs},'_',subjStr{subjectI},'_',num2str(RSA_timepoints(1)),'--',num2str(RSA_timepoints(end)),'ms_',num2str(fps),'fps.avi']),'colormap',cols,'compression','None','fps',fps);
    
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
        movieFrames(timepointI)=im2frame(cRDM_sc_cmap,cols);        
    end
    if monitor
        figure(202); clf; h=figure(202);
        movie(h,movieFrames,1,25);
    end
    fps=25;
    save(fullfile(resultsPath,['RDMmovieFrames_srd_origScale_',responseVals(propI).strings{valueIs},'_',subjStr{subjectI}]),'movieFrames','RSA_timepoints','extreme');
    movie2avi(movieFrames,fullfile(resultsPath,['RDMmovie_srd_origScale_',responseVals(propI).strings{valueIs},'_',subjStr{subjectI},'_',num2str(RSA_timepoints(1)),'--',num2str(RSA_timepoints(end)),'ms_',num2str(fps),'fps.avi']),'colormap',cols,'compression','None','fps',fps);
    
end % subjectI