function START_C1_estimateRDMs
% >>> needs final round of debugging

% Version 1: correlation distance

% Author: Marieke Mur; last edit 06-10-2015


%% preparation
clear; close all;

resultsPath='/imaging/mm07/mITReprDynamics/analysis/results';
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));


%% control variables
subjStr={'Stew' 'Wigg'};

RSA_timepoints=-100:700; % ms
RSA_tempSm_window=21; % ms (should be odd)

catOrder_str={'faces','fruit','places','bodyparts','objects'};
save_str='';
changeCatOrder=0;


%% estimate RDMs
for subjectI=1:numel(subjStr)
    % load extracted trial-average responses
    load(fullfile(resultsPath,['SPDexemplar_4RSA_',subjStr{subjectI}]));
    load(fullfile(resultsPath,['diagnostics_',subjStr{subjectI}]),'nStims__neuron');    

    % select neurons that have data for all 100 stimuli
    neuronSelect_allStim_LOG=false(size(nStims__neuron));
    neuronSelect_allStim_LOG(nStims__neuron==100)=true;   
    spdTrialAvg_4RSAcd__neuron_stim_timepoints=spdTrialAvg_4RSA__neuron_stim_timepoints(neuronSelect_allStim_LOG,:,:);   
    nNeurons=sum(neuronSelect_allStim_LOG);
    
    % RDM at each timepoint 
    RSA_timepoints_recod=RSA_timepoints+1000; % spd data ranges from 0 to 5000 ms; 1000 = stimulus onset
    windowWidthOneSide=floor(RSA_tempSm_window/2);    
    
    for timepointI=1:numel(RSA_timepoints)        
        % average trial-average spike-density functions within a temporal window
        windowStart=RSA_timepoints_recod(timepointI)-windowWidthOneSide; windowEnd=RSA_timepoints_recod(timepointI)+windowWidthOneSide;
        spdTrialAvg_4RSAcd_tempSm__neuron_stim_timepoints(:,:,timepointI)=mean(spdTrialAvg_4RSAcd__neuron_stim_timepoints(:,:,windowStart:windowEnd),3);
        
        % compute dissimilarities (correlation distance)
        patterns=[spdTrialAvg_4RSAcd_tempSm__neuron_stim_timepoints(:,:,timepointI)]';
        if changeCatOrder
            catOrderChange_vec=[1:20,61:80,21:40,81:100,41:60]; 
            catOrder_str={'faces','bodyparts','fruit','objects','places'};
            save_str='_catOrderChanged';
            patterns=patterns(catOrderChange_vec,:);
        end      
        RDM_ltv=pdist(patterns,'correlation'); RDM_sq=squareform(RDM_ltv);
        RDMs_cd(:,:,timepointI)=RDM_sq;
    end % timepointI
    
    % save
    save(fullfile(resultsPath,['SPDexemplar_4RSAcd_',subjStr{subjectI}]),'neuronSelect_allStim_LOG','spdTrialAvg_4RSAcd__neuron_stim_timepoints');
    save(fullfile(resultsPath,['SPDexemplar_4RSAcd_tempSm_',subjStr{subjectI}]),'RSA_tempSm_window','spdTrialAvg_4RSAcd_tempSm__neuron_stim_timepoints');
    save(fullfile(resultsPath,['RDMs_cd_',subjStr{subjectI},save_str]),'RDMs_cd','RSA_timepoints','RSA_tempSm_window','catOrder_str','neuronSelect_allStim_LOG','nNeurons');   
    
    clear spdTrialAvg_tempSm__neuron_stim_timepoints RDMs_cd 
end % subjectI



%% show RDMs
load(fullfile(resultsPath,['RDMs_cd_',subjStr{2}]));
cols=RDMcolormap;

timeline=ones(1,numel(RSA_timepoints));
anchorCols=[1 1 1; 0 0 0];
cols_blackWhite=colorScale(anchorCols,100,0);    

for timepointI=1:numel(RSA_timepoints)
    cRDM=RDMs_cd(:,:,timepointI);
    cRDM_rtf=scale01(rankTransform_equalsStayEqual(cRDM,1));
    cRDM_rtf_cmap=(cRDM_rtf*255)+1;
    movieFrames(timepointI)=im2frame(cRDM_rtf_cmap,cols);
    timeline(timepointI)=100;
    movieFrames_timeline(timepointI)=im2frame(timeline,cols_blackWhite);
end

figure(201); clf; h=figure(201);
movie(h,movieFrames,1,25);
figure(301); clf; h=figure(301);
movie(h,movieFrames_timeline,1,25);

fps=500;
save(fullfile(resultsPath,['RDMmovieFrames_cd_',subjStr{subjectI}]),'movieFrames','RSA_timepoints');
movie2avi(movieFrames,fullfile(resultsPath,['RDMmovie_cd_',subjStr{subjectI},'_',num2str(RSA_timepoints(1)),'--',num2str(RSA_timepoints(end)),'ms_',num2str(fps),'fps.avi']),'colormap',cols,'compression','None','fps',fps);








