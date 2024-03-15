function START_D5_MDSplots

% Dissimilarity measure = single-neuron spike rate distance (START_C4a, START_D4)

% Author: Marieke Mur; last edit: 04-12-2017 


%% preparation
clear; close all;

resultsPath = '/imaging/mm07/mITReprDynamics/analysis/results';
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox'));


%% control variables
subjStr = {'Stew' 'Wigg'};
nStimuli = 100;


%% plot MDS solutions
for subjectI = 1:numel(subjStr)
    load(fullfile(resultsPath,['RDMs_srd_weighted_',subjStr{subjectI}]),'RDMs_timepoints');
    pats_mds_2D_tminus1 = [];
    
    for timepointI = 1:size(RDMs_timepoints,3)
        cRDM = RDMs_timepoints(:,:,timepointI);
        
        % predefined categories
        categoryLabels = {'faces','fruits','places','bodyparts','objects'};
        categoryColours = [.9 .15 .15; 1 0.65 0; .3 .3 .9; .65 0 .65; .3 .8 .3];
        c_faces = [ones(20,1); zeros(80,1)];
        c_fruits = [zeros(20,1); ones(20,1); zeros(60,1)];
        c_places = [zeros(40,1); ones(20,1); zeros(40,1)];
        c_bodyparts = [zeros(60,1); ones(20,1); zeros(20,1)];
        c_objects = [zeros(80,1); ones(20,1)];
        C = [c_faces, c_fruits, c_places, c_bodyparts, c_objects];
        
        % set MDS control variables
        criterion = 'metricstress';
        nDims = 2;
        
        % preparation
        figI_MDS = 325; figure(figI_MDS); clf(figI_MDS);
        figI_MDS_shepard = 350; figure(figI_MDS_shepard); clf(figI_MDS_shepard);
        figI_MDS_procrustes = 425; figure(figI_MDS_procrustes); clf(figI_MDS_procrustes);
        figI_MDS_shepard_procrustes = 450; figure(figI_MDS_shepard_procrustes); clf(figI_MDS_shepard_procrustes);
        
        % MDS
        cRDM_neg = cRDM < 0;
        nNegSRDs = sum(cRDM_neg(:))/2;
        cRDM(cRDM_neg) = 0;
        [pats_mds_2D,stress,disparities] = mdscale(cRDM,nDims,'criterion',criterion);
        
        % plot category-coloured dots
        figure(figI_MDS);
        plotDots(pats_mds_2D,C,categoryColours);
        ys(timepointI,:) = ylim; xs(timepointI,:) = xlim;
        ylim([-5 5]); xlim([-5 5]);
        text(3,4,['time = ',num2str(timepointI - 101),'ms']);
        axis square off;
        movie_mds(timepointI) = getframe;
        
        % shepard-plots
        shepardPlot(cRDM,disparities,pdist(pats_mds_2D),figI_MDS_shepard,['\fontsize{12}shepard plot\fontsize{10}','(',criterion,')']);
        
        % procrustes alignment to MDS solution of previous timepoint
        if isempty(pats_mds_2D_tminus1)
            [pats_mds_2D_procrustes,mu,sigma] = zscore(pats_mds_2D,1);
            [nPats,nCoords] = size(pats_mds_2D_procrustes);
            pats_mds_2D_procrustes = pats_mds_2D_procrustes.*repmat(sigma,[nPats 1]) + repmat(mu,[nPats 1]);
        else
            [pats_mds_2D_tminus1,mu_tminus1,sigma_tminus1] = zscore(pats_mds_2D_tminus1,1);
            [pats_mds_2D,mu,sigma] = zscore(pats_mds_2D,1);
            [d,pats_mds_2D_procrustes] = procrustes(pats_mds_2D_tminus1,pats_mds_2D);
            [nPats,nCoords] = size(pats_mds_2D_procrustes);
            pats_mds_2D_procrustes = pats_mds_2D_procrustes.*repmat(sigma,[nPats 1]) + repmat(mu,[nPats 1]);
        end
        
        % plot category-coloured dots
        figure(figI_MDS_procrustes);
        plotDots(pats_mds_2D_procrustes,C,categoryColours);
        ys_procrustes(timepointI,:) = ylim; xs_procrustes(timepointI,:) = xlim;
        ylim([-5 5]); xlim([-5 5]);
        text(3,4,['time = ',num2str(timepointI - 101),'ms']);
        axis square off;
        movie_mds_procrustes(timepointI) = getframe;
        
        % shepard-plots
        shepardPlot(cRDM,disparities,pdist(pats_mds_2D_procrustes),figI_MDS_shepard_procrustes,['\fontsize{12}shepard plot\fontsize{10}','(',criterion,')']);
    
        pats_mds_2D_tminus1 = pats_mds_2D_procrustes;
    end % timepointI
    
    save(fullfile(resultsPath,'MDSmovie'),'movie_mds');
    save(fullfile(resultsPath,'MDSmovie_procrustes'),'movie_mds_procrustes');
    
    v1 = VideoWriter('MDSmovie.avi');
    open(v1);
    writeVideo(v1,movie_mds);
        
    v2 = VideoWriter('MDSmovie_procrustes.avi');
    open(v2);
    writeVideo(v2,movie_mds_procrustes);
    
    v3 = VideoWriter('MDSmovie_procrustes_highRes.avi','Uncompressed AVI');
    open(v3);
    writeVideo(v3,movie_mds_procrustes);
end % subjectI



