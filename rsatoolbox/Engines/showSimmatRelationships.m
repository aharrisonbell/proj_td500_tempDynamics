function showSimmatRelationships(simmats,options)
% uses multidimensional scaling to simultaneously relate all simmats passed
% in struct simmats in terms of their similarity (i.e. second-order
% similarity). additionally visualizes the dissimilarity of each simmat to
% the first one (by default) or the one specified as the reference simmat
% in the options. 


%% preparations
simmats=squareSimmats(simmats);
[n,n]=size(squareSimmat(simmats(1).simmat));

if ~exist('options','var'), options.rankTransform=true; end

if ~isfield(options,'rankTransform'), options.rankTransform=false; end
if ~isfield(options,'distanceMeasure'), options.distanceMeasure='euclidean'; end % alternatives: 'correlation','spearman','seuclidean','mahalanobis'
if ~isfield(options,'isomap'), options.isomap=false; end
if ~isfield(options,'nmMDS'), options.nmMDS=true; end
if ~isfield(options,'rubberbands'), options.rubberbands=true; end
if ~isfield(options,'barGraph'), options.barGraph=false; end
if ~isfield(options,'barGraphWithBootStrapErrorBars'), options.barGraphWithBootStrapErrorBars=true; end
if ~isfield(options,'figI'), options.figI=[3600 2 1 1; 3600 2 1 2]; end
if ~isfield(options,'conditionSetIndexVec')||(isfield(options,'conditionSetIndexVec')&&isempty(options.conditionSetIndexVec)), options.conditionSetIndexVec=ones(1,n); end
if ~isfield(options,'criterion'), options.criterion='metricstress'; end
if ~isfield(options,'referenceSimmatI'), options.referenceSimmatI=1; end


options.conditionSetIndexVec=options.conditionSetIndexVec(:)'; % convert to row

distanceMeasureNote='';

nSimmats=size(simmats,2);
figIsI=1; 

%simmatCorrMat(simmats,1);


%% grant index 1 to the reference simmat
simmat_temp=simmats(1);
simmats(1)=simmats(options.referenceSimmatI);
simmats(options.referenceSimmatI)=simmat_temp;


%% reduce simmats to rows and cols defined (non-nan) in all of them
[simmats,validConditionsLOG]=reduceSimmatsToValidConditionSet(simmats);
conditionSetIs_vector=options.conditionSetIndexVec;
conditionSetIs_vector=conditionSetIs_vector(validConditionsLOG);

simmatCorrMat(simmats,1,'spearman');
% simmatCorrMat(simmats,1,'pearson'); % changed by mm on 28 april 2009


%% prepare for bootstrapping with reduced conditions sets
[simmatMask,condSet1_LOG,condSet2_LOG,nCondSets,nCond1,nCond2]=convertToSimmatMask(conditionSetIs_vector);
simmatMask=tril(simmatMask,-1); % retain only upper triangular part (rest set to zero)
%show(simmatMask);


%% optional rank transform and simmat vectorization (simmatMask)
% showSimmats(simmats,2,0);
simmat_rowvecs=nan(nSimmats,sum(simmatMask(:)));
for simmatI=1:nSimmats
    cSimmat_maskContents=simmats(simmatI).simmat(simmatMask);
    if options.rankTransform
        scale01=true;
        cSimmat_maskContents=rankTransform_randomOrderAmongEquals(cSimmat_maskContents,scale01);
    end
    simmats(simmatI).simmat(simmatMask)=cSimmat_maskContents; % copy back into simmats
    simmat_rowvecs(simmatI,:)=cSimmat_maskContents(:);        %  ...and simmat_rowvecs
end
% showSimmats(simmats,2,0);
% showSimmats(simmats,2,1);
%simmatCorrMat(simmats,4);
%show(corrcoef(simmat_rowvecs'));

if options.rankTransform
    rankTransformString='rank-transf. (rand. among eq.)';
else
    rankTransformString='non-rank-transf.';
end


%% compute dissimilarities
try
    D=pdist(simmat_rowvecs,options.distanceMeasure); % compute pairwise distances
catch
    D=pdist(simmat_rowvecs,'euclidean'); % compute pairwise distances
    distanceMeasureNote=[' (',options.distanceMeasure,' dist. FAILED)'];
    options.distanceMeasure='euclidean'
end

% alternatives for options.distanceMeasure: 'correlation', 'euclidean', 'seuclidean', 'mahalanobis', 'spearman', 'cosine', 'cityblock', 'hamming' and others
D=squareSimmat(D);

%show(D);


%% isomap
if options.isomap
    % plots similarity matrices as points in an isomap representation of
    % similarity structure of a set of similarity matrices (2nd-order similarity).
    isomap_options.dims=[1:2];
    isomap_options.display=1;
    [Y, R, E] = isomap_nk(D,'k',6,isomap_options); %[Y, R, E] = isomap(D, n_fcn, n_size, options);
    %[Y, R, E] = isomap(D,'epsilon',2); %[Y, R, E] = isomap(D, n_fcn, n_size, options);

    % plot the isomap arrangement
    coords=Y.coords{2}'; % 2 dimensional coords of isomap arrangement

    selectPlot(options.figI(figIsI,:)); cla; figIsI=figIsI+1;
    if options.rubberbands
        rubberbandGraphPlot(coords,D);
        %         rubberbandGraphPlot(coords,D,'color');
    end
    
    downShift=0;
    for simmatI=1:nSimmats
        plot(coords(simmatI,1),coords(simmatI,2),'o','MarkerFaceColor',simmats(simmatI).color,'MarkerEdgeColor','none','MarkerSize',8); hold on;
        text(coords(simmatI,1),coords(simmatI,2)-downShift,{'',simmats(simmatI).name},'Color',simmats(simmatI).color,'HorizontalAlignment','Center','FontSize',12,'FontName','Arial');
    end
    axis equal tight; zoomOut(10);
    set(gca,'xtick',[],'ytick',[]);
    %xlabel('isomap of simmat space');
    title({'\fontsize{12}isomap\fontsize{9}', ['(',options.distanceMeasure,' dist.',distanceMeasureNote,', ',rankTransformString,' simmats)']});

    shepardPlot(D,[],pdist(coords),986,['isomap(simmats)']);
    
end % isomap

%% nonmetric multidimensional scaling
if options.nmMDS
    % plots similarity matrices as points in an nmMDS representation of
    % similarity structure of a set of similarity matrices (2nd-order similarity).
    nDims=2;
    
    selectPlot(options.figI(figIsI,:)); cla; figIsI=figIsI+1;
    
    try
        [coords, stress, disparities] = mdscale(D, nDims,'criterion',options.criterion);

        % plot the mds arrangement
        if options.rubberbands
            rubberbandGraphPlot(coords,D);
            %         rubberbandGraphPlot(coords,D,'color');
        end

        downShift=0;
        for simmatI=1:nSimmats
            plot(coords(simmatI,1),coords(simmatI,2),'o','MarkerFaceColor',simmats(simmatI).color,'MarkerEdgeColor','none','MarkerSize',8); hold on;
            text(coords(simmatI,1),coords(simmatI,2)-downShift,{'',simmats(simmatI).name},'Color',simmats(simmatI).color,'HorizontalAlignment','Center','FontSize',12,'FontName','Arial');
        end
        axis equal tight; zoomOut(10);
        set(gca,'xtick',[],'ytick',[]);
        %xlabel('nmMDS of simmat space');
        title({['\fontsize{12}\bfdissimilarity-matrix MDS\rm (',options.criterion,')'],['\fontsize{9}(',options.distanceMeasure,' dist.',distanceMeasureNote,', ',rankTransformString,' simmats)']});
        axis off;
        
        % shepardPlot(D,disparities,pdist(coords),987,'non-metric MDS(simmats)');
        shepardPlot(D,disparities,pdist(coords),987,[options.criterion,' MDS(simmats)']); % changed by mm 28 april 2009
    catch
        title('MDS failed.')
    end
end % nmMDS


%% horizontal bar graph of the distances of each simmat from the reference simmat
% (by default the first simmat receives this treatment)
if options.barGraph
    selectPlot(options.figI(figIsI,:)); figIsI=figIsI+1;
    distsToData=D(1,1:end);
    [sortedDistsToData,sortedIndices]=sort(distsToData);

    hb=barh(sortedDistsToData(end:-1:2));
    set(gca,'YTick',[1:nSimmats-1]);
    set(gca,'YTickLabel',{simmats(sortedIndices(end:-1:2)).name},'FontUnits','normalized','FontSize',1/nSimmats);

    if strcmp(options.distanceMeasure,'euclidean')
        xlabel(['euclidean dist.',distanceMeasureNote,' from ',simmats(1).name,' in ',rankTransformString,' simmat space']);
    elseif strcmp(options.distanceMeasure,'spearman')
        xlabel(['(1 - Spearman''s rank corr.) dist.',distanceMeasureNote,' from ',simmats(1).name,' in ',rankTransformString,' simmat space']);
    else
        xlabel([options.distanceMeasure,' dist.',distanceMeasureNote,' from ',simmats(1).name,' in ',rankTransformString,' simmat space']);
    end

    axis tight;
end % bar graph


%% bar graph with bootstrap error bars of the distances of each simmat from the first
if options.barGraphWithBootStrapErrorBars

    % preparations
    if strcmp(options.distanceMeasure,'euclidean')
        titleString=['euclidean dist.',distanceMeasureNote,' from \bf',simmats(1).name,'\rm in ',rankTransformString,' simmat space'];
    elseif strcmp(options.distanceMeasure,'spearman')
        titleString=['(1 - Spearman''s rank corr.) dist.',distanceMeasureNote,' from \bf',simmats(1).name,'\rm in ',rankTransformString,' simmat space'];
    else
        titleString=[options.distanceMeasure,' dist.',distanceMeasureNote,' from \bf',simmats(1).name,'\rm in ',rankTransformString,' simmat space'];
    end
    
    
    % bootstrap resample the similarity matrices
    nBootstrapResamplings=100;

    resampledSimmats_utv=condSetBootstrapOfSimmats_condSetRestriction(simmats,nBootstrapResamplings,conditionSetIs_vector);
    % uses bootstrap resampling of the conditions set to resample a set of simmats.
    % the resampled simmats are returned in upper triangular form (rows), stacked
    % along the 3rd (index of input simmat) and 4th (resampling index)
    % dimensions (for compatibility with square simmats).
    
    % compute second-order similarity matrix for each bootstrap resampling
    bootstrapDistsToRefSimmat=nan(nBootstrapResamplings,nSimmats-1);
    
    for boostrapResamplingI=1:nBootstrapResamplings
         cResampling_simmats_utv=resampledSimmats_utv(:,:,:,boostrapResamplingI);
         cResampling_simmats_utv(:,isnan(mean(cResampling_simmats_utv,3)),:)=[];
         
         %          showSimmats(cResampling_simmats_utv(:,:,[1,5]));
         %          simmatCorrMat(cResampling_simmats_utv(:,:,[1,5]),334);

         % line'em up along the first dimension to please pdist
         simmats_rowvecs=permute(cResampling_simmats_utv,[3 2 1]);
                  
         cResamplingD=squareSimmat(pdist(simmats_rowvecs,options.distanceMeasure)); % compute pairwise distances

         % DEBUG: 
         %          show([cResamplingD,D],1) % TEST PASSED
         %          cResamplingD(1,5)
         %          D(1,5)
         %          pause
         
         bootstrapDistsToRefSimmat(boostrapResamplingI,:)=cResamplingD(1,2:end); % bootstrap distances to the reference simmat (simmat #1)
         % alternatives for options.distanceMeasure: 'correlation', 'euclidean', 'seuclidean', 'mahalanobis', 'spearman', 'cosine', 'cityblock', 'hamming' and others
    end
    
    % sort by the means (bootstrapped: bs)
    [sortedMeanDistsToData_bs,sortedIndices_bs]=sort(mean(bootstrapDistsToRefSimmat,1));
    sortedBootstrapDistsToRefSimmat=bootstrapDistsToRefSimmat(:,sortedIndices_bs);
    
    
    
    % plot bar graph with bootstrap standard-error bars
    col=[0 0 0];
    barNOTplot=true;
    [mean_Y,se_mean_Y]=plotMeanWithStandardDeviationBars(1:nSimmats-1,sortedBootstrapDistsToRefSimmat,options.figI(figIsI,:),titleString,col,barNOTplot)
    %figIsI=figIsI+1;

    % for original simmats (not bootstrapped: nbs)
    distsToData_nbs=D(1,2:end);
    %[sortedDistsToData_nbs,sortedIndices_nbs]=sort(distsToData_nbs);
    plot(1:nSimmats-1,distsToData_nbs(sortedIndices_bs),'.','Color',[.8 .8 .8]);
    %sortedBootstrapDistsToRefSimmat=bootstrapDistsToRefSimmat(:,sortedIndices_nbs);

    
    set(gca,'XTick',[]);
    
    for simmatI=1:nSimmats-1
        col=simmats(sortedIndices_bs(simmatI)+1).color;
        text(simmatI,0,['\bf',simmats(sortedIndices_bs(simmatI)+1).name],'Rotation',90,'Color',col);
    end

    axis_Xmin(0); axis_Xmax(nSimmats);

end % bar graph with bootstrap error bars





