function rsvp_JN2_fMRI_analysis(prefix);
% function rsvp_JN2_fMRI_analysis(prefix);
% by AHB, February 2011
% Program to load, analyze, and create figures from FMRI data linked to
% RSVP study
% Note: this version is specifically written to address concerns raised in
% the 2nd J Neurosci review 
% PREFIX = Full name of monkey subject - data to load

warning off; close all;
hmiconfig=generate_hmi_configplex;

% PREFIX CORRECTIONS
if strcmp(prefix,'S')==1, prefix='Stewie'; numpatch=4; blocklength=16; end
if strcmp(prefix,'W')==1, prefix='Wiggum'; numpatch=5; blocklength=19; end
if strcmp(prefix,'Stewie')==1,
    sheetname='RSVP Cells_S';
elseif strcmp(prefix,'Wiggum')==1,
    sheetname='RSVP Cells_W';
end

% LOAD DEFAULTS
fmridata_dir=['W:\RSVP_fMRI_Data\',prefix,'_Group',filesep,'Matlab_Analysis',filesep];
figpath = ['\\.psf\Home\Documents\Manuscripts\RSVP_JN_RevisionNo2\Figure_Source_Images\'];
maxroi=6; % maximum number of ROIs for each category
categories={'Faces','Bodyparts','Objects','Places'};
fmri_rsp=struct('responses',zeros(1,4));
bufferlead=2;

% ANALYSIS #1 - Create Average Time Series for each ROI
cd(fmridata_dir)
for cc=1:4, % 1 loop per category 
    figure; clf; cla; set(gcf,'Units','Normalized'); set(gcf,'Position',[0.1 0.1 0.5 0.8]); set(gca,'FontName','Helvetica');
    for rr=1:maxroi, % 1 loop per ROI
        data1=load([prefix,'_AvgTimeSeries_',char(categories(cc)),'3_ROI',num2str(rr),'_norm.1D']);
        data2=load([prefix,'_AvgTimeSeries_',char(categories(cc)),'3_ROI',num2str(rr),'_beta.1D']);
        data3=load([prefix,'_AvgTimeSeries_',char(categories(cc)),'Alt_ROI',num2str(rr),'_norm.1D']);
        data4=load([prefix,'_AvgTimeSeries_',char(categories(cc)),'Alt_ROI',num2str(rr),'_beta.1D']);
        subplot(3,2,rr); hold on
        if isempty(data1)~=1, plot(data1,'k--','LineWidth',2);
            %blocklength=length(data1)/11; % NOTE THIS MAY NEED TO CHANGE FOR WIGGUM
            fmri_rsp(cc).responses(rr,1)=max(data1(17+bufferlead:17+blocklength)); % face response
            fmri_rsp(cc).responses(rr,2)=max(data1(49+bufferlead:49+blocklength)); % bodypart response
            fmri_rsp(cc).responses(rr,3)=max(data1(81+bufferlead:81+blocklength)); % object response
            fmri_rsp(cc).responses(rr,4)=max(data1(113+bufferlead:113+blocklength)); % place response
        end
        if isempty(data2)~=1, plot(data2,'k-','LineWidth',2);
            %blocklength=length(data2)/11; % NOTE THIS MAY NEED TO CHANGE FOR WIGGUM
            fmri_rsp(cc).responses(rr,5)=max(data2(17+bufferlead:17+blocklength)); % face response
            fmri_rsp(cc).responses(rr,6)=max(data2(49+bufferlead:49+blocklength)); % bodypart response
            fmri_rsp(cc).responses(rr,7)=max(data2(81+bufferlead:81+blocklength)); % object response
            fmri_rsp(cc).responses(rr,8)=max(data2(113+bufferlead:113+blocklength)); % place response
        end
        if isempty(data3)~=1, plot(data3,'b--','LineWidth',2);
            %blocklength=length(data3)/11; % NOTE THIS MAY NEED TO CHANGE FOR WIGGUM
            fmri_rsp(cc).responses(rr,9)=max(data3(17+bufferlead:17+blocklength)); % face response
            fmri_rsp(cc).responses(rr,10)=max(data3(49+bufferlead:49+blocklength)); % bodypart response
            fmri_rsp(cc).responses(rr,11)=max(data3(81+bufferlead:81+blocklength)); % object response
            fmri_rsp(cc).responses(rr,12)=max(data3(113+bufferlead:113+blocklength)); % place response
        end
        if isempty(data4)~=1, plot(data4,'b-','LineWidth',2);
            %blocklength=length(data4)/11; % NOTE THIS MAY NEED TO CHANGE FOR WIGGUM
            fmri_rsp(cc).responses(rr,13)=max(data4(17+bufferlead:17+blocklength)); % face response
            fmri_rsp(cc).responses(rr,14)=max(data4(49+bufferlead:49+blocklength)); % bodypart response
            fmri_rsp(cc).responses(rr,15)=max(data4(81+bufferlead:81+blocklength)); % object response
            fmri_rsp(cc).responses(rr,16)=max(data4(113+bufferlead:113+blocklength)); % place response
        end
        set(gca,'XTick',0:blocklength:length(data1))
        for bl=0:blocklength:length(data1), plot([bl bl],[-5 5],'k:','LineWidth',0.5); end
        xlabel('TRs'); ylabel('MION Response (% Change)')
        text(1,2.6,'','FontSize',8,'FontName','Helvetica')
        text(17,2.6,'Fa','FontSize',8,'FontName','Helvetica')
        text(33,2.6,'','FontSize',8,'FontName','Helvetica')
        text(49,2.6,'BP','FontSize',8,'FontName','Helvetica')
        text(65,2.6,'','FontSize',8,'FontName','Helvetica')
        text(81,2.6,'Ob','FontSize',8,'FontName','Helvetica')
        text(97,2.6,'','FontSize',8,'FontName','Helvetica')
        text(113,2.6,'Pl','FontSize',8,'FontName','Helvetica')
        text(129,2.6,'','FontSize',8,'FontName','Helvetica')
        text(145,2.6,'Sc','FontSize',8,'FontName','Helvetica')
        text(161,2.6,'','FontSize',8,'FontName','Helvetica')
        xlim([0 blocklength*9]); ylim([-2 3]); % remove two blocks (scram + baseline)
        title([prefix,' ',char(categories(cc)),' - ROI #',num2str(rr)],'FontSize',10,'FontWeight','Bold','Interpreter','None')
    end
    print(gcf,[figpath,prefix,'',char(categories(cc)),'TimeSeries.ai'],'-dill') % generates an AI file of the figure
    print(gcf,[figpath,prefix,'_',char(categories(cc)),'TimesSeries.jpeg'],'-djpeg') % generates an JPEG file of the figure
end

% ANALYSIS #2 - Bar Graphs showing average fMRI response for each category, for each patch (ROI)
% There are two ways of doing this:
% 1) single value per category (take max/mean from average response profile) 
% 2) average values from each block (take max/mean from each block, i.e., 74ish per monkey)

% 1) Figure based on AVERAGE values (Beta)
figure; clf; cla; set(gcf,'Units','Normalized'); set(gcf,'Position',[0.1 0.1 0.5 0.8]); set(gca,'FontName','Helvetica')
if strcmp('Stewie',prefix)==1, fmricols=5:8; else fmricols=13:16; end
for cc=1:4,
    subplot(4,1,cc)
    bar(1:4,fmri_rsp(cc).responses(1:4,fmricols),'group')
    title([prefix,' ',char(categories(cc))],'FontSize',10,'FontWeight','Bold','Interpreter','None')
    ylim([-2 3])
end
print(gcf,[figpath,prefix,'',char(categories(cc)),'TimeSeriesBarGraphs_BetaAvg.ai'],'-dill') % generates an AI file of the figure
print(gcf,[figpath,prefix,'_',char(categories(cc)),'TimesSeriesBarGraphs_BetaAvg.jpeg'],'-djpeg') % generates an JPEG file of the figure

% 2) Figure based on RAW values (Beta)
fmristruct=r500pn_ProcessfMRIData(prefix,fmridata_dir,blocklength);
figure; clf; cla; set(gcf,'Units','Normalized'); set(gcf,'Position',[0.1 0.1 0.5 0.8]); set(gca,'FontName','Helvetica')
if strcmp(prefix,'Wiggum')==1,
    subplot(4,1,1); hold on
    %bar(1:6,fmristruct.facesALTmean(:,1:4),'group');
    bar(1:6*4,reshape(fmristruct.facesALTmean',1,24))
    errorbar(1:6*4,reshape(fmristruct.facesALTmean',1,24),reshape(fmristruct.facesALTsem',1,24))
    title([prefix,' - Faces']); ylim([-2 5])
    subplot(4,1,2); hold on
    %bar(1:6,fmristruct.placesALTmean(:,1:4),'group');
    bar(1:6*4,reshape(fmristruct.placesALTmean',1,24))
    errorbar(1:6*4,reshape(fmristruct.placesALTmean',1,24),reshape(fmristruct.placesALTsem',1,24))
    title([prefix,' - Places']); ylim([-2 5])
    subplot(4,1,3); hold on
    %bar(1:6,fmristruct.objectsALTmean(:,1:4),'group');
    bar(1:6*4,reshape(fmristruct.objectsALTmean',1,24))
    errorbar(1:6*4,reshape(fmristruct.objectsALTmean',1,24),reshape(fmristruct.objectsALTsem',1,24))
    title([prefix,' - Objects']); ylim([-2 5])
    subplot(4,1,4); hold on
    %bar(1:6,fmristruct.bodypartsALTmean(:,1:4),'group');
    bar(1:6*4,reshape(fmristruct.bodypartsALTmean',1,24))
    errorbar(1:6*4,reshape(fmristruct.bodypartsALTmean',1,24),reshape(fmristruct.bodypartsALTsem',1,24))
    title([prefix,' - Bodyparts']); ylim([-2 5])
else
    subplot(4,1,1); hold on
    %bar(1:6,fmristruct.faces3mean(:,1:4),'group');
    bar(1:6*4,reshape(fmristruct.faces3mean',1,24))
    errorbar(1:6*4,reshape(fmristruct.faces3mean',1,24),reshape(fmristruct.faces3sem',1,24))
    title([prefix,' - Faces']); ylim([-2 5])

    subplot(4,1,2); hold on
    %bar(1:6,fmristruct.places3mean(:,1:4),'group');
    bar(1:6*4,reshape(fmristruct.places3mean',1,24))
    errorbar(1:6*4,reshape(fmristruct.places3mean',1,24),reshape(fmristruct.places3sem',1,24))
    title([prefix,' - Places']); ylim([-2 5])

    subplot(4,1,3); hold on
    %bar(1:6,fmristruct.objects3mean(:,1:4),'group');
    bar(1:6*4,reshape(fmristruct.objects3mean',1,24))
    errorbar(1:6*4,reshape(fmristruct.objects3mean',1,24),reshape(fmristruct.objects3sem',1,24))
    title([prefix,' - Objects']); ylim([-2 5])

    subplot(4,1,4); hold on
    %bar(1:6,fmristruct.bodyparts3mean(:,1:4),'group');
    bar(1:6*4,reshape(fmristruct.bodyparts3mean',1,24))
    errorbar(1:6*4,reshape(fmristruct.bodyparts3mean',1,24),reshape(fmristruct.bodyparts3sem',1,24))
    title([prefix,' - Bodyparts']); ylim([-2 5])

end
print(gcf,[figpath,prefix,'',char(categories(cc)),'TimeSeriesBarGraphs_BetaRaw.ai'],'-dill') % generates an AI file of the figure
print(gcf,[figpath,prefix,'_',char(categories(cc)),'TimesSeriesBarGraphs_BetaRaw.jpeg'],'-djpeg') % generates an JPEG file of the figure


% ANALYSIS #3 - Create Average Spike Density Function for each Patch (ROI)
[grp,grpf,grpbp,grpob,grppl,grpfnf,grpbpnf,grpobnf,grpplnf]=r500pn_GridLocations(prefix);
if exist([hmiconfig.rsvpanal,filesep,prefix(1),'_data.mat'])==2,
    load([hmiconfig.rsvpanal,filesep,prefix(1),'_data.mat']);
else
    [data,numgrids,counts_matrix,allunits,unit_index,unitdata]=plx500_prepproject1data(hmiconfig,sheetname);
end
for pp=1:5, % scroll through each patch
  r500pn_AvgSpdenPatch(grp(pp).grids,unit_index,prefix,figpath,pp);
end

% ANALYSIS #4 - Correlate FMRI response with SPIKING PROPERTIES:
% a) Spiking RawSI vs. FMRI RawSI
% b) Spiking Activity vs. FMRI Response
% c) Distribution vs. 

corrdata=struct('excite_Activity',[],'excite_NormActivity',[],'excite_RawSI',[],'excite_CatSI',[],...
    'sensory_Activity',[],'sensory_NormActivity',[],'sensory_RawSI',[],'sensory_CatSI',[],...
    'fmri_Response',[],'fmri_RawSI',[],'fmri_CatSI',[]); % 1 complete entry / patch

% Extract Data
for np=1:numpatch, % scroll through each patch
    % extract spiking activity
    % Select all SENSORY neurons
    temp_index=find(ismember(unit_index.GridLoc,grp(np).grids)==1 & strcmp(unit_index.SensoryConf,'Sensory')==1);
    corrdata(np).sensory_Activity=unitdata.cat_avg(temp_index,[1 4 5 3]); % faces,bp,objects,places
    corrdata(np).sensory_NormActivity=unitdata.norm_cat_avg(temp_index,[1 4 5 3]); % faces,bp,objects,places
    corrdata(np).sensory_RawSI=unitdata.excite_rawsi_nofruit(temp_index);
    corrdata(np).sensory_CatSI=unitdata.cat_si_nofruit(temp_index,[1 3 4 2]);
    clear temp_index
    corrdata(np).sensory_mean_Activity=mean(corrdata(np).sensory_Activity,1); 
    corrdata(np).sensory_mean_NormActivity=mean(corrdata(np).sensory_NormActivity,1);
    corrdata(np).sensory_mean_RawSI=mean(corrdata(np).sensory_RawSI);
    corrdata(np).sensory_mean_CatSI=mean(corrdata(np).sensory_CatSI,1);
    
    % Select all SENSORY neurons
    temp_index=find(ismember(unit_index.GridLoc,grp(np).grids)==1 & strcmp(unit_index.ExciteConf,'Excite')==1);
    corrdata(np).excite_Activity=unitdata.cat_avg(temp_index,[1 4 5 3]); % faces,bp,objects,places
    corrdata(np).excite_NormActivity=unitdata.norm_cat_avg(temp_index,[1 4 5 3]); % faces,bp,objects,places
    corrdata(np).excite_RawSI=unitdata.excite_rawsi_nofruit(temp_index);
    corrdata(np).excite_CatSI=unitdata.cat_si_nofruit(temp_index,[1 3 4 2]);
    clear temp_index
    corrdata(np).excite_mean_Activity=mean(corrdata(np).excite_Activity,1);
    corrdata(np).excite_mean_NormActivity=mean(corrdata(np).excite_NormActivity,1);
    corrdata(np).excite_mean_RawSI=mean(corrdata(np).excite_RawSI);
    corrdata(np).excite_mean_CatSI=mean(corrdata(np).excite_CatSI,1);
    
    % Calculate Distribution of Neurons per patch
    distrib=[];
    distrib(1)=length(find(ismember(unit_index.GridLoc,grp(np).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
    distrib(2)=length(find(ismember(unit_index.GridLoc,grp(np).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
    distrib(3)=length(find(ismember(unit_index.GridLoc,grp(np).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
    distrib(4)=length(find(ismember(unit_index.GridLoc,grp(np).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));
    corrdata(np).distribution(1)=distrib(1)/sum(distrib)*100;  % 
    corrdata(np).distribution(2)=distrib(2)/sum(distrib)*100;  % 
    corrdata(np).distribution(3)=distrib(3)/sum(distrib)*100;  % 
    corrdata(np).distribution(4)=distrib(4)/sum(distrib)*100;  % 
    clear distrib
    
    % extract FMRI data
    corrdata(np).fmri_Response=r500pn_ExtractfMRIdata(prefix,fmristruct,grp(np).ROI_cat,grp(np).ROI); % faces, bp, objects, places
    corrdata(np).fmri_RawSI=r500pn_calcRawSI(corrdata(np).fmri_Response);
    corrdata(np).fmri_CatSI=r500pn_calcCatSI(corrdata(np).fmri_Response); 
    
end
save([hmiconfig.rsvpanal,filesep,prefix,'_CorrelationData.mat'],'corrdata');


% Generate Correlation Plots
% Compress Data
total_patches=length(corrdata)
corrdata_mini=struct('panelA_x',[],'panelB_x',[],'panelC_x',[],'panelA_y',[],'panelB_y',[],'panelC_y',[]);
for tp=1:total_patches,
    corrdata_mini.panelA_x(tp,1:4)=corrdata(tp).fmri_CatSI;
    corrdata_mini.panelA_y(tp,1:4)=corrdata(tp).excite_mean_CatSI;
    corrdata_mini.panelB_x(tp,1:4)=corrdata(tp).fmri_Response;
    corrdata_mini.panelB_y(tp,1:4)=corrdata(tp).excite_mean_NormActivity;
    corrdata_mini.panelC_x(tp,1:4)=corrdata(tp).fmri_Response;
    corrdata_mini.panelC_y(tp,1:4)=corrdata(tp).distribution;
end

% reshape
panelA_x=reshape(corrdata_mini.panelA_x,4*size(corrdata_mini.panelA_x,1),1);
panelA_y=reshape(corrdata_mini.panelA_y,4*size(corrdata_mini.panelA_y,1),1);
panelB_x=reshape(corrdata_mini.panelB_x,4*size(corrdata_mini.panelB_x,1),1);
panelB_y=reshape(corrdata_mini.panelB_y,4*size(corrdata_mini.panelB_y,1),1);
panelC_x=reshape(corrdata_mini.panelC_x,4*size(corrdata_mini.panelC_x,1),1);
panelC_y=reshape(corrdata_mini.panelC_y,4*size(corrdata_mini.panelC_y,1),1);


figure; clf; cla; set(gcf,'Units','Normalized'); set(gcf,'Position',[0.1 0.1 0.5 0.8]); set(gca,'FontName','Helvetica')
subplot(1,3,1); hold on % compare FMRI SI with Cat SI
plot(panelA_x,panelA_y,'ks','MarkerFaceColor',[0 0 0])
polypar=polyfit(panelA_x,panelA_y,1);
linefit=polyval(polypar,[-0.5 0.5 -0.5 0.5]);
plot([-0.5 0.5 -0.5 0.5],linefit,'k-','LineWidth',1.5);
[rho,pval]=corr(panelA_x,panelA_y);
xlabel('fMRI Selectivity','FontSize',8); ylabel('Neuronal Selectivity','FontSize',8)
xlim([-0.5 0.5]); ylim([-0.5 0.5]); axis square
title({'fMRI vs. Neuronal Selectivity',['rho=',num2str(rho,'%0.2g'),', p=',num2str(pval,'%0.2g')]})

subplot(1,3,2); hold on % compare FMRI activation with Normalized Firing Rate
plot(panelB_x,panelB_y,'ks','MarkerFaceColor',[0 0 0])
polypar=polyfit(panelB_x,panelB_y,1);
linefit=polyval(polypar,[0 3 0.5 1.0]);
plot([0 3 0.5 1.0],linefit,'k-','LineWidth',1.5);
[rho,pval]=corr(panelB_x,panelB_y);
xlabel('fMRI Activation (% Signal Change)','FontSize',8); ylabel('Normalized Neuronal Response (sp/s)','FontSize',8)
xlim([0 3]); ylim([0.5 1.0]); axis square
title({'fMRI vs. Neuronal Response',['rho=',num2str(rho,'%0.2g'),', p=',num2str(pval,'%0.2g')]})

subplot(1,3,3); hold on % compare FMRI activation with Proportion
plot(panelC_x,panelC_y,'ks','MarkerFaceColor',[0 0 0])
polypar=polyfit(panelC_x,panelC_y,1);
linefit=polyval(polypar,[0 3 0 100]);
plot([0 3 0 100],linefit,'k-','LineWidth',1.5);
[rho,pval]=corr(panelC_x,panelC_y);
xlabel('fMRI Activation (% Signal Change)','FontSize',8); ylabel('Neuronal Distribution','FontSize',8)
xlim([0 3]); ylim([0 100]); axis square
title({'fMRI vs. Neuronal Distribution',['rho=',num2str(rho,'%0.2g'),', p=',num2str(pval,'%0.2g')]})


if strcmp(prefix,'B')==1, % group analysis for correlations
   rsvp_JN2_fMRI_analysis_group;
end






return


















%%% Load Spike Data
disp('Loading data for all neurons...')
[data,numgrids,counts_matrix,allunits,unit_index,unitdata]=plx500_prepproject1data(hmiconfig,sheetname);
save([hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'Project1Data_',monkeyname,'.mat'],'data','unit_index','unitdata');






% Figure 2 - Category Selectivity for Each ROI
disp('Figure 2 Category Selectivity for Each ROI')
figure; clf; cla; set(gcf,'Units','Normalized','Position',[0.05 0.15 0.9 0.8]); set(gca,'FontName','Arial','FontSize',8)
%%% MAY EVENTUALLY CHANGE THIS ANALYSIS TO MEASURE RESPONSE FROM EACH
%%% BLOCK, NOT FROM THE AVERAGE TIME SERIES
pref=[1 2 3 4]; categories={'Faces','BodyParts','Objects','Places'}; 
for cc=1:4,
    subplot(2,2,cc)
    numrois=size(fmri_rsp(cc).responses,1); bardata=zeros(numrois,1);
    for rr=1:numrois,
        Rp=fmri_rsp(cc).responses(rr,cc);
        Rnp=sum(fmri_rsp(cc).responses(rr,pref~=cc))/3;
        bardata(rr)=(Rp-Rnp)/(Rp+Rnp);
    end

    bar([bardata ; mean(bardata)]);
    xlabel('ROI #'); ylabel(['Average ',char(categories(cc)),' SI'])
    title([char(categories(cc)),' ROIs']); ylim([0 0.75]);
    text(length(bardata)+1,0.7,'AVERAGE','HorizontalAlignment','Center');
end
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVPfMRI_Fig2_ROITimeSeries_',monkeyname,'_',char(categories(cc)),'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVPfMRI_Fig2_ROITimeSeries_',monkeyname,'_',char(categories(cc)),'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)


% Figure 3 - Average Time Series for Each Patch within recording grid
disp('Figure 3 Average Time Series for Each Patch within recording grid')
figure; clf; cla; set(gcf,'Units','Normalized','Position',[0.05 0.15 0.9 0.8]); set(gca,'FontName','Arial','FontSize',8)
fmri_faces=struct('responses',[],'average',[]);
fmri_bodyparts=struct('responses',[],'average',[]);
fmri_objects=struct('responses',[],'average',[]);
fmri_places=struct('responses',[],'average',[]);

for pp=1:5,
    subplot(4,5,pp); hold on; % Faces IN NEAR OUT
    temp_pointer=find(ismember(unit_index.GridLoc,grpfnf(pp).grids)==1 & strcmp(unit_index.SensoryConf,'Sensory')==1);
    temp_coords=unique(unit_index.APcoords(temp_pointer,:),'Rows');
    temp_data=[];
    for tc=1:size(temp_coords,1),
        temp_afni_coords(tc,1:3)=rsvp_apml2xyz(monkinitial,temp_coords(tc,:));
    end
    temp_fmri_locs=unique(temp_afni_coords,'Rows');
    for tf=1:size(temp_fmri_locs,1),
        newdata=squeeze(fmridata(temp_fmri_locs(tf,1),temp_fmri_locs(tf,2),temp_fmri_locs(tf,3),:))';
        temp_data=[temp_data;newdata];
        fmri_faces(pp).responses(tf,1)=max(newdata(17+bufferlead:17+blocklength)); % face response
        fmri_faces(pp).responses(tf,2)=max(newdata(49+bufferlead:49+blocklength)); % bodypart response
        fmri_faces(pp).responses(tf,3)=max(newdata(81+bufferlead:81+blocklength)); % object response
        fmri_faces(pp).responses(tf,4)=max(newdata(113+bufferlead:113+blocklength)); % place response
    end
    newdata=mean(temp_data,1);
    plot(newdata);
    fmri_faces(pp).average(1)=max(newdata(17+bufferlead:17+blocklength)); % face response
    fmri_faces(pp).average(2)=max(newdata(49+bufferlead:49+blocklength)); % bodypart response
    fmri_faces(pp).average(3)=max(newdata(81+bufferlead:81+blocklength)); % object response
    fmri_faces(pp).average(4)=max(newdata(113+bufferlead:113+blocklength)); % place response
    plot([16 16],[-2 3],'k:','LineWidth',0.5); plot([32 32],[-2 3],'k:','LineWidth',0.5);
    plot([48 48],[-2 3],'k:','LineWidth',0.5); plot([64 64],[-2 3],'k:','LineWidth',0.5);
    plot([80 80],[-2 3],'k:','LineWidth',0.5); plot([96 96],[-2 3],'k:','LineWidth',0.5);
    plot([112 112],[-2 3],'k:','LineWidth',0.5); plot([128 128],[-2 3],'k:','LineWidth',0.5);
    xlim([0 145]); ylim([-2 3]); ylabel('% Signal Change'); xlabel('#TRs');
    title(['Faces - Patch #',num2str(pp),'(',num2str(size(temp_data,1)),')'],'FontSize',8,'FontWeight','Bold')
    clear temp_fmri_locs temp_afni_coords temp_pointer temp_coords newdata
end

for pp=1:3,
    subplot(4,3,pp+3); hold on; % Bodyparts IN NEAR OUT
    temp_pointer=find(ismember(unit_index.GridLoc,grpbpnf(pp).grids)==1 & strcmp(unit_index.SensoryConf,'Sensory')==1);
    temp_coords=unique(unit_index.APcoords(temp_pointer,:),'Rows');
    temp_data=[];
    for tc=1:size(temp_coords,1),
        temp_afni_coords(tc,1:3)=rsvp_apml2xyz(monkinitial,temp_coords(tc,:));
    end
    temp_fmri_locs=unique(temp_afni_coords,'Rows');
    for tf=1:size(temp_fmri_locs,1),
        newdata=squeeze(fmridata(temp_fmri_locs(tf,1),temp_fmri_locs(tf,2),temp_fmri_locs(tf,3),:))';
        temp_data=[temp_data;newdata];
        fmri_bodyparts(pp).responses(tf,1)=max(newdata(17+bufferlead:17+blocklength)); % face response
        fmri_bodyparts(pp).responses(tf,2)=max(newdata(49+bufferlead:49+blocklength)); % bodypart response
        fmri_bodyparts(pp).responses(tf,3)=max(newdata(81+bufferlead:81+blocklength)); % object response
        fmri_bodyparts(pp).responses(tf,4)=max(newdata(113+bufferlead:113+blocklength)); % place response
    end
    newdata=mean(temp_data,1);
    plot(newdata);
    fmri_bodyparts(pp).average(1)=max(newdata(17+bufferlead:17+blocklength)); % face response
    fmri_bodyparts(pp).average(2)=max(newdata(49+bufferlead:49+blocklength)); % bodypart response
    fmri_bodyparts(pp).average(3)=max(newdata(81+bufferlead:81+blocklength)); % object response
    fmri_bodyparts(pp).average(4)=max(newdata(113+bufferlead:113+blocklength)); % place response
    plot([16 16],[-2 3],'k:','LineWidth',0.5); plot([32 32],[-2 3],'k:','LineWidth',0.5);
    plot([48 48],[-2 3],'k:','LineWidth',0.5); plot([64 64],[-2 3],'k:','LineWidth',0.5);
    plot([80 80],[-2 3],'k:','LineWidth',0.5); plot([96 96],[-2 3],'k:','LineWidth',0.5);
    plot([112 112],[-2 3],'k:','LineWidth',0.5); plot([128 128],[-2 3],'k:','LineWidth',0.5);
    xlim([0 145]); ylim([-2 3]); ylabel('% Signal Change'); xlabel('#TRs');
    title(['BodyParts - Patch #',num2str(pp),'(',num2str(size(temp_data,1)),')'],'FontSize',8,'FontWeight','Bold')
    clear temp_fmri_locs temp_afni_coords temp_pointer temp_coords newdata
end

for pp=1:3,
    subplot(4,3,pp+6); hold on; % Objects IN NEAR OUT
    temp_pointer=find(ismember(unit_index.GridLoc,grpobnf(pp).grids)==1 & strcmp(unit_index.SensoryConf,'Sensory')==1);
    temp_coords=unique(unit_index.APcoords(temp_pointer,:),'Rows');
    temp_data=[];
    for tc=1:size(temp_coords,1),
        temp_afni_coords(tc,1:3)=rsvp_apml2xyz(monkinitial,temp_coords(tc,:));
    end
    temp_fmri_locs=unique(temp_afni_coords,'Rows');
    for tf=1:size(temp_fmri_locs,1),
        newdata=squeeze(fmridata(temp_fmri_locs(tf,1),temp_fmri_locs(tf,2),temp_fmri_locs(tf,3),:))';
        temp_data=[temp_data;newdata];
        fmri_objects(pp).responses(tf,1)=max(newdata(17+bufferlead:17+blocklength)); % face response
        fmri_objects(pp).responses(tf,2)=max(newdata(49+bufferlead:49+blocklength)); % bodypart response
        fmri_objects(pp).responses(tf,3)=max(newdata(81+bufferlead:81+blocklength)); % object response
        fmri_objects(pp).responses(tf,4)=max(newdata(113+bufferlead:113+blocklength)); % place response
    end
    newdata=mean(temp_data,1);
    plot(newdata);
    fmri_objects(pp).average(1)=max(newdata(17+bufferlead:17+blocklength)); % face response
    fmri_objects(pp).average(2)=max(newdata(49+bufferlead:49+blocklength)); % bodypart response
    fmri_objects(pp).average(3)=max(newdata(81+bufferlead:81+blocklength)); % object response
    fmri_objects(pp).average(4)=max(newdata(113+bufferlead:113+blocklength)); % place response
    plot([16 16],[-2 3],'k:','LineWidth',0.5); plot([32 32],[-2 3],'k:','LineWidth',0.5);
    plot([48 48],[-2 3],'k:','LineWidth',0.5); plot([64 64],[-2 3],'k:','LineWidth',0.5);
    plot([80 80],[-2 3],'k:','LineWidth',0.5); plot([96 96],[-2 3],'k:','LineWidth',0.5);
    plot([112 112],[-2 3],'k:','LineWidth',0.5); plot([128 128],[-2 3],'k:','LineWidth',0.5);
    xlim([0 145]); ylim([-2 3]); ylabel('% Signal Change'); xlabel('#TRs');
    title(['Objects - Patch #',num2str(pp),'(',num2str(size(temp_data,1)),')'],'FontSize',8,'FontWeight','Bold')
    clear temp_fmri_locs temp_afni_coords temp_pointer temp_coords newdata
end

for pp=1:3,
    subplot(4,3,pp+9); hold on; % Places IN NEAR OUT
    temp_pointer=find(ismember(unit_index.GridLoc,grpplnf(pp).grids)==1 & strcmp(unit_index.SensoryConf,'Sensory')==1);
    temp_coords=unique(unit_index.APcoords(temp_pointer,:),'Rows');
    temp_data=[];
    for tc=1:size(temp_coords,1),
        temp_afni_coords(tc,1:3)=rsvp_apml2xyz(monkinitial,temp_coords(tc,:));
    end
    temp_fmri_locs=unique(temp_afni_coords,'Rows');
    for tf=1:size(temp_fmri_locs,1),
        newdata=squeeze(fmridata(temp_fmri_locs(tf,1),temp_fmri_locs(tf,2),temp_fmri_locs(tf,3),:))';
        temp_data=[temp_data;newdata];
        fmri_places(pp).responses(tf,1)=max(newdata(17+bufferlead:17+blocklength)); % face response
        fmri_places(pp).responses(tf,2)=max(newdata(49+bufferlead:49+blocklength)); % bodypart response
        fmri_places(pp).responses(tf,3)=max(newdata(81+bufferlead:81+blocklength)); % object response
        fmri_places(pp).responses(tf,4)=max(newdata(113+bufferlead:113+blocklength)); % place response
    end
    newdata=mean(temp_data,1);
    plot(newdata);
    fmri_places(pp).average(1)=max(newdata(17+bufferlead:17+blocklength)); % face response
    fmri_places(pp).average(2)=max(newdata(49+bufferlead:49+blocklength)); % bodypart response
    fmri_places(pp).average(3)=max(newdata(81+bufferlead:81+blocklength)); % object response
    fmri_places(pp).average(4)=max(newdata(113+bufferlead:113+blocklength)); % place response
    plot([16 16],[-2 3],'k:','LineWidth',0.5); plot([32 32],[-2 3],'k:','LineWidth',0.5);
    plot([48 48],[-2 3],'k:','LineWidth',0.5); plot([64 64],[-2 3],'k:','LineWidth',0.5);
    plot([80 80],[-2 3],'k:','LineWidth',0.5); plot([96 96],[-2 3],'k:','LineWidth',0.5);
    plot([112 112],[-2 3],'k:','LineWidth',0.5); plot([128 128],[-2 3],'k:','LineWidth',0.5);
    xlim([0 145]); ylim([-2 3]); ylabel('% Signal Change'); xlabel('#TRs');
    title(['Places - Patch #',num2str(pp),'(',num2str(size(temp_data,1)),')'],'FontSize',8,'FontWeight','Bold')
    clear temp_fmri_locs temp_afni_coords temp_pointer temp_coords newdata
end
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVPfMRI_Fig3_PatchTimeSeries_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVPfMRI_Fig3_PatchTimeSeries_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)



% Figure 4 - Category Selectivity for Each Patch within Grid (and Correlation)
disp('Figure 4 Category Selectivity for Each ROI')
figure; clf; cla; set(gcf,'Units','Normalized','Position',[0.05 0.15 0.9 0.8]); set(gca,'FontName','Arial','FontSize',8)
%%% MAY EVENTUALLY CHANGE THIS ANALYSIS TO MEASURE RESPONSE FROM EACH
%%% BLOCK, NOT FROM THE AVERAGE TIME SERIES
pref=[1 2 3 4]; categories={'Faces','BodyParts','Objects','Places'};
clear bardata
subplot(2,2,1); % faces
for rr=1:5,
    Rp=fmri_faces(rr).average(1);
    Rnp=sum(fmri_faces(rr).average(pref~=1))/3;
    bardata(rr)=(Rp-Rnp)/(Rp+Rnp);
end
bar(bardata);
xlabel('Location'); ylabel('Average Faces SI')
title(['Face Patches']); ylim([-0.5 0.5]); set(gca,'XTickLabel',{'AntIn','AntNear','PostIn','PostNear','Far'})
clear bardata

subplot(2,2,2); % bodyparts
for rr=1:3,
    Rp=fmri_bodyparts(rr).average(2);
    Rnp=sum(fmri_bodyparts(rr).average(pref~=2))/3;
    bardata(rr)=(Rp-Rnp)/(Rp+Rnp);
end
bar(bardata);
xlabel('Location'); ylabel('Average BP SI')
title(['BodyPart Patches']); ylim([-0.5 0.5]); set(gca,'XTickLabel',{'In','Near','Far'})
clear bardata

subplot(2,2,3); % objects
for rr=1:3,
    Rp=fmri_objects(rr).average(3);
    Rnp=sum(fmri_objects(rr).average(pref~=3))/3;
    bardata(rr)=(Rp-Rnp)/(Rp+Rnp);
end
bar(bardata);
xlabel('Location'); ylabel('Average Objects SI')
title(['Object Patches']); ylim([-0.5 0.5]); set(gca,'XTickLabel',{'In','Near','Far'})
clear bardata

subplot(2,2,4); % places
for rr=1:3,
    Rp=fmri_places(rr).average(4);
    Rnp=sum(fmri_places(rr).average(pref~=4))/3;
    bardata(rr)=(Rp-Rnp)/(Rp+Rnp);
end
bar(bardata);
xlabel('Location'); ylabel('Average Places SI')
title(['Place Patches']); ylim([-.5 0.5]); set(gca,'XTickLabel',{'In','Near','Far'})
clear bardata
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVPfMRI_Fig4_GridSI_',monkeyname,'_',char(categories(cc)),'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVPfMRI_Fig4_GridSI_',monkeyname,'_',char(categories(cc)),'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)





return

%%%% NESTED FUNCTIONS %%%%


function AFNIcoords=rsvp_apml2xyz(monkinitial,APcoords);
% Converts grid coords (AP,ML) to hard-coded X,Y,Z slice coords for
% comparing neuronal data to EPI timeseries data.
% Note: this program uses hard-coded values that were determined by
% comparing electrode scans and EPI/anatomicals
% monkinitial = 'S' or 'W'
% ap = anterior/posterior grid location (5-19mm)
% ml = medial/lateral grid location (13-27mm)

ap=APcoords(1); ml=APcoords(2);
if monkinitial=='S' % Stewie % Updated June 15, 2010
    APrange=5:19;
    MLrange=13:27;
    sag_range=[0 0 0 17 17 16 15 14 13 13 12 11 10 9 0]; % y slice# (Medial to Lateral)
    axi_range=[29 29 28 28 26 26 26 25 25 25 24 24 23 23 23]; % z slice# (Posterior to Anterior)
    cor_range=[19 19 18 18 17 17 16 16 15 15 14 14 13 13 12]; % x slice# (Posterior to Anterior)
    AFNIcoords(1)=sag_range(find(MLrange==ml)); % sagital
    AFNIcoords(2)=axi_range(find(APrange==ap)); % axial
    AFNIcoords(3)=cor_range(find(APrange==ap)); % coronal
    
    %% Correct for Matrix Orientation
    AFNIcoords(1)=64-AFNIcoords(1);
    AFNIcoords(2)=64-AFNIcoords(2);
    AFNIcoords(3)=(20+5)-AFNIcoords(3);
elseif monkinitial=='W', % wiggum
    APrange=5:19;
    MLrange=14:29;
    cor_range=[20 19 18 17 16 15 14 13 12 11 10 9 8 7 7]; % x (must match to AFNI output)
    sag_range=[38 39 40 41 42 43 44 45 46 47 48 49 50 51 52]; % y
    axi_range=[28 28 28 28 28 26 26 21 23 22 23 23 22 22 20]; % z
    
    man_range=[22 20 20 24 0 0 0;22 23 21 22 0 0 0;23 23 22 22 22 0 0;23 23 23 23 23 0 0;23 23 23 23 23 23 0;23 22 22 22 22 24 0;23 21 21 21 23 23 23;24 21 21 21 21 21 21;25 26 26 26 26 26 26;...
        25 26 26 26 26 26 26;26 28 28 28 28 29 29;27 28 28 28 28 28 28;28 28 28 28 28 28 28;28 28 28 28 28 28 28;29 28 28 28 28 28 28];
    xrange=[7 7 8 9 10 11 12 13 14 15 16 17 18 19 20];
    yrange=[43 44 45 46 47 48 49];
    
    AFNIcoords(3)=cor_range(APrange==ap); % coronal
    AFNIcoords(1)=sag_range(MLrange==ml); % sagital
    AFNIcoords(2)=axi_range(APrange==ap); % axial
    
    
    
    %% Note: Wiggum's data were clipped Ant/Post - 5 from Ant, 10, from
    %% post.  Therefore:
    AFNIcoords(3)=AFNIcoords(3)-5;
end
return
     