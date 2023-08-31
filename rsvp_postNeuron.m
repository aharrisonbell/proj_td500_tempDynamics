function rsvp_postNeuron(monkinitial);
%%%%%%%%%%%%%%%%%%%
% rsvp_postNeuron %
%%%%%%%%%%%%%%%%%%%
% written by AHB, Nov2009,
% This program includes old and new analyses for an individual subject
% since the neuron submission.
% MONKINITIAL (required) = 'W' or 'S'

%%% SETUP DEFAULTS
warning off; close all;
hmiconfig=generate_hmi_configplex; % generates and loads config file
if nargin==0, error('You must specify a monkey (''S''/''W'')'); elseif nargin==1, neurtype='both'; end

%%%% Grids defined according to fMRI category-selective regions
% if monkinitial=='S',
%     monkeyname='Stewie'; sheetname='RSVP Cells_S';
%     % Grid location groups for comparison
%     grp(1).grids={'A7L2','A7L1','A7R1'}; % BodyPart Selective
%     grp(2).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Face Selective
%     grp(3).grids={'A4L2','A4L1','A4R1'}; % Place Selective
%     grp(4).grids={'A2L5','A0L6','A0L2','A0L0','A0R0','P1L1','P2L3','P3L5','P3L4','P5L2','P5L3','P6L3'}; % Object Selective
%     grp(5).grids={'P6L2','P6L1','P7L2'}; % Face Selective
% elseif monkinitial=='W',
%     monkeyname='Wiggum'; sheetname='RSVP Cells_W';
%     % Grid location groups for comparison
%     grp(1).grids={'A6R2','A5R0','A4R3'}; % Bodypart Selective
%     grp(2).grids={'R0','A2R1','A2R3','A2R5'}; % Face Selective
%     grp(3).grids={'P1R0','P1R3'}; % Bodypart Selective
%     grp(4).grids={'P3R0','P3R2','P5R0'}; % Place Selective
%     grp(5).grids={'P3R0','P3R2','P5R0'}; % Place Selective
% end

%%% Updated after Neuron Submission
if monkinitial=='S',
    monkeyname='Stewie'; sheetname='RSVP Cells_S';
    % Grid location groups for comparison
    grp(1).grids={'A7L2','A7L1','A6L3'}; % BodyPart Selective
    grp(2).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Face Selective
    grp(3).grids={'A4L2','A4L1','A4R1'}; % No Category Selectivity
    grp(4).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4', 'P4L2','P4L4','P5L3','P6L3'}; % Object Selective
    grp(5).grids={'P6L2','P6L1','P7L2'}; % Face Selective
    
    sampledlocations={'A7R1','A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'};
    % Grid location groups for comparison
    grpf(1).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Anterior Face Selective
    grpf(2).grids={'P6L2','P6L1','P7L2'}; % Posterior Face Selective
    grpf(3).grids={'A6L2','A6L0','A5L2','A5L1','A5L0','P6L2','P6L1','P7L2'}; % Inside All
    grpf(4).grids={'A7L2','A7L1','A7R1','A4L2','A4L1','A4R1','A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P5L3','P6L3'}; % Outside All Face Selective
    grpbp(1).grids={'A7L2','A7L1'}; % Inside Bodypart Selective
    grpbp(2).grids={'A7R1','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Outside All Bodypart Selective
    grpob(1).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L2'}; % Inside Object Selective
    grpob(2).grids={'A7R1','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','P5L0','P6L1','P6L3','P7L2'}; % Outside All Object Selective
    grppl(1).grids={'A4L2','A4L1','A4R1','A7R1'}; % Non Category Selective
    grppl(2).grids={'A7L2','A7L1','A6L2','A6L0','A5L2','A5L1','A5L0','A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L3','P6L2','P6L1','P7L2'}; % Category Selective
 
    grpfnf(1).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Anterior Face Selective
    grpfnf(2).grids={'A7L2','A7L1','A7R1','A4L2','A4L1','A4R1'}; % Near Face Selective (Anterior)
    grpfnf(3).grids={'P6L2','P6L1','P7L2'}; % Posterior Face Selective
    grpfnf(4).grids={'P4L2','P4L4','P5L2','P5L3','P6L3','P5L0'}; % Near Face Selective (Posterior)
    grpfnf(5).grids={'A7R1','A7L1','A7L2','A6L3','A4R1','A4L1','A4L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Outside Face Selective (All)
    grpbpnf(1).grids={'A7L2','A7L1','A6L3'}; % Inside Bodypart Selective
    grpbpnf(2).grids={'A7R1','A6L2','A6L0','A5L2','A5L1','A5L0','A4L2','A4L1','A4R1',}; % Near Bodypart Selective
    grpbpnf(3).grids={'A7R1','A6L0','A6L2','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Outside BodySelective
    grpobnf(1).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L2'}; % Inside Object Selective
    grpobnf(2).grids={'A4L3','A4L2','A4L1','A4R1','P6L3','P6L1','P7L2'}; % Near Object Selective
    grpobnf(3).grids={'A7R1','A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4R1','A4L1','A4L2','A4L3','P5L0','P6L1','P6L3','P7L2'}; % Outside All Object Selective
    grpplnf(1).grids={'A4L2','A4L1','A4R1','A7R1'}; % Non Category Selective
    grpplnf(2).grids={'A7L2','A7L1','A6L2','A6L0','A5L2','A5L1','A5L0','A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L3','P6L2','P6L1','P7L2'}; % Category Selective
    grpplnf(3).grids={'A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Category Selective

    
    % In Near Far
    grpfnf(1).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Anterior Face Selective
    grpfnf(2).grids={'A7L2','A7L1','A7R1','A4L2','A4L1','A4R1'}; % Near Face Selective (Anterior)
    grpfnf(3).grids={'P6L2','P6L1','P7L2'}; % Posterior Face Selective
    grpfnf(4).grids={'P4L2','P4L4','P5L2','P5L3','P6L3','P5L0'}; % Near Face Selective (Posterior)
    grpfnf(5).grids={'A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5'}; % Outside Face Selective (All)
    grpbpnf(1).grids={'A7L2','A7L1','A6L3'}; % Inside Bodypart Selective
    grpbpnf(2).grids={'A7R1','A6L2','A6L0','A5L2','A5L1','A5L0','A4L2','A4L1','A4R1',}; % Near Bodypart Selective
    grpbpnf(3).grids={'A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Outside BodySelective
    grpobnf(1).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L2'}; % Inside Object Selective
    grpobnf(2).grids={'A4L3','A4L2','A4L1','A4R1','P6L3','P6L1','P7L2'}; % Near Object Selective
    grpobnf(3).grids={'A7R1','A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','P5L0'}; % Outside All Object Selective
    grpplnf(1).grids={'A4L2','A4L1','A4R1','A7R1'}; % Non Category Selective
    grpplnf(2).grids={'A7L2','A7L1','A6L2','A6L0','A5L2','A5L1','A5L0','A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4','P4L2','P4L4','P5L3','P6L3','P6L2','P6L1','P7L2'}; % Category Selective
    grpplnf(3).grids={'A7L1','A7L2','A6L0','A6L2','A6L3','A5L0','A5L1','A5L2','A4L3','A2L5','A0L0','A0L2','A0L6','P1L1','P2L3','P3L4','P3L5','P4L2','P4L4','P5L0','P5L3','P6L1','P6L2','P6L3','P7L2'}; % Category Selective
elseif monkinitial=='W',
    monkeyname='Wiggum'; sheetname='RSVP Cells_W';
    % Grid location groups for comparison OLD
    %grp(1).grids={'A6R2','A5R0','A4R3'}; % Bodypart Selective
    %grp(2).grids={'AR0','A2R1','A2R3','A2R5'}; % Face Selective
    %grp(3).grids={'P1R0','P1R3'}; % Bodypart Selective
    %grp(4).grids={'P3R0','P3R2','P5R0'}; % Place Selective
    %grp(5).grids={'P3R0','P3R2','P5R0'}; % Place Selective
    % Grid location groups for comparison NEW
    
    % Current Locations:
    % 'A0R0','A1R0','A2R1','A2R3','A2R5','A3R0','A3R2','A4R3','A5R0','A6R2'
    % 'A7L1','P1R0','P1R3','P3R0','P3R2','P5R0','P6R0'
    %grp(1).grids={'A7L2','A7L1','A7L0','A7R0','A7R1','A6L2','A6L1','A6L0','A6R0','A6R1'}; % Face Selective
    %grp(2).grids={'A7R2','A6R2','A5R3','A5R4'}; % Object Selective
    %grp(3).grids={'A1L1','A1L0','A1R0','A1R1','A1R2','A0L1','A0L0','A0R0','A0R1','A0R2','P1R0','P1L0','P1R1','P1R2','P1R3'}; % Bodypart Selective
    %grp(4).grids={'P3L1','P3L0','P3R0','P3R1','P3R2','P4L0','P4R0','P4R1','P4R2','P4R3'}; % Face Selective
    %grp(5).grids={'P6L1','P6L0','P6R0','P6R1','P6R2','P7L1','P7L0','P7R0','P7R1','P7R2'}; % Object Selective
     
    % updated May 10, 2010
    grp(1).grids={'A7R0','A7L1','A6L1'}; %  Places
    grp(2).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Faces
    grp(3).grids={'P1R0','P1R3'}; % Bodyparts
    grp(4).grids={'P3R0','P3R2'}; % Objects
    grp(5).grids={'P7R0','P7R2'}; %  Face Selective
    
    sampledlocations={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'};
    % Grid location groups for comparison
    grpf(1).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Anterior Face Selective
    grpf(2).grids={'P7R0','P7R2'}; % Posterior Face Selective
    grpf(3).grids={'A0R0','A0R1','A1R0','A2R1','A3R0','P7R0','P7R2'}; % Inside All
    grpf(4).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R2','A2R3','A2R5','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0'}; % Outside Face Selective
    grpbp(1).grids={'P1R0','P1R3'}; % Inside Bodypart Selective
    grpbp(2).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Bodypart Selective
    grpob(1).grids={'P3R0','P3R2'}; % Inside Object Selective
    grpob(2).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Object Selective
    grppl(1).grids={'A7R0','A7L1','A6L1'}; % Inside Place Category Selective
    grppl(2).grids={'A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Place Selective

    grpfnf(1).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Anterior Face Selective
    grpfnf(2).grids={'A5R0','A4R3','A3R2','A2R3','A2R5','P1R0','P1R3'}; % Near Face Selective (Anterior)
    grpfnf(3).grids={'P7R0','P7R2'}; % Posterior Face Selective
    grpfnf(4).grids={'P4R1','P5R0','P6R0','P7R0','P7R2'}; % Near Face Selective (Posterior)
    grpfnf(5).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R2','A2R3','A2R5','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0'}; % Outside all Face Selective
    grpbpnf(1).grids={'P1R0','P1R3'}; % Inside Bodypart Selective
    grpbpnf(2).grids={'A2R1','A1R0','A0R0','A0R1','P3R0','P3R2'}; % Near Bodypart Selective
    grpbpnf(3).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Bodypart Selective 
    grpobnf(1).grids={'P3R0','P3R2'}; % Inside Object Selective
    grpobnf(2).grids={'A0R0','A0R1','P1R0','P1R3','P4R1','P5R0'}; % Near Object Selective
    grpobnf(3).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside Object Selective
    grpplnf(1).grids={'A7R0','A7L1','A6L1'}; % Inside Place Category Selective
    grpplnf(2).grids={'A6R2','A5R0','A4R3'}; % Near Place Category Selective
    grpplnf(3).grids={'A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Category Selective

    % In Near Far
    grpfnf(1).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Anterior Face Selective
    grpfnf(2).grids={'A5R0','A4R3','A3R2','A2R3','A2R5','P1R0','P1R3'}; % Near Face Selective (Anterior)
    grpfnf(3).grids={'P7R0','P7R2'}; % Posterior Face Selective
    grpfnf(4).grids={'P4R1','P5R0','P6R0'}; % Near Face Selective (Posterior)
    grpfnf(5).grids={'A7L1','A7R0','A6L1','A6R2','P3R0','P3R2'}; % Outside all Face Selective
    grpbpnf(1).grids={'P1R0','P1R3'}; % Inside Bodypart Selective
    grpbpnf(2).grids={'A2R1','A1R0','A0R0','A0R1','P3R0','P3R2'}; % Near Bodypart Selective
    grpbpnf(3).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R3','A2R5','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Bodypart Selective 
    grpobnf(1).grids={'P3R0','P3R2'}; % Inside Object Selective
    grpobnf(2).grids={'A0R0','A0R1','P1R0','P1R3','P4R1','P5R0'}; % Near Object Selective
    grpobnf(3).grids={'A7L1','A7R0','A6L1','A6R2','A5R0','A4R3','A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','P6R0','P7R0','P7R2'}; % Outside Object Selective    
    grpplnf(1).grids={'A7R0','A7L1','A6L1'}; % Inside Place Category Selective
    grpplnf(2).grids={'A6R2','A5R0','A4R3'}; % Near Place Category Selective
    grpplnf(3).grids={'A3R0','A3R2','A2R1','A2R3','A2R5','A1R0','A0R0','A0R1','P1R0','P1R3','P3R0','P3R2','P4R1','P5R0','P6R0','P7R0','P7R2'}; % Outside All Category Selective
end


fontsize_sml=7; fontsize_med=8; fontsize_lrg=9;
minunitnum=5; % minimum number of units for site to be included in colourmaps
%hmiconfig.printer=1;
disp('*******************************************************************')
disp('* rsvp_postNeuron.m - Generates figures listed under Project 1 in *')
disp('*     RSVP500_Outline.docx.                                       *')
disp('*     SPECIFIC FOR AFTER NEURON SUBMISSION.                       *')
disp('*******************************************************************')

disp('Loading data for all neurons...')
neurlabel='Both';
[data,numgrids,counts_matrix,allunits,unit_index,unitdata]=plx500_prepproject1data(hmiconfig,sheetname);
save([hmiconfig.rootdir,'rsvp500_project1',filesep,'Project1DataNT_',monkeyname,'_',neurtype,'.mat'],'data','unit_index','unitdata');


%%%%%%% GENERATE NEW FIGURES %%%%%%%
% Figure 1 - Methods
% No matlab needed

% Figure 2 - Example Neurons & Tsao-style figures (Neuron Summary)
% Neuron summary will consist of individual PCOLOR plots (a la Tsao et al., 2006), ordered according to location
% Also neuronal responses will be normalized to average firing rate?
disp('Figure 2 Pcolor Plots')
pdata=[];
pdata=unitdata.stimresponse(:,[1:20,61:100,41:60]); % rearrange individual stim responses to F,Bp,O,P
%pbase=unitdata.baseline(:,[1:20,61:100,41:60]);
pbase=unitdata.baseline;
% Step2 - Remove baseline and normalize responses (try to MAX)
for rr=1:size(pdata,1), pdatanobase(rr,:)=pdata(rr,:)-pbase(rr); end
for rr=1:size(pdata,1), pdatanorm(rr,:)=pdatanobase(rr,:)/max(pdatanobase(rr,:)); end

figure; clf; cla; % selectivity index histograms
set(gcf,'Units','Normalized','Position',[0.05 0.25 0.8 0.6]); set(gca,'FontName','Arial','FontSize',8);
subplot(2,7,1); % Excitatory Responses (sorted according to category)
f1=find(strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 );
f2=find(strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1);
f3=find(strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1);
f4=find(strcmp(unit_index.prefcat_excite_nofruit,'Places')==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1);
pointer=[f1;f2;f3;f4];
pcolor(1:80,1:length(pointer),pdatanorm(pointer,1:80))
shading flat
caxis([-.75 .75])
hold on
%plot([1 80],[length(f1) length(f1)],'w-','LineWidth',1)
%plot([1 80],[length([f1;f2]) length([f1;f2])],'w-','LineWidth',1)
%plot([1 80],[length([f1;f2;f3]) length([f1;f2;f3])],'w-','LineWidth',1)
%plot([1 80],[length([f1;f2;f3;f4]) length([f1;f2;f3;f4])],'w-','LineWidth',1)
title(['Excitatory Responses (Neurons) sorted according Category Preferences (n=',num2str(size(pdatanorm,1)),')'])
set(gca,'FontSize',7); axis off; axis ij;
bigpointer=[];
for gg=1:5,
    clear f1 f2 f3 f4 pointer
    f1=find(strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & ismember(unit_index.GridLoc,grp(gg).grids)==1);
    f2=find(strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & ismember(unit_index.GridLoc,grp(gg).grids)==1);
    f3=find(strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & ismember(unit_index.GridLoc,grp(gg).grids)==1);
    f4=find(strcmp(unit_index.prefcat_excite_nofruit,'Places')==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & ismember(unit_index.GridLoc,grp(gg).grids)==1);
    pointer=[f1;f2;f3;f4];
    subplot(2,7,1+gg)
    pcolor(1:80,1:length(pointer),pdatanorm(pointer,1:80))
    shading flat
    caxis([-.75 .75])
    hold on
    plot([1 80],[length(f1) length(f1)],'w-','LineWidth',1)
    plot([1 80],[length([f1;f2]) length([f1;f2])],'w-','LineWidth',1)
    plot([1 80],[length([f1;f2;f3]) length([f1;f2;f3])],'w-','LineWidth',1)
    plot([1 80],[length([f1;f2;f3;f4]) length([f1;f2;f3;f4])],'w-','LineWidth',1)
    title(['Grid ',num2str(gg),' n=',num2str(length(pointer))])
    set(gca,'FontSize',7); axis off; axis ij;
    bigpointer=[bigpointer;pointer];
    ylim([0 100]);
end
subplot(2,7,7); % Suppressed Responses (sorted according to category)
pcolor(1:80,1:length(bigpointer),pdatanorm(bigpointer,1:80))
shading flat
caxis([-.75 .75])
hold on
set(gca,'FontSize',7); axis off; axis ij;

subplot(2,7,8); % Suppressed Responses (sorted according to category)
f1=find(strcmp(unit_index.prefcat_inhibit_nofruit,'Faces')==1 & ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 );
f2=find(strcmp(unit_index.prefcat_inhibit_nofruit,'BodyParts')==1 & ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1);
f3=find(strcmp(unit_index.prefcat_inhibit_nofruit,'Objects')==1 & ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1);
f4=find(strcmp(unit_index.prefcat_inhibit_nofruit,'Places')==1 & ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1);
pointer=[f1;f2;f3;f4];
pcolor(1:80,1:length(pointer),pdatanorm(pointer,1:80))
shading flat
caxis([-.75 .75])
hold on
%plot([1 80],[length(f1) length(f1)],'w-','LineWidth',1)
%plot([1 80],[length([f1;f2]) length([f1;f2])],'w-','LineWidth',1)
%plot([1 80],[length([f1;f2;f3]) length([f1;f2;f3])],'w-','LineWidth',1)
%plot([1 80],[length([f1;f2;f3;f4]) length([f1;f2;f3;f4])],'w-','LineWidth',1)
title(['Suppressed Responses (Neurons) sorted according Category Preferences (n=',num2str(size(pdatanorm,1)),')'])
set(gca,'FontSize',7); axis off; axis ij;
bigpointer=[];
for gg=1:5,
    clear f1 f2 f3 f4 pointer
    f1=find(strcmp(unit_index.prefcat_inhibit_nofruit,'Faces')==1 & ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 & ismember(unit_index.GridLoc,grp(gg).grids)==1);
    f2=find(strcmp(unit_index.prefcat_inhibit_nofruit,'BodyParts')==1 & ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 & ismember(unit_index.GridLoc,grp(gg).grids)==1);
    f3=find(strcmp(unit_index.prefcat_inhibit_nofruit,'Objects')==1 & ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 & ismember(unit_index.GridLoc,grp(gg).grids)==1);
    f4=find(strcmp(unit_index.prefcat_inhibit_nofruit,'Places')==1 & ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 & ismember(unit_index.GridLoc,grp(gg).grids)==1);
    pointer=[f1;f2;f3;f4];
    subplot(2,7,8+gg)
    pcolor(1:80,1:length(pointer),pdatanorm(pointer,1:80))
    shading flat
    caxis([-.75 .75])
    hold on
    plot([1 80],[length(f1) length(f1)],'w-','LineWidth',1)
    plot([1 80],[length([f1;f2]) length([f1;f2])],'w-','LineWidth',1)
    plot([1 80],[length([f1;f2;f3]) length([f1;f2;f3])],'w-','LineWidth',1)
    plot([1 80],[length([f1;f2;f3;f4]) length([f1;f2;f3;f4])],'w-','LineWidth',1)
    title(['Grid ',num2str(gg),' n=',num2str(length(pointer))])
    set(gca,'FontSize',7); axis off; axis ij;
    bigpointer=[bigpointer;pointer];
    ylim([0 100]);
end
subplot(2,7,14); % Suppressed Responses (sorted according to category)
pcolor(1:80,1:length(bigpointer),pdatanorm(bigpointer,1:80))
shading flat
caxis([-.75 .75])
hold on
set(gca,'FontSize',7); axis off; axis ij;
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig2_PcolorPlots_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig2_PcolorPlots_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)


% Figure 3 - Descriptive Stats, independent of location (?)
disp('Figure 3 Descriptive Stats, Independent of Location')
figure; clf; cla; % selectivity index histograms
set(gcf,'Units','Normalized','Position',[0.05 0.25 0.8 0.6]); set(gca,'FontName','Arial','FontSize',8);
subplot(3,5,1); piedata=[]; % Neuron Distribution
piedata(1)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1));
piedata(2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1));
piedata(3)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1));
piedata(4)=length(find(strcmp(unit_index.excitetype_nofruit,'Non-Responsive')==1));
pie(piedata,...
    {['n=',num2str(piedata(1)),'(',num2str(piedata(1)/sum(piedata)*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(2)),'(',num2str(piedata(2)/sum(piedata)*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(3)),'(',num2str(piedata(3)/sum(piedata)*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(4)),'(',num2str(piedata(4)/sum(piedata)*100,'%1.2g'),'%)']})
title(['(A) Excite/Inhibit/Both (n=',num2str(sum(piedata)),')'],'FontSize',fontsize_med,'FontWeight','Bold')
legend('E','B','I','NS','Location','Best'); set(gca,'FontSize',7)
subplot(3,5,2); bardata=[]; % breakdown of category selectivity according to response type
bardata(1,1)=length(find(strcmp(unit_index.selective_nofruit,'Selective')==1 & strcmp(unit_index.excitetype_nofruit,'Excite')==1));
bardata(1,2)=length(find(strcmp(unit_index.selective_nofruit,'Not Selective')==1 & strcmp(unit_index.excitetype_nofruit,'Excite')==1));
bardata(2,1)=length(find(strcmp(unit_index.selective_nofruit,'Selective')==1 & strcmp(unit_index.excitetype_nofruit,'Both')==1));
bardata(2,2)=length(find(strcmp(unit_index.selective_nofruit,'Not Selective')==1 & strcmp(unit_index.excitetype_nofruit,'Both')==1));
bardata(3,1)=length(find(strcmp(unit_index.selective_nofruit,'Selective')==1 & strcmp(unit_index.excitetype_nofruit,'Inhibit')==1));
bardata(3,2)=length(find(strcmp(unit_index.selective_nofruit,'Not Selective')==1 & strcmp(unit_index.excitetype_nofruit,'Inhibit')==1));
bardata(1,3)=bardata(1,1)/sum(bardata(1,1:2)); bardata(1,4)=bardata(1,2)/sum(bardata(1,1:2));
bardata(2,3)=bardata(2,1)/sum(bardata(2,1:2)); bardata(2,4)=bardata(2,2)/sum(bardata(2,1:2));
bardata(3,3)=bardata(3,1)/sum(bardata(3,1:2)); bardata(3,4)=bardata(3,2)/sum(bardata(3,1:2));
bar(1:3,bardata(:,3:4),'stack')
text(1,.25,['n=',num2str(bardata(1,1))],'FontSize',6,'HorizontalAlignment','Center')
text(1,.75,['n=',num2str(bardata(1,2))],'FontSize',6,'HorizontalAlignment','Center')
text(2,.25,['n=',num2str(bardata(2,1))],'FontSize',6,'HorizontalAlignment','Center')
text(2,.75,['n=',num2str(bardata(2,2))],'FontSize',6,'HorizontalAlignment','Center')
text(3,.25,['n=',num2str(bardata(3,1))],'FontSize',6,'HorizontalAlignment','Center')
text(3,.75,['n=',num2str(bardata(3,2))],'FontSize',6,'HorizontalAlignment','Center')
title(['(B) Selective/Non-Selective (n=',num2str(sum(sum(bardata(:,1:2)))),')'],'FontSize',fontsize_med,'FontWeight','Bold')
legend('S','NS','Location','Best'); set(gca,'FontSize',7)
set(gca,'XTick',1:3,'XTickLabel',{'E','B','I'}); ylabel('% of Neurons'); axis square
subplot(3,5,3); piedata=[]; % pie chart, excitatory preferences, only those that are significant for their excite category
pointer=find(unitdata.stats_prefexcite_v_others_nofruit<0.05);
tempdata=unitdata.prefcat_excite_nofruit(pointer);
piedata(1)=length(find(strcmp(tempdata,'Faces')==1));
piedata(2)=length(find(strcmp(tempdata,'BodyParts')==1));
piedata(3)=length(find(strcmp(tempdata,'Objects')==1));
piedata(4)=length(find(strcmp(tempdata,'Places')==1));
pie(piedata,...
    {['n=',num2str(piedata(1)),'(',num2str((piedata(1)/sum(piedata))*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(2)),'(',num2str((piedata(2)/sum(piedata))*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(3)),'(',num2str((piedata(3)/sum(piedata))*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(4)),'(',num2str((piedata(4)/sum(piedata))*100,'%1.2g'),'%)']})
title(['(C'') Neurons with Sig Excitatory Selectivity (n=',num2str(sum(piedata)),')'],'FontSize',fontsize_med,'FontWeight','Bold')
legend('Fa','Bp','Ob','Pl','Location','Best'); set(gca,'FontSize',7)
subplot(3,5,4); piedata=[]; % pie chart, excitatory preferences, only those that are significant for their excite category
pointer=find(unitdata.stats_prefinhibit_v_others_nofruit<0.05);
tempdata=unitdata.prefcat_excite_nofruit(pointer);
piedata(1)=length(find(strcmp(tempdata,'Faces')==1));
piedata(2)=length(find(strcmp(tempdata,'BodyParts')==1));
piedata(3)=length(find(strcmp(tempdata,'Objects')==1));
piedata(4)=length(find(strcmp(tempdata,'Places')==1));
pie(piedata,...
    {['n=',num2str(piedata(1)),'(',num2str((piedata(1)/sum(piedata))*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(2)),'(',num2str((piedata(2)/sum(piedata))*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(3)),'(',num2str((piedata(3)/sum(piedata))*100,'%1.2g'),'%)'] ...
    ['n=',num2str(piedata(4)),'(',num2str((piedata(4)/sum(piedata))*100,'%1.2g'),'%)']})
title(['(D'') Neurons with Sig Suppressed Selectivity (n=',num2str(sum(piedata)),')'],'FontSize',fontsize_med,'FontWeight','Bold')
legend('Fa','Bp','Ob','Pl','Location','Best'); set(gca,'FontSize',7)

subplot(3,5,5); pmap=zeros(5,5); % both neurons
catnames={'Faces','BodyParts','Objects','Places'};
for y=1:4, % each column (inhibit responses)
    for x=1:4, % each row (excite responses)
        pmap(y,x)=length(find(strcmp(unit_index.selective_nofruit,'Selective')==1 & strcmp(unit_index.excitetype_nofruit,'Both')==1 & ...
            strcmp(unit_index.prefcat_excite_nofruit,catnames(x))==1 & strcmp(unit_index.prefcat_inhibit_nofruit,catnames(y))==1));
    end
end
pmap=[pmap;[0 0 0 0 0]]; pmap=[pmap,[0 0 0 0 0 0]'];
tt=sum(sum(pmap)); pmap2=pmap/tt;
pcolor(pmap2); shading flat; set(gca,'YDir','reverse');
axis square; set(gca,'CLim',[0 .20]); 
mp=colormap; mp(1,:)=[0.7529 0.7529 0.7529]; colormap(mp)
set(gca,'XTick',1.5:4.5,'XTickLabel',catnames,'YTick',1.5:4.5,'YTickLabel',catnames,'FontSize',7)
ylabel('Preferred Inhibited Category','fontsize',7);
xlabel('Preferred Excited Category','fontsize',7);
colorbar('SouthOutside','FontSize',6)
title(['(E) Breakdown of Excite/Inhibited Responses of BOTH Neurons (n=',num2str(tt),')'],'FontSize',fontsize_med,'FontWeight','Bold')
subplot(3,2,3); % Normalized Responses (Excitatory) (All Neurons)
r500pn_normalized_rsp(unit_index,unitdata,{'Excite','Both'},sampledlocations)
title('Normalized Responses (All Excitatory Responses')
subplot(3,2,4); % Normalized Responses (Suppressed) (All Neurons)
r500pn_normalized_rsp(unit_index,unitdata,{'Inhibit','Both'},sampledlocations)
title('Normalized Responses (All Suppressed Responses')
subplot(3,2,5); bardata=[]; % Average Category SI
pointer1=find(ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1);
pointer2=find(ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1);
pointer3=find(ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1);
pointer4=find(ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1);
[bardata(1,1),bardata(1,2)]=mean_sem(unitdata.cat_si_nofruit(pointer1,1));
[bardata(2,1),bardata(2,2)]=mean_sem(unitdata.cat_si_nofruit(pointer2,3));
[bardata(3,1),bardata(3,2)]=mean_sem(unitdata.cat_si_nofruit(pointer3,4));
[bardata(4,1),bardata(4,2)]=mean_sem(unitdata.cat_si_nofruit(pointer4,2));
pointer1a=find(ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Faces')==1);
pointer2a=find(ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'BodyParts')==1);
pointer3a=find(ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Objects')==1);
pointer4a=find(ismember(unit_index.excitetype_nofruit,{'Inhibit','Both'})==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Places')==1);
[bardata(5,1),bardata(5,2)]=mean_sem(unitdata.cat_si_nofruit(pointer1a,1));
[bardata(6,1),bardata(6,2)]=mean_sem(unitdata.cat_si_nofruit(pointer2a,3));
[bardata(7,1),bardata(7,2)]=mean_sem(unitdata.cat_si_nofruit(pointer3a,4));
[bardata(8,1),bardata(8,2)]=mean_sem(unitdata.cat_si_nofruit(pointer4a,2));
bar(1:8,bardata(:,1)); hold on; errorbar(1:8,bardata(:,1),bardata(:,2))
set(gca,'XTick',0.5:8.5,'XTickLabel',{'F','Bp','Ob','Pl','F','Bp','Ob','Pl'},'FontSize',7)
title(['(F) Average Category SI'],'FontSize',fontsize_med,'FontWeight','Bold')
ylim([-0.5 0.5])
[p,h]=ranksum(unitdata.cat_si_nofruit(pointer1,1),unitdata.cat_si_nofruit(pointer2,3));
text(1.5,0.45,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
[p,h]=ranksum(unitdata.cat_si_nofruit(pointer2,3),unitdata.cat_si_nofruit(pointer3,4));
text(2.5,0.45,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
[p,h]=ranksum(unitdata.cat_si_nofruit(pointer3,4),unitdata.cat_si_nofruit(pointer4,2));
text(3.5,0.45,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
[p,h]=ranksum(unitdata.cat_si_nofruit(pointer1,1),unitdata.cat_si_nofruit(pointer4,2));
text(2,0.48,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
[p,h]=ranksum(unitdata.cat_si_nofruit(pointer1a,1),unitdata.cat_si_nofruit(pointer2a,3));
text(5.5,-0.45,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
[p,h]=ranksum(unitdata.cat_si_nofruit(pointer2a,3),unitdata.cat_si_nofruit(pointer3a,4));
text(6.5,-0.45,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
[p,h]=ranksum(unitdata.cat_si_nofruit(pointer3a,4),unitdata.cat_si_nofruit(pointer4a,2));
text(7.5,-0.45,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
[p,h]=ranksum(unitdata.cat_si_nofruit(pointer1a,1),unitdata.cat_si_nofruit(pointer4a,2));
text(6,-0.48,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
subplot(3,2,6) % latency, for each pref cat (CatSel, NonCatSel, All) - excitatory only
r500pn_latency(unit_index,unitdata,sampledlocations);
title('Response Latency')
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig3_GeneralProperties_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig3_GeneralProperties_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)


% Figure 4 - Samples and Properties, per Patch
% Note: To generate average Spike Density Functions, use
% rsvp_postNeuron_avgSpden('MonkInitial')
disp('Figure 4 Properties, per Patch')
figure; clf; cla; % selectivity index histograms
set(gcf,'Units','Normalized','Position',[0.05 0.1 0.8 0.8]); set(gca,'FontName','Arial','FontSize',8);
for pp=1:5,
    subplot(5,2,(pp*2)-1)
    r500pn_patch_rsp(unit_index,unitdata,grp(pp).grids)
    title(['Normalized Responses - Patch #',num2str(pp)])
    subplot(5,2,(pp*2))
    r500pn_latency(unit_index,unitdata,grp(pp).grids);
    title(['Latencies (Excite) - Patch #',num2str(pp)])
end
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig4_PatchGeneralProperties_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig4_PatchProperties_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)


% Figure 5 - Category Proportion for each Patch
disp('Figure 5 Category Proportion for each Patch')
figure; clf; cla; set(gcf,'Units','Normalized','Position',[0.05 0.15 0.9 0.8]); set(gca,'FontName','Arial','FontSize',8);
for pp=1:5, % Excitatory Preferences
    % one loop per patch
    subplot(4,5,pp)
    [ns,nc,prop]=extractPropGrid_Excite(data,grp(pp).grids); 
    bar(prop);
    x2a=chi2_test(prop,[25 25 25 25]);
    x2=chi2_test(nc,ones(1,4)*(ns/4));
    set(gca,'FontName','Arial','FontSize',7,'XTick',1:4,'XTickLabel',{'F','Bp','Ob','Pl'})
    ylabel('Average % CatPref Neurons','FontSize',8); ylim([0 50]); axis square
    title(['Grid ',num2str(pp),', Excite/CatPrefs'],'FontSize',7,'FontWeight','Bold')
    text(3,45,['p(X2)=',num2str(x2,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center')
    text(1,40,['n=',num2str(nc(1))],'FontSize',7,'HorizontalAlignment','Center')
    text(2,40,['n=',num2str(nc(2))],'FontSize',7,'HorizontalAlignment','Center')
    text(3,40,['n=',num2str(nc(3))],'FontSize',7,'HorizontalAlignment','Center')
    text(4,40,['n=',num2str(nc(4))],'FontSize',7,'HorizontalAlignment','Center')
end
for pp=1:5, % Excitatory Preferences
    % one loop per patch
    subplot(4,5,pp+5)
    [ns,nc,prop]=extractPropGrid_Excite(data,grp(pp).grids); 
    pie(prop);
    set(gca,'FontName','Arial','FontSize',7); legend('F','Bp','Ob','Pl')
end
for pp=1:5, % Inhibited Preferences
     % one loop per patch
    subplot(4,5,pp+10)
    [ns,nc,prop]=extractPropGrid_Inhibit(data,grp(pp).grids); 
    bar(prop);
    x2=chi2_test(prop,[25 25 25 25]);
    set(gca,'FontName','Arial','FontSize',7,'XTick',1:4,'XTickLabel',{'F','Bp','Ob','Pl'})
    ylabel('Average % CatPref Neurons','FontSize',8); ylim([0 50]); axis square
    title(['Patch ',num2str(pp),', Inhibit/CatPrefs'],'FontSize',7,'FontWeight','Bold')
    text(3,45,['p(X2)=',num2str(x2,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center')
    text(1,40,['n=',num2str(nc(1))],'FontSize',7,'HorizontalAlignment','Center')
    text(2,40,['n=',num2str(nc(2))],'FontSize',7,'HorizontalAlignment','Center')
    text(3,40,['n=',num2str(nc(3))],'FontSize',7,'HorizontalAlignment','Center')
    text(4,40,['n=',num2str(nc(4))],'FontSize',7,'HorizontalAlignment','Center')
end
for pp=1:5, % Inhibited Preferences
     % one loop per patch
    subplot(4,5,pp+15)
    [ns,nc,prop]=extractPropGrid_Inhibit(data,grp(pp).grids); 
    pie(prop);
    set(gca,'FontName','Arial','FontSize',7); legend('F','Bp','Ob','Pl')
end
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig5_CatPrefperpatch_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig5_CatPrefperpatch_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)


% Figure 6 - Distribution of Category-Selective Neurons (In/Near/Out)
% What is the difference between face neurons IN the patch vs. OUT vs. NEAR the patch?
disp('Figure 6 Distribution of Category-Selective Neurons (In/Near/Out)')
figure; clf; cla; 
set(gcf,'Units','Normalized','Position',[0.05 0.15 0.9 0.8])
set(gca,'FontName','Arial','FontSize',8)
subplot(4,3,1); ns=[]; nc=[]; prop=[];
for pp=1:5, % Excitatory Preferences for faces
    % one loop per patch
    [ns(pp),nc(pp,:),prop(pp,:)]=extractPropGrid_Excite(data,grpfnf(pp).grids); 
end
bar(prop(:,1));
%x2=chi2_test(prop,[20 20 20 20 20]);
set(gca,'FontName','Arial','FontSize',7,'XTick',1:5,'XTickLabel',{'Ant(in)','Ant(out)','Post(in)','Post(out)','Far'})
ylabel('% Face Neurons','FontSize',8); ylim([0 85]);
% Another ChiSquare
expected=round(ns*(sum(nc(:,1))/sum(ns)));
observed=nc(:,1);
x2=chi2_test(observed,expected);
title({'Face Neurons vs. Location (Excitatory)',['x2=',num2str(x2,'%1.2g')]},'FontSize',10,'FontWeight','Bold')
%text(3,65,['p(X2)=',num2str(x2,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center')
text(1,60,['n=',num2str(nc(1,1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,60,['n=',num2str(nc(2,1))],'FontSize',7,'HorizontalAlignment','Center')
text(3,60,['n=',num2str(nc(3,1))],'FontSize',7,'HorizontalAlignment','Center')
text(4,60,['n=',num2str(nc(4,1))],'FontSize',7,'HorizontalAlignment','Center')
text(5,60,['n=',num2str(nc(5,1))],'FontSize',7,'HorizontalAlignment','Center')
text(1,65,[num2str(prop(1,1),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(2,65,[num2str(prop(2,1),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(3,65,[num2str(prop(3,1),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(4,65,[num2str(prop(4,1),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(5,65,[num2str(prop(5,1),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
%%% Excitatory Responses
subplot(4,3,2) % CatSI Analysis
% Compare ability of patches to discriminate faces vs. other categories
SI1=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpfnf(1).grids,'Faces');
SI2=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpfnf(2).grids,'Faces');
SI3=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpfnf(3).grids,'Faces');
SI4=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpfnf(4).grids,'Faces');
SI5=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpfnf(5).grids,'Faces');
hold on
bar([mean(SI1) mean(SI2) mean(SI3) mean(SI4) mean(SI5)])
errorbar(1:5,[mean(SI1) mean(SI2) mean(SI3) mean(SI4) mean(SI5)],[sem(SI1) sem(SI2) sem(SI3) sem(SI4) sem(SI5)])
set(gca,'FontName','Arial','FontSize',7,'XTick',1:5,'XTickLabel',{'Ant(in)','Ant(out)','Post(in)','Post(out)','Far'})
ylabel('Average Face SI','FontSize',8); ylim([0 .50]);
text(1,.38,['n=',num2str(length(SI1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,.38,['n=',num2str(length(SI2))],'FontSize',7,'HorizontalAlignment','Center')
text(3,.38,['n=',num2str(length(SI3))],'FontSize',7,'HorizontalAlignment','Center')
text(4,.38,['n=',num2str(length(SI4))],'FontSize',7,'HorizontalAlignment','Center')
text(5,.38,['n=',num2str(length(SI5))],'FontSize',7,'HorizontalAlignment','Center')
title({'SI Analysis (FaceNeuronsOnly) Excite',monkeyname},'FontWeight','Bold','FontSize',7);
try [p,h]=ranksum(SI1,SI2); text(1.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum(SI2,SI3); text(2.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum(SI3,SI4); text(3.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum(SI4,SI5); text(4.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum([SI1;SI2],[SI3;SI4]); text(2,0.48,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum([SI1;SI3],[SI5]); text(4,0.48,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
subplot(4,3,3) % stim select in face patches
SI1=extractStimSelect_Grid(unit_index,unitdata,grpfnf(1).grids);
SI2=extractStimSelect_Grid(unit_index,unitdata,grpfnf(2).grids);
SI3=extractStimSelect_Grid(unit_index,unitdata,grpfnf(3).grids);
SI4=extractStimSelect_Grid(unit_index,unitdata,grpfnf(4).grids);
SI5=extractStimSelect_Grid(unit_index,unitdata,grpfnf(5).grids);
catlabels={'Ant(in)','Ant(out)','Post(in)','Post(out)','Far'};
bardata=[SI1;SI2;SI3;SI4;SI5];
bar(bardata(:,1),'stack')
set(gca,'FontName','Arial','FontSize',8,'XTick',1:5,'XTickLabel',catlabels)
ylabel({'% Within Category Selectivity','(of Total Sensory Neurons)'},'FontSize',8); ylim([0 75]);
title({'Stimulus Selectivity per Patch',[char(catlabels(pp)),' - ',monkeyname]},'FontWeight','Bold','FontSize',7);
%for tl=1:5,
%    text(tl,sum(bardata(tl,1:2)),[num2str(bardata(tl,3),'%1.2g'),'%'])
%end
subplot(4,3,4); clear ns nc prop
for pp=1:3, % Excitatory Preferences for bodyparts
    % one loop per patch
    [ns(pp),nc(pp,:),prop(pp,:)]=extractPropGrid_Excite(data,grpbpnf(pp).grids); 
end
bar(prop(:,2));
%x2=chi2_test(prop,[20 20 20 20 20]);
set(gca,'FontName','Arial','FontSize',7,'XTick',1:3,'XTickLabel',{'In','Near','Far'})
ylabel('% Bodypart Neurons','FontSize',8); ylim([0 85]);

% Another ChiSquare
expected=[];
expected=round(ns*(sum(nc(:,1))/sum(ns)));
observed=nc(:,1);
x2=chi2_test(observed,expected);

title({'Bodypart Neurons vs. Location (Excitatory)',['x2=',num2str(x2,'%1.2g')]},'FontSize',10,'FontWeight','Bold')
%text(3,65,['p(X2)=',num2str(x2,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center')
text(1,60,['n=',num2str(nc(1,1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,60,['n=',num2str(nc(2,1))],'FontSize',7,'HorizontalAlignment','Center')
text(3,60,['n=',num2str(nc(3,1))],'FontSize',7,'HorizontalAlignment','Center')
text(1,65,[num2str(prop(1,2),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(2,65,[num2str(prop(2,2),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(3,65,[num2str(prop(3,2),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
%%% Excitatory Responses
subplot(4,3,5) % CatSI Analysis
% Compare ability of patches to discriminate faces vs. other categories
SI1=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[3],grpbpnf(1).grids,'BodyParts');
SI2=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[3],grpbpnf(2).grids,'BodyParts');
SI3=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[3],grpbpnf(3).grids,'BodyParts');
hold on
bar([mean(SI1) mean(SI2) mean(SI3)])
errorbar(1:3,[mean(SI1) mean(SI2) mean(SI3)],[sem(SI1) sem(SI2) sem(SI3)])
set(gca,'FontName','Arial','FontSize',7,'XTick',1:3,'XTickLabel',{'In','Near','Far'})
ylabel('Average Bodypart SI','FontSize',8); ylim([0 .50]);
text(1,.38,['n=',num2str(length(SI1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,.38,['n=',num2str(length(SI2))],'FontSize',7,'HorizontalAlignment','Center')
text(3,.38,['n=',num2str(length(SI3))],'FontSize',7,'HorizontalAlignment','Center')
title({'SI Analysis (BPNeuronsOnly) Excite',monkeyname},'FontWeight','Bold','FontSize',7);
try [p,h]=ranksum(SI1,SI2); text(1.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum(SI2,SI3); text(2.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum([SI1],[SI3]); text(2,0.48,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
subplot(4,3,6) % stim select in bodypart patches
SI1=extractStimSelect_Grid(unit_index,unitdata,grpbpnf(1).grids);
SI2=extractStimSelect_Grid(unit_index,unitdata,grpbpnf(2).grids);
SI3=extractStimSelect_Grid(unit_index,unitdata,grpbpnf(3).grids);
catlabels={'Body In','Body Near','Body Far'};
bardata=[SI1;SI2;SI3];
bar(bardata(:,2),'stack')
set(gca,'FontName','Arial','FontSize',8,'XTick',1:3,'XTickLabel',catlabels)
ylabel({'% Within Category Selectivity','(of Total Sensory Neurons)'},'FontSize',8); ylim([0 75]);
title({'Stimulus Selectivity per Patch',[char(catlabels(pp)),' - ',monkeyname]},'FontWeight','Bold','FontSize',7);
%for tl=1:3,
%    text(tl,sum(bardata(tl,1:2)),[num2str(bardata(tl,3),'%1.2g'),'%'])
%end
subplot(4,3,7); clear ns nc prop
for pp=1:3, % Excitatory Preferences for objects
    % one loop per patch
    [ns(pp),nc(pp,:),prop(pp,:)]=extractPropGrid_Excite(data,grpobnf(pp).grids); 
end
bar(prop(:,3));
%x2=chi2_test(prop,[20 20 20 20 20]);
set(gca,'FontName','Arial','FontSize',7,'XTick',1:3,'XTickLabel',{'In','Near','Far'})
ylabel('% Object Neurons','FontSize',8); ylim([0 85]);

% Another ChiSquare
expected=round(ns*(sum(nc(:,1))/sum(ns)));
observed=nc(:,1);
x2=chi2_test(observed,expected);
title({'Object Neurons vs. Location (Excitatory)',['x2=',num2str(x2,'%1.2g')]},'FontSize',10,'FontWeight','Bold')
%text(3,65,['p(X2)=',num2str(x2,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center')
text(1,60,['n=',num2str(nc(1,1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,60,['n=',num2str(nc(2,1))],'FontSize',7,'HorizontalAlignment','Center')
text(3,60,['n=',num2str(nc(3,1))],'FontSize',7,'HorizontalAlignment','Center')
text(1,65,[num2str(prop(1,3),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(2,65,[num2str(prop(2,3),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(3,65,[num2str(prop(3,3),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
%%% Excitatory Responses
subplot(4,3,8) % CatSI Analysis
% Compare ability of patches to discriminate faces vs. other categories
SI1=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[4],grpobnf(1).grids,'Objects');
SI2=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[4],grpobnf(2).grids,'Objects');
SI3=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[4],grpobnf(3).grids,'Objects');
hold on
bar([mean(SI1) mean(SI2) mean(SI3)])
errorbar(1:3,[mean(SI1) mean(SI2) mean(SI3)],[sem(SI1) sem(SI2) sem(SI3)])
set(gca,'FontName','Arial','FontSize',7,'XTick',1:3,'XTickLabel',{'In','Near','Far'})
ylabel('Average Object SI','FontSize',8); ylim([0 .50]);
text(1,.38,['n=',num2str(length(SI1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,.38,['n=',num2str(length(SI2))],'FontSize',7,'HorizontalAlignment','Center')
text(3,.38,['n=',num2str(length(SI3))],'FontSize',7,'HorizontalAlignment','Center')
title({'SI Analysis (ObjectNeuronsOnly) Excite',monkeyname},'FontWeight','Bold','FontSize',7);
try [p,h]=ranksum(SI1,SI2); text(1.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum(SI2,SI3); text(2.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum([SI1],[SI3]); text(2,0.48,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
subplot(4,3,9) % stim select in object patches
SI1=extractStimSelect_Grid(unit_index,unitdata,grpobnf(1).grids);
SI2=extractStimSelect_Grid(unit_index,unitdata,grpobnf(2).grids);
SI3=extractStimSelect_Grid(unit_index,unitdata,grpobnf(3).grids);
catlabels={'Object In','Object Near','Object Far'};
bardata=[SI1;SI2;SI3];
bar(bardata(:,3),'stack')
set(gca,'FontName','Arial','FontSize',8,'XTick',1:3,'XTickLabel',catlabels)
ylabel({'% Within Category Selectivity','(of Total Sensory Neurons)'},'FontSize',8); ylim([0 75]);
title({'Stimulus Selectivity per Patch',[char(catlabels(pp)),' - ',monkeyname]},'FontWeight','Bold','FontSize',7);
%for tl=1:3,
%    text(tl,sum(bardata(tl,1:2)),[num2str(bardata(tl,3),'%1.2g'),'%'])
%end
subplot(4,3,10); clear ns nc prop
for pp=1:3, % Excitatory Preferences for faces
    % one loop per patch
    [ns(pp),nc(pp,:),prop(pp,:)]=extractPropGrid_Excite(data,grpplnf(pp).grids); 
end
bar(prop(:,4));
%x2=chi2_test(prop,[20 20 20 20 20]);
set(gca,'FontName','Arial','FontSize',7,'XTick',1:3,'XTickLabel',{'In','Near','Far'})
ylabel('% Place Neurons','FontSize',8); ylim([0 85]); 

% Another ChiSquare
expected=round(ns*(sum(nc(:,1))/sum(ns)));
observed=nc(:,1);
x2=chi2_test(observed,expected);
title({'Place Neurons vs. Location (Excitatory)',['x2=',num2str(x2,'%1.2g')]},'FontSize',10,'FontWeight','Bold')
%text(3,65,['p(X2)=',num2str(x2,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center')
text(1,60,['n=',num2str(nc(1,1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,60,['n=',num2str(nc(2,1))],'FontSize',7,'HorizontalAlignment','Center')
text(3,60,['n=',num2str(nc(3,1))],'FontSize',7,'HorizontalAlignment','Center')
text(1,65,[num2str(prop(1,4),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(2,65,[num2str(prop(2,4),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
text(3,65,[num2str(prop(3,4),'%g'),'%'],'FontSize',7,'HorizontalAlignment','Center')
%%% Excitatory Responses
subplot(4,3,11) % CatSI Analysis
% Compare ability of patches to discriminate faces vs. other categories
SI1=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[2],grpplnf(1).grids,'Places');
SI2=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[2],grpplnf(2).grids,'Places');
SI3=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[2],grpplnf(3).grids,'Places');
hold on
bar([mean(SI1) mean(SI2) mean(SI3)])
errorbar(1:3,[mean(SI1) mean(SI2) mean(SI3)],[sem(SI1) sem(SI2) sem(SI3)])
set(gca,'FontName','Arial','FontSize',7,'XTick',1:3,'XTickLabel',{'In','Near','AllOut'})
ylabel('Average Place SI','FontSize',8); ylim([0 .50]);
text(1,.38,['n=',num2str(length(SI1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,.38,['n=',num2str(length(SI2))],'FontSize',7,'HorizontalAlignment','Center')
text(3,.38,['n=',num2str(length(SI3))],'FontSize',7,'HorizontalAlignment','Center')
title({'SI Analysis (PlaceNeuronsOnly) Excite',monkeyname},'FontWeight','Bold','FontSize',7);
try [p,h]=ranksum(SI1,SI2); text(1.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum(SI2,SI3); text(2.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum([SI1],[SI3]); text(2,0.48,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
subplot(4,3,12) % stim select in place patches
SI1=extractStimSelect_Grid(unit_index,unitdata,grpplnf(1).grids);
SI2=extractStimSelect_Grid(unit_index,unitdata,grpplnf(2).grids);
SI3=extractStimSelect_Grid(unit_index,unitdata,grpplnf(3).grids);
catlabels={'Place In','Place In','Places Far'};
bardata=[SI1;SI2;SI3];
bar(bardata(:,4),'stack')
set(gca,'FontName','Arial','FontSize',8,'XTick',1:2,'XTickLabel',catlabels)
ylabel({'% Within Category Selectivity','(of Total Sensory Neurons)'},'FontSize',8); ylim([0 75]);
title({'Stimulus Selectivity per Patch',[char(catlabels(pp)),' - ',monkeyname]},'FontWeight','Bold','FontSize',7);
%for tl=1:3,
%    text(tl,sum(bardata(tl,1:2)),[num2str(bardata(tl,3),'%1.2g'),'%'])
%end
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig6_DistInNearOut_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig6_DistInNearOut_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)



% Figure 7 - Inside vs. Outside, comparison
disp('Figure 7 Inside vs. Outside, comparison')
figure; clf; cla; set(gcf,'Units','Normalized','Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial','FontSize',8)
subplot(4,1,1); hold on; % average CSI (excitatory responses)
[bardata(1,1) bardata(1,2)]=mean_sem(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpf(3).grids,'Faces'));
[bardata(2,1) bardata(2,2)]=mean_sem(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpf(4).grids,'Faces'));
[bardata(3,1) bardata(3,2)]=mean_sem(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[3],grpbp(1).grids,'BodyParts'));
[bardata(4,1) bardata(4,2)]=mean_sem(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[3],grpbp(2).grids,'BodyParts'));
[bardata(5,1) bardata(5,2)]=mean_sem(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[4],grpob(1).grids,'Objects'));
[bardata(6,1) bardata(6,2)]=mean_sem(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[4],grpob(2).grids,'Objects'));
[bardata(7,1) bardata(7,2)]=mean_sem(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[2],grppl(1).grids,'Places'));
[bardata(8,1) bardata(8,2)]=mean_sem(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[2],grppl(2).grids,'Places'));
bardata(1,3)=ranksum(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpf(3).grids,'Faces'),extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpf(4).grids,'Faces'));
bardata(3,3)=ranksum(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[3],grpbp(1).grids,'BodyParts'),extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[3],grpbp(2).grids,'BodyParts'));
bardata(5,3)=ranksum(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[4],grpob(1).grids,'Objects'),extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[4],grpob(2).grids,'Objects'));
bardata(7,3)=ranksum(extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[2],grppl(1).grids,'Places'),extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[2],grppl(2).grids,'Places'));
bar(1:8,bardata(:,1)); hold on
errorbar(1:8,bardata(:,1),bardata(:,2))
ylabel('Average CSI'); xlabel('Location')
for cc=1:2:7,
    text(cc,0.45,['p=',num2str(bardata(cc,3),'%1.2g')])
end

subplot(4,1,2); bardata=[];
bardata(1,1)= length(find(ismember(unit_index.GridLoc,grpf(3).grids)==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
bardata(2,1)= length(find(ismember(unit_index.GridLoc,grpf(4).grids)==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
bardata(3,1)=length(find(ismember(unit_index.GridLoc,grpbp(1).grids)==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
bardata(4,1)=length(find(ismember(unit_index.GridLoc,grpbp(2).grids)==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
bardata(5,1)=length(find(ismember(unit_index.GridLoc,grpob(1).grids)==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
bardata(6,1)=length(find(ismember(unit_index.GridLoc,grpob(2).grids)==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
bardata(7,1)=length(find(ismember(unit_index.GridLoc,grppl(1).grids)==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));
bardata(8,1)=length(find(ismember(unit_index.GridLoc,grppl(2).grids)==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));

bardata(1,2)= length(find(ismember(unit_index.GridLoc,grpf(3).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
bardata(2,2)= length(find(ismember(unit_index.GridLoc,grpf(4).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
bardata(3,2)=length(find(ismember(unit_index.GridLoc,grpbp(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
bardata(4,2)=length(find(ismember(unit_index.GridLoc,grpbp(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
bardata(5,2)=length(find(ismember(unit_index.GridLoc,grpob(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
bardata(6,2)=length(find(ismember(unit_index.GridLoc,grpob(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
bardata(7,2)=length(find(ismember(unit_index.GridLoc,grppl(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));
bardata(8,2)=length(find(ismember(unit_index.GridLoc,grppl(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));

bardata(1,3)=bardata(1,1)/bardata(1,2)*100;
bardata(2,3)=bardata(2,1)/bardata(2,2)*100;
bardata(3,3)=bardata(3,1)/bardata(3,2)*100;
bardata(4,3)=bardata(4,1)/bardata(4,2)*100;
bardata(5,3)=bardata(5,1)/bardata(5,2)*100;
bardata(6,3)=bardata(6,1)/bardata(6,2)*100;
bardata(7,3)=bardata(7,1)/bardata(7,2)*100;
bardata(8,3)=bardata(8,1)/bardata(8,2)*100;
bar(1:8,bardata(:,3)); hold on
ylabel('% Across Category Selectivity'); xlabel('Location')

subplot(4,1,3); bardata=[];
bardata(1,1)= length(find(ismember(unit_index.GridLoc,grpf(3).grids)==1  & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & unitdata.anova_within_group(:,1)<0.05));
bardata(2,1)= length(find(ismember(unit_index.GridLoc,grpf(4).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & unitdata.anova_within_group(:,1)<0.05));
bardata(3,1)=length(find(ismember(unit_index.GridLoc,grpbp(1).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & unitdata.anova_within_group(:,4)<0.05));
bardata(4,1)=length(find(ismember(unit_index.GridLoc,grpbp(2).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & unitdata.anova_within_group(:,4)<0.05));
bardata(5,1)=length(find(ismember(unit_index.GridLoc,grpob(1).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & unitdata.anova_within_group(:,5)<0.05));
bardata(6,1)=length(find(ismember(unit_index.GridLoc,grpob(2).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & unitdata.anova_within_group(:,5)<0.05));
bardata(7,1)=length(find(ismember(unit_index.GridLoc,grppl(1).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & unitdata.anova_within_group(:,3)<0.05));
bardata(8,1)=length(find(ismember(unit_index.GridLoc,grppl(2).grids)==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & unitdata.anova_within_group(:,3)<0.05));

bardata(1,2)= length(find(ismember(unit_index.GridLoc,grpf(3).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
bardata(2,2)= length(find(ismember(unit_index.GridLoc,grpf(4).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
bardata(3,2)=length(find(ismember(unit_index.GridLoc,grpbp(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
bardata(4,2)=length(find(ismember(unit_index.GridLoc,grpbp(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
bardata(5,2)=length(find(ismember(unit_index.GridLoc,grpob(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
bardata(6,2)=length(find(ismember(unit_index.GridLoc,grpob(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
bardata(7,2)=length(find(ismember(unit_index.GridLoc,grppl(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));
bardata(8,2)=length(find(ismember(unit_index.GridLoc,grppl(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));

bardata(1,3)=bardata(1,1)/bardata(1,2)*100;
bardata(2,3)=bardata(2,1)/bardata(2,2)*100;
bardata(3,3)=bardata(3,1)/bardata(3,2)*100;
bardata(4,3)=bardata(4,1)/bardata(4,2)*100;
bardata(5,3)=bardata(5,1)/bardata(5,2)*100;
bardata(6,3)=bardata(6,1)/bardata(6,2)*100;
bardata(7,3)=bardata(7,1)/bardata(7,2)*100;
bardata(8,3)=bardata(8,1)/bardata(8,2)*100;
bar(1:8,bardata(:,3)); hold on
ylabel('% Within Category Selectivity'); xlabel('Location')

subplot(4,1,4); bardata=[];
[bardata(1,1) bardata(1,2)]=mean_sem(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpf(3).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1))));
[bardata(2,1) bardata(2,2)]=mean_sem(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpf(4).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1))));
[bardata(3,1) bardata(3,2)]=mean_sem(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpbp(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1))));
[bardata(4,1) bardata(4,2)]=mean_sem(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpbp(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1))));
[bardata(5,1) bardata(5,2)]=mean_sem(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpob(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1))));
[bardata(6,1) bardata(6,2)]=mean_sem(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpob(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1))));
[bardata(7,1) bardata(7,2)]=mean_sem(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grppl(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1))));
[bardata(8,1) bardata(8,2)]=mean_sem(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grppl(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1))));
bardata(1,3)=ranksum(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpf(3).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1))),...
    nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpf(4).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1))));
bardata(3,3)=ranksum(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpbp(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1))),...
    nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpbp(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1))));
bardata(5,3)=ranksum(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpob(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1))),...
    nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grpob(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1))));
bardata(7,3)=ranksum(nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grppl(1).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1))),...
    nonzeros(unitdata.latency(find(ismember(unit_index.GridLoc,grppl(2).grids)==1 & strcmp(unit_index.SensoryAuto,'Sensory')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1))));
bar(1:8,bardata(:,1)); hold on
errorbar(1:8,bardata(:,1),bardata(:,2))
ylabel('Average ROL'); xlabel('Location')
for cc=1:2:7,
    text(cc,145,['p=',num2str(bardata(cc,3),'%1.2g')])
end
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig7_CompareInOut_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig7_CompareInOut_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)




disp('Figure 8 Face Processing, Ant vs. Post comparison')
figure; clf; cla; set(gcf,'Units','Normalized','Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial','FontSize',8)
SI1=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpfnf(1).grids,'Faces');
SI2=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpfnf(3).grids,'Faces');
hold on
bar([mean(SI1) mean(SI2)])
errorbar(1:2,[mean(SI1) mean(SI2)],[sem(SI1) sem(SI2)])
set(gca,'FontName','Arial','FontSize',7,'XTick',1:2,'XTickLabel',{'Ant(in)','Post(in)'})
ylabel('Average Face SI','FontSize',8); ylim([0 .50]);
text(1,.38,['n=',num2str(length(SI1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,.38,['n=',num2str(length(SI2))],'FontSize',7,'HorizontalAlignment','Center')
title({'SI Analysis (FaceNeuronsOnly) Excite',monkeyname},'FontWeight','Bold','FontSize',7);
try [p,h]=ranksum(SI1,SI2); text(1.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig8_Faces_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig8_Faces_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)




%%% Table 1
f = figure('Position',[200 200 400 150]);
cnames = {'CatSelect','NotCatSelect','TOTAL'};
rnames = {'ExciteOnly','Faces','BodyParts','Objects','Places',...
    'SuppressedOnly','Faces','BodyParts','Objects','Places',...
    'BothExcite','Faces','BodyParts','Objects','Places',...
    'BothSuppress','Faces','BodyParts','Objects','Places',...
    'NonResponsive','TOTAL'};
tabledata=zeros(22,3);
% 
tabledata(1,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1));
tabledata(2,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
tabledata(3,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
tabledata(4,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
tabledata(5,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));
tabledata(1,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1));
tabledata(2,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
tabledata(3,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
tabledata(4,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
tabledata(5,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Excite')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));
tabledata(1,3)=sum(tabledata(1,1:2));
tabledata(2,3)=sum(tabledata(2,1:2));
tabledata(3,3)=sum(tabledata(3,1:2));
tabledata(4,3)=sum(tabledata(4,1:2));
tabledata(5,3)=sum(tabledata(5,1:2));

tabledata(6,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1));
tabledata(7,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Faces')==1));
tabledata(8,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'BodyParts')==1));
tabledata(9,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Objects')==1));
tabledata(10,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Places')==1));
tabledata(6,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1));
tabledata(7,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Faces')==1));
tabledata(8,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_inhibit_nofruit,'BodyParts')==1));
tabledata(9,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Objects')==1));
tabledata(10,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Inhibit')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Places')==1));
tabledata(6,3)=sum(tabledata(6,1:2));
tabledata(7,3)=sum(tabledata(7,1:2));
tabledata(8,3)=sum(tabledata(8,1:2));
tabledata(9,3)=sum(tabledata(9,1:2));
tabledata(10,3)=sum(tabledata(10,1:2));

tabledata(11,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1));
tabledata(12,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
tabledata(13,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
tabledata(14,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
tabledata(15,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));
tabledata(11,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1));
tabledata(12,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1));
tabledata(13,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_excite_nofruit,'BodyParts')==1));
tabledata(14,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_excite_nofruit,'Objects')==1));
tabledata(15,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_excite_nofruit,'Places')==1));
tabledata(11,3)=sum(tabledata(11,1:2));
tabledata(12,3)=sum(tabledata(12,1:2));
tabledata(13,3)=sum(tabledata(13,1:2));
tabledata(14,3)=sum(tabledata(14,1:2));
tabledata(15,3)=sum(tabledata(15,1:2));

tabledata(16,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1));
tabledata(17,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Faces')==1));
tabledata(18,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'BodyParts')==1));
tabledata(19,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Objects')==1));
tabledata(20,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')==1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Places')==1));
tabledata(16,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1));
tabledata(17,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Faces')==1));
tabledata(18,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_inhibit_nofruit,'BodyParts')==1));
tabledata(19,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Objects')==1));
tabledata(20,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Both')==1 & strcmp(unit_index.SelectiveConf,'Selective')~=1 & strcmp(unit_index.prefcat_inhibit_nofruit,'Places')==1));
tabledata(16,3)=sum(tabledata(16,1:2));
tabledata(17,3)=sum(tabledata(17,1:2));
tabledata(18,3)=sum(tabledata(18,1:2));
tabledata(19,3)=sum(tabledata(19,1:2));
tabledata(20,3)=sum(tabledata(20,1:2));

tabledata(21,1)=length(find(strcmp(unit_index.excitetype_nofruit,'Non-Responsive')==1));
tabledata(21,2)=length(find(strcmp(unit_index.excitetype_nofruit,'Non-Responsive')==1));
tabledata(21,3)=length(find(strcmp(unit_index.excitetype_nofruit,'Non-Responsive')==1));

tabledata(22,1)=sum(tabledata([1 6 11],1));
tabledata(22,2)=sum(tabledata([1 6 11],2));
tabledata(22,3)=sum(tabledata([1 6 11 21],3));


t = uitable(tabledata,cnames);


return

% Figure 8 - Face Processing - Anterior vs. Posterior Face-Patches
disp('Figure 8 Face Processing')
figure; clf; cla; set(gcf,'Units','Normalized','Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial','FontSize',8)
subplot(2,2,1); clear ns nc prop; hold on; % basic comparison = % face neurons
for pp=1:2, % Excitatory Preferences for faces
    % one loop per patch
    [ns(pp),nc(pp,:),prop(pp,:)]=extractPropGrid_Excite(data,grpf(pp).grids); 
end
for pp=1:2, % Excitatory Preferences for faces
    % one loop per patch
    [ns(pp+2),nc(pp+2,:),prop(pp+2,:)]=extractPropGrid_Inhibit(data,grpf(pp).grids); 
end
bar(prop(:,1));
set(gca,'FontName','Arial','FontSize',7,'XTick',1:4,'XTickLabel',{'Ant Excite','Post Excite','Ant Supp','Post Supp'})
ylabel('% Face Neurons','FontSize',8); ylim([0 75]); axis square
title(['Face Neuron Distribution (Ant vs. Post)'],'FontSize',10,'FontWeight','Bold')
text(1,60,['n=',num2str(nc(1,1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,60,['n=',num2str(nc(2,1))],'FontSize',7,'HorizontalAlignment','Center')
text(3,60,['n=',num2str(nc(3,1))],'FontSize',7,'HorizontalAlignment','Center')
text(4,60,['n=',num2str(nc(4,1))],'FontSize',7,'HorizontalAlignment','Center')
subplot(2,2,2); 
SI1=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpf(1).grids,'Faces');
SI2=extractCatSI_Grid_Excite_prefCat(unit_index,unitdata,[1],grpf(2).grids,'Faces');
SI3=extractCatSI_Grid_Inhibit_prefCat(unit_index,unitdata,[1],grpf(1).grids,'Faces');
SI4=extractCatSI_Grid_Inhibit_prefCat(unit_index,unitdata,[1],grpf(2).grids,'Faces');
hold on
bar([mean(SI1) mean(SI2) mean(SI3) mean(SI4)])
errorbar(1:4,[mean(SI1) mean(SI2) mean(SI3) mean(SI4)],[sem(SI1) sem(SI2) sem(SI3) sem(SI4)])
set(gca,'FontName','Arial','FontSize',8,'XTick',1:4,'XTickLabel',{'Ant Excite','Post Excite','Ant Supp','Post Supp'})
ylabel('Average Face SI','FontSize',8); ylim([-.50 .50]);
text(1,.38,['n=',num2str(length(SI1))],'FontSize',7,'HorizontalAlignment','Center')
text(2,.38,['n=',num2str(length(SI2))],'FontSize',7,'HorizontalAlignment','Center')
text(3,.38,['n=',num2str(length(SI3))],'FontSize',7,'HorizontalAlignment','Center')
text(4,.38,['n=',num2str(length(SI4))],'FontSize',7,'HorizontalAlignment','Center')
title({'SI Analysis (FaceNeuronsOnly)',monkeyname},'FontWeight','Bold','FontSize',7); axis square
try [p,h]=ranksum(SI1,SI2); text(1.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
try [p,h]=ranksum(SI3,SI4); text(3.5,0.44,['p=',num2str(p,'%1.2g')],'FontSize',7,'HorizontalAlignment','Center'); end
subplot(2,2,3) % stim select in face patches
SI1=extractStimSelect_Grid(unit_index,unitdata,grpf(1).grids);
SI2=extractStimSelect_Grid(unit_index,unitdata,grpf(2).grids);
catlabels={'Ant','Post'};
bardata=[SI1;SI2];
bar(1:2,bardata(:,1),'stack')
set(gca,'FontName','Arial','FontSize',8,'XTick',1:2,'XTickLabel',catlabels)
ylabel({'% Within Category Selectivity','(of Total Sensory Neurons)'},'FontSize',8); ylim([0 20]); axis square
title({'Stimulus Selectivity per Patch',[char(catlabels(pp)),' - ',monkeyname]},'FontWeight','Bold','FontSize',7);
subplot(2,2,4); hold on; clear bardata pointer1 pointer2 % latency
pointer1=find(ismember(unit_index.GridLoc,grpf(1).grids)==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1);
pointer2=find(ismember(unit_index.GridLoc,grpf(2).grids)==1 & ismember(unit_index.excitetype_nofruit,{'Excite','Both'})==1 & strcmp(unit_index.prefcat_excite_nofruit,'Faces')==1);
[bardata(1,1) bardata(1,2)]=mean_sem(nonzeros(unitdata.latency(pointer1,1:20)));
[bardata(2,1) bardata(2,2)]=mean_sem(nonzeros(unitdata.latency(pointer2,1:20)));
bardata(1,3)=length(pointer1); bardata(2,3)=length(pointer2);
bar(1:2,bardata(:,1))
errorbar(1:2,bardata(:,1),bardata(:,2))
set(gca,'XTick',1:2,'XTickLabel',{'Ant','Post'},'FontSize',7)
ylim([90 150]); axis square
for x=1:2,
    text(x,140,['n=',num2str(bardata(x,3))],'FontSize',7,'HorizontalAlignment','Center')
end
[p,h]=ranksum(nonzeros(unitdata.latency(pointer1,1:20)),nonzeros(unitdata.latency(pointer2,1:20)));
text(1.5,130,['p=',num2str(p,'%1.2g')],'FontSize',6,'HorizontalAlignment','Center')
ylabel('Average Response Latency (ms)');
jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig7_FaceProcessing_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_Fig7_FaceProcessing_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)




%%% NESTED FUNCTIONS %%%
function [numsensory,numcat,output]=extractPropGrid_Excite(data,gridlocs); % output will be multiple values (1/gridloc)
% Updated : No Fruit
catnames={'Face','Place','Bodypart','Object'};
numgrids=size(data,2); numsensory=0; numcat=0;
for gg=1:numgrids,
    if ismember(data(gg).gridloc(1,1),gridlocs)==1,
        numsensory=numsensory+sum(data(gg).counts_nofruit); numcat=numcat+data(gg).counts_nofruit;
    end
end
output=numcat/numsensory*100;
% resort
numcat=[numcat(1) numcat(3) numcat(4) numcat(2)]; output=[output(1) output(3) output(4) output(2)];
return

function [numsensory,numcat,output]=extractPropGrid_Inhibit(data,gridlocs); % output will be multiple values (1/gridloc)
% Updated : No Fruit
catnames={'Face','Place','Bodypart','Object'};
numgrids=size(data,2); numsensory=0; numcat=0;
for gg=1:numgrids,
    if ismember(data(gg).gridloc(1,1),gridlocs)==1,
        numsensory=numsensory+sum(data(gg).countsI_nofruit); numcat=numcat+data(gg).countsI_nofruit;
    end
end
output=numcat/numsensory*100;
% resort
numcat=[numcat(1) numcat(3) numcat(4) numcat(2)]; output=[output(1) output(3) output(4) output(2)];
return

function [numsensory,numcat,output]=extractPropGrid_Both(data,gridlocs); % output will be multiple values (1/gridloc)
% Updated : No Fruit
catnames={'Face','Place','Bodypart','Object'};
numgrids=size(data,2); numsensory=0; numcat=0;
for gg=1:numgrids,
    if ismember(data(gg).gridloc(1,1),gridlocs)==1,
        numsensory=numsensory+sum(data(gg).countsB_nofruit); numcat=numcat+data(gg).countsB_nofruit;
    end
end
output=numcat/numsensory*100;
% resort
numcat=[numcat(1) numcat(3) numcat(4) numcat(2)]; output=[output(1) output(3) output(4) output(2)];
return

function output=extractCatSI_Grid_Excite_prefCat(uindex,udata,catcol,gridlocs,prefcat)
pointer1=find(strcmp(uindex.SensoryConf,'Sensory')==1);
pointer2=find(ismember(uindex.excitetype_nofruit,{'Excite' 'Both'})==1);
pointer3=find(ismember(uindex.GridLoc,gridlocs)==1);
pointer4=find(strcmp(udata.prefcat_excite_nofruit,prefcat)==1);
pointerT1=intersect(pointer1,pointer2);
pointerT2=intersect(pointer3,pointer4);
pointer=intersect(pointerT1,pointerT2);
output=udata.cat_si_nofruit(pointer,catcol);
return

function output=extractCatSI_Grid_Inhibit_prefCat(uindex,udata,catcol,gridlocs,prefcat)
pointer1=find(strcmp(uindex.SensoryConf,'Sensory')==1);
pointer2=find(ismember(uindex.excitetype_nofruit,{'Inhibit' 'Both'})==1);
pointer3=find(ismember(uindex.GridLoc,gridlocs)==1);
pointer4=find(strcmp(udata.prefcat_inhibit_nofruit,prefcat)==1);
pointerT1=intersect(pointer1,pointer2);
pointerT2=intersect(pointer3,pointer4);
pointer=intersect(pointerT1,pointerT2);
output=udata.cat_si_nofruit(pointer,catcol);
return

function bardata=extractStimSelect_Grid(uindex,udata,gridlocs);
pointer1=find(strcmp(uindex.SensoryConf,'Sensory')==1);
pointer2=find(ismember(uindex.GridLoc,gridlocs)==1);
pointerT=intersect(pointer1,pointer2);
catnames={'Faces','BodyParts','Objects','Places'};
catcols=[1 4 5 3];
for cat=1:4,
    pointer3=find(strcmp(catnames(cat),udata.prefcat_excite_nofruit)==1|strcmp(catnames(cat),udata.prefcat_inhibit_nofruit)==1);
    pointer4=find(udata.anova_within_group(:,catcols(cat))<0.05);
    totalpointer=intersect(pointerT,pointer3);
    totalnum=length(totalpointer);
    stimnum=length(intersect(totalpointer,pointer4));
    bardata(cat)=(stimnum/totalnum)*100;
end
return

function bardata=extractStimSelect_Gridall(uindex,udata,gridlocs);
pointer1=find(strcmp(uindex.SensoryConf,'Sensory')==1);
pointer2=find(ismember(uindex.GridLoc,gridlocs)==1);
pointerT=intersect(pointer1,pointer2);
catnames={'Faces','BodyParts','Objects','Places'};
catcols=[1 4 5 3];
for ct=1:4,
    pointer3=find(strcmp(catnames(ct),udata.prefcat_excite_nofruit)==1|strcmp(catnames(ct),udata.prefcat_inhibit_nofruit)==1)'; % catselectneurons
    pointer4=find(udata.anova_within_group(:,catcols(ct))<0.05); % stim select
    bardata(ct,1)=length(intersect(pointerT,intersect(pointer4,pointer3))); % stim select
    bardata(ct,2)=length(intersect(pointerT,pointer3)); % total sensory select for cat
    bardata(ct,3)=length(pointerT); % total sensory
    bardata(ct,4)=(bardata(ct,1)/bardata(ct,3))*100;
    bardata(ct,5)=((bardata(ct,2)/bardata(ct,3))*100)-bardata(ct,4);
    bardata(ct,6)=bardata(ct,1)/bardata(ct,2)*100;
end
return

