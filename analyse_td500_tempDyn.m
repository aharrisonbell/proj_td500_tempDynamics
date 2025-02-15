%% analyse_td500_tempDyn.m
% Master file for temporal dynamics project, based on RSVP500 data
% in collaboration with Marieke Mur, UWO, Canada
% v.1.0 - Feb 16, 2024 - for use by KCL students
% v.1.1 - updated Sept 2024
clc; clear; close all;
program_version='1.1.2024-09-26';

%% TASK DESCRIPTION:

%% SETUP DEFAULTS
dbstop if error;
global exptdata
if ispc
    rootdir='C:\Users\ahbel\OneDrive - King''s College London\ephysProjects';
else % MAC
    rootdir='~/OneDrive - King''s College London/ephysProjects';
end
addpath(genpath([rootdir,filesep,'T500_TemporalDynamics_Study',filesep,'proj_td500_tempDynamics']));
addpath(genpath([rootdir,filesep,'commonProjectFunctions']));

[td500specs.fList,td500specs.pList] = matlab.codetools.requiredFilesAndProducts('analyse_td500_tempDyn.m');

exptdata=generate_td500_config; % generates and loads config file
exptdata.lastModified=datetime('today');
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

% Preprocessing Parameters (most established by generate_td500_config.m)
qualityCutoff=2; % select only neurons that have at least this grade

% Analysis Parameters (many established by generate_td500_config.m)
exptdata.xrange_psths=-250:500; % window surrounding stimulus onset (in ms)
exptdata.reprocess=1; % recreate trialised files even if they already exist (preprocess) (includes prepping megaMatrices)
exptdata.reviewNeurons=1; % determines whether to generate figure for each neuron
exptdata.behav_reprocess=1; % set to one if you want to regenerate output files for existing files

% Save EXPTDATA structure
save([exptdata.projectdir,filesep,exptdata.analysisName,'_exptdata.mat'],'exptdata');
diary([exptdata.projectdir,lower(exptdata.analysisName),'_',exptdata.analysisName,'.txt']);


% --------------------------------------------------------------------------------------------------------------------
%% STAGE 1: Compile data from both monkeys 
disp('+------------------------------------------------------------------------+')
disp('| td500_paper.m - Main analysis program for TEMPORAL/POPULATION DYNAMICS |')
disp('|                 analysis of RSVP/FACES datasets.                       |')
disp('+------------------------------------------------------------------------+')
disp(['Program Version: ',program_version])

load([exptdata.projectdir, filesep, 'rsvp500_massiveMatrix.mat'])
% load([exptdata.analysisdir, filesep, 'faces570_massiveMatrix.mat'])


% Add any additional fields to RESPSTRUCTSINGLE
%td500_AddFields('S',[1 2]); % add new fields
%td500_AddFields('W',[1 2]); % add new fields













% Load Data
load([exptdata.datadir,'td500data_Stewie.mat']);
S_xldata=UnitData; S_grids=grids;
load([exptdata.datadir,'td500data_Wiggum.mat']);
W_xldata=UnitData; W_grids=grids;

% Combine structs
NeuronData=mergeStructs(S_xldata, W_xldata);
Grids=mergeStructs(S_grids,W_grids);
save([exptdata.datadir,'td500data_BothMonkeys.mat'],'NeuronData','grids');

% Filter Neurons
qualityPointer=find(NeuronData.quality>=qualityCutoff & ismember(NeuronData.confRespDir,{'Excite','Both'}));

% Generate Figures - See td500_paper.docx for details

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1.1 - Latency 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(find(ismember(analyses,[1])==1,1))==0,
    td500_latency;
end
disp(monkeyname)

return

cd(['~/Library/CloudStorage/OneDrive-King''sCollegeLondon/', 'MATLAB_Data/'])
load("rsvp500_massiveMatrix.mat")





% SETUP DEFAULTS
close all; dbstop if error
global NeuronData Grids monkeyname lsnconfig qualityPointer
lsnconfig=generate_td500_config; % generates and loads config file
%warning off
if nargin==0, analyses=0:9; end
program_version='2014-05-01';
qualityCutoff=2; % select only neurons that have at least this grade

disp('+------------------------------------------------------------------------+')
disp('| td500_paper.m - Main analysis program for TEMPORAL/POPULATION DYNAMICS |')
disp('|                 analysis of RSVP dataset.                              |')
disp('+------------------------------------------------------------------------+')
disp(['Program Version: ',program_version])

lsnconfig.figurepath='~/Documents/_Current_Projects/TemporalDynamics500/figure_source_images';
cd(lsnconfig.figurepath)

% Add any additional fields to RESPSTRUCTSINGLE
%td500_AddFields('S',[1 2]); % add new fields
%td500_AddFields('W',[1 2]); % add new fields

% Compile Data
%td500_compiledata('Stewie','RSVP_Cells_S');
%td500_compiledata('Wiggum','RSVP_Cells_W');

% Load Data
load([lsnconfig.datadir,'td500data_Stewie.mat']);
S_xldata=UnitData; S_grids=grids;
load([lsnconfig.datadir,'td500data_Wiggum.mat']);
W_xldata=UnitData; W_grids=grids;

% Combine structs
NeuronData=mergeStructs(S_xldata, W_xldata);
Grids=mergeStructs(S_grids,W_grids);
save([lsnconfig.datadir,'td500data_BothMonkeys.mat'],'NeuronData','grids');

% Filter Neurons
qualityPointer=find(NeuronData.quality>=qualityCutoff & ismember(NeuronData.confRespDir,{'Excite','Both'}));

% Generate Figures - See td500_paper.docx for details

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Figure 1.1 - Latency 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if isempty(find(ismember(analyses,[1])==1,1))==0,
    td500_latency;
end
disp(monkeyname)

return