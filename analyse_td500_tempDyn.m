%% analyse_td500_tempDyn.m
% Master file for temporal dynamics project, based on RSVP500 data
% in collaboration with Marieke Mur, UWO, Canada
% v.1.0 - Feb 16, 2024 - for use by KCL students

%% TASK DESCRIPTION:

clc; clear; close all;
dbstop if error;
global exptdata

program_version='2024-02-16';

if ispc
    rootdir='C:\Users\ahbel\OneDrive - King''s College London\ephysProjects';
else % MAC
    rootdir='~/OneDrive - King''s College London/ephysProjects';
end
addpath(userpath);
addpath(genpath([rootdir,filesep,'currentProjects']));
addpath(genpath([rootdir,filesep,'currentProjects',filesep,'proj_td500_tempDynamics']));
addpath(genpath([rootdir,filesep,'currentProjects',filesep,'commonProjectFunctions']));
addpath(genpath([rootdir,filesep,'Common_Functions']));

[td500specs.fList,td500specs.pList] = matlab.codetools.requiredFilesAndProducts('analyse_td500_tempDyn.m');

%exptdata=generate_td500_config; % generates and loads config file
exptdata.lastModified=datetime('today');
warning('off','MATLAB:MKDIR:DirectoryExists');
warning('off', 'MATLAB:table:ModifiedAndSavedVarnames')

% Modifiable defaults and directories (Change These)

exptdata.analysisName='TD500_TEMPDYN_Study'; % used for savenames, figures, etc. (pick whatever you want; will be used for filenames)


exptdata.projectdir=[exptdata.analysisdir,exptdata.analysisName,filesep]; mkdir(exptdata.projectdir); %


exptdata.figuredir500=[exptdata.projectdir,'figures',filesep,'td500_tempDyn',filesep]; mkdir(exptdata.figuredir500);
mkdir([exptdata.figuredir500,'matlabFigFiles']); mkdir([exptdata.figuredir500,'neuronPrintouts']);
diary([exptdata.projectdir,lower(exptdata.analysisName),'_',exptdata.analysisName,'.txt']);

% Preprocessing Parameters


% Analysis Parameters
exptdata.xrange_psths=-250:500; % window surrounding stimulus onset (in ms)
exptdata.reprocess=1; % recreate trialised files even if they already exist (preprocess) (includes prepping megaMatrices)
exptdata.reviewNeurons=1; % determines whether to generate figure for each neuron
exptdata.behav_reprocess=1; % set to one if you want to regenerate output files for existing files

% Save EXPTDATA structure
save([exptdata.analysisdir,filesep,exptdata.analysisName,'_exptdata.mat'],'exptdata');


% --------------------------------------------------------------------------------------------------------------------
%% STAGE 1: Compile and analyse behavioural data from both monkeys 
% (RUN THIS BEFORE NEURAL PREPROCESSING; does not require access to external drive)

% SETUP DEFAULTS
close all; dbstop if error
global NeuronData Grids monkeyname exptdata qualityPointer

%warning off

qualityCutoff=2; % select only neurons that have at least this grade

disp('+------------------------------------------------------------------------+')
disp('| td500_paper.m - Main analysis program for TEMPORAL/POPULATION DYNAMICS |')
disp('|                 analysis of RSVP dataset.                              |')
disp('+------------------------------------------------------------------------+')
disp(['Program Version: ',program_version])

exptdata.figurepath='~/Documents/_Current_Projects/TemporalDynamics500/figure_source_images';
cd(exptdata.figurepath)

% Add any additional fields to RESPSTRUCTSINGLE
%td500_AddFields('S',[1 2]); % add new fields
%td500_AddFields('W',[1 2]); % add new fields




cd(['~/Library/CloudStorage/OneDrive-King''sCollegeLondon/', 'MATLAB_Data/'])
load("rsvp500_massiveMatrix.mat")








% Compile Data
%td500_compiledata('Stewie','RSVP_Cells_S');
%td500_compiledata('Wiggum','RSVP_Cells_W');

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

