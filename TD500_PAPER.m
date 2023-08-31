function TD500_PAPER(analyses)
%%%%%%%%%%%%%%%
% TD500_PAPER %
%%%%%%%%%%%%%%%
% written by AHB, May 2014 
% Code adjusted to work under Mac OSX, Matlab v2013b
% Program to compile select RSVP500 neuronal data and perform TEMPORAL AND POPULATION DYNAMICS analysis
%
% Note: Program conventions:
%   all *.m files with lowercase letters are NESTED functions
%   all *.m files with UPPERCASE LETTERS are main programs

% ADD PATHS
addpath(userpath)
addpath(genpath('~/Documents/MATLAB/Common_Functions'))
addpath(genpath('~/Documents/MATLAB/chronux'));
addpath(genpath('~/Documents/MATLAB/td500_code'),'-BEGIN');

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