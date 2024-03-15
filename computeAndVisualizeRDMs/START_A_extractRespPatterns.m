function START_A_extractRespPatterns

% DESCRIPTION
% This script extracts response patterns for the ephys data described in 
% Bell et al. 2009 (https://doi.org/10.1523/JNEUROSCI.5865-10.2011). We
% extract response patterns from the spike-density functions (SPDs) 
% provided by Bell. We extract response patterns at millisecond resolution;
% the SPD data ranges from 0 to 5000 ms (1000 = stimulus onset).

% OUTPUT
% (1) Single-trial response patterns and associated neuron and trial info. 
% (2) Trial-average single-image response patterns, only including visually
%     responsive neurons for which we have trial info and AP index and only
%     including trials with nonzero mean and variance. We also save info 
%     for the included neurons and trials. 

% AUTHOR 
% Marieke Mur; last edit: 04-10-2015 (DD-MM-YYYY)


%% preparation
clear; close all;

dataPath='/imaging/ab03/PlexonData/rsvp500spks';
resultsPath='/imaging/mm07/projects/mITReprDynamics/analysis/results';
addpath(genpath('/imaging/mm07/programs/matlab/rsatoolbox')); % first MATLAB release of RSA toolbox
try mkdir(resultsPath); end


%% control variables
subjStr={'Stew' 'Wigg'};
nStimuli=100;


%% extract response patterns (use spike-density functions for now)
for subjectI=1:numel(subjStr)
    files=get_files(dataPath,[subjStr{subjectI},'*data.mat']);
    if strcmp(subjStr{subjectI},'Wigg'), files=files(1:end-6,:); end % last six files seem different (filename does not contain neuron info)
    nFiles=size(files,1);
    
    responseProp_str={'visResp' 'catSel' 'stimSel' 'excInh' 'quality'};
    responseProp__neuron_prop=zeros(nFiles/2,numel(responseProp_str));
    
    trialInfo_LOG=false(nFiles/2,1);
    APindex_LOG=false(nFiles/2,1);
    spd_LOG=false(nFiles/2,1);
    
    for fileI=1:nFiles
        
        if fileI==nFiles, disp('last file'); end
        
        if mod(fileI,2)==0 % even file number (= responsedata)
            fname=files(fileI,:);
            load(strtrim(fname));
            neuronI=fileI/2;
            fileIs_resp(neuronI)=fileI;
            
            % extract trial sequence
            if isfield(respstructsingle,'trial_id')
                trialInfo__neuron__trial_stimCat{neuronI}=respstructsingle.trial_id;
                trialInfo_LOG(neuronI)=true;
            else
                trialInfo__neuron__trial_stimCat{neuronI}=[];
                disp([files(fileI,:),': missing trial sequence']);
            end
            % some info on respstructsingle.trial_id (see Andrew's email):
            % column 1 = stimulus number for each trial
            %            1:20 faces
            %            21:40 fruit
            %            41:60 places
            %            61:80 bodyparts
            %            81:100 objects
            % column 2 = stimulus category
            %            1 - faces
            %            2 - fruit
            %            3 - places
            %            4 - bodyparts
            %            5 - objects
            
            % extract AP index
            if isfield(respstructsingle,'APIndex')
                APindex{neuronI}=respstructsingle.APIndex{1};
                APindex_LOG(neuronI)=true;
            else
                APindex{neuronI}=[];
                disp([files(fileI,:),': missing AP index']);
            end
            
            % extract response properties
            if isfield(respstructsingle,'conf_neurtype')
                % extract visual responsiveness (no/yes)
                if strcmp(respstructsingle.conf_neurtype,'Sensory'), responseProp__neuron_prop(neuronI,1)=1; end
                % extract category selectivity (no/faces/fruit/places/bodyparts/objects)
                if strcmp(respstructsingle.conf_preferred_cat,'Faces'), responseProp__neuron_prop(neuronI,2)=1;
                elseif strcmp(respstructsingle.conf_preferred_cat,'Fruit'), responseProp__neuron_prop(neuronI,2)=2;
                elseif strcmp(respstructsingle.conf_preferred_cat,'Places'), responseProp__neuron_prop(neuronI,2)=3;
                elseif strcmp(respstructsingle.conf_preferred_cat,'BodyParts'), responseProp__neuron_prop(neuronI,2)=4;
                elseif strcmp(respstructsingle.conf_preferred_cat,'Objects'), responseProp__neuron_prop(neuronI,2)=5;
                end
                % extract stimulus selectivity (no/yes)
                if strcmp(respstructsingle.conf_selective,'Selective'), responseProp__neuron_prop(neuronI,3)=1; end
                % extract excited/inhibited (exc/inh/both)
                if strcmp(respstructsingle.conf_excite,'Excite'), responseProp__neuron_prop(neuronI,4)=1;
                elseif strcmp(respstructsingle.conf_excite,'Inhibit'), responseProp__neuron_prop(neuronI,4)=2;
                elseif strcmp(respstructsingle.conf_excite,'Both'), responseProp__neuron_prop(neuronI,4)=3;
                end
                % extract quality assessment (0-3 - 3 is good?)
                responseProp__neuron_prop(neuronI,5)=respstructsingle.quality;
            else
                disp([files(fileI,:),': missing visual responsiveness']);
            end
            
        else % odd filenumber (= graphdata)
            
            % check whether pairs of files (graph and response data) refer to the same neuron
            filename_graphdata=files(fileI,:);
            filename_responsedata=files(fileI+1,:);
            hyphenIs=strfind(filename_graphdata,'-');
            if ~strcmp(filename_graphdata(1:hyphenIs(2)-1),filename_responsedata(1:hyphenIs(2)-1)),
                disp('WARNING: pairs of filenames do not seem to refer to the same neuron!');
            end
            
            % extract single-trial spike-density functions for each neuron
            fname=files(fileI,:);
            load(strtrim(fname));
            neuronI=ceil(fileI/2);
            fileIs_graph(neuronI)=fileI;
            spd__neuron__trial_timepoints{neuronI}=graphstructsingle.spden_trial;
            
            % graphstructsingle.spden_trial always exists but is sometimes empty (check that here)
            if ~isempty(graphstructsingle.spden_trial), spd_LOG(neuronI)=true;
            else disp([files(fileI,:),': missing spike-density functions']);
            end
            
        end
        
    end % fileI
    visResp_LOG=logical(responseProp__neuron_prop(:,1));    
    save(fullfile(resultsPath,['trialInfo_',subjStr{subjectI}]),'trialInfo__neuron__trial_stimCat','trialInfo_LOG','files');
    save(fullfile(resultsPath,['respProp_',subjStr{subjectI}]),'APindex','APindex_LOG','responseProp__neuron_prop','responseProp_str','visResp_LOG');
    save(fullfile(resultsPath,['SPDsingleTrial_',subjStr{subjectI}]),'spd__neuron__trial_timepoints','spd_LOG','-v7.3');
    
    % (1) select visually-responsive neurons for which we have trial information (stimulus and category) and AP index    
    neuronSelect_LOG=logical(floor((visResp_LOG+trialInfo_LOG+APindex_LOG+spd_LOG)/4));    
    spd__visRespNeuron__trial_timepoints=spd__neuron__trial_timepoints(neuronSelect_LOG);
    trialInfo__visRespNeuron__trial_stimCat=trialInfo__neuron__trial_stimCat(neuronSelect_LOG);    
    APindex_4RSA=APindex(neuronSelect_LOG);
    responseProp_4RSA__neuron_prop=responseProp__neuron_prop(neuronSelect_LOG,:);    
    
    % (2) select trials with nonzero mean and variance
    for neuronI=1:sum(neuronSelect_LOG)
        spd_cNeuron__trial_timepoints=spd__visRespNeuron__trial_timepoints{neuronI};
        trialInfo_cNeuron__trial_stimCat=trialInfo__visRespNeuron__trial_stimCat{neuronI};
 
        % make trial-selection vector
        spd_cNeuron__trial_timepointAvg=mean(spd_cNeuron__trial_timepoints,2);
        spd_cNeuron__trial_timepointVar=var(spd_cNeuron__trial_timepoints,0,2);
        avg_LOG=false(size(spd_cNeuron__trial_timepoints,1),1); avg_LOG(spd_cNeuron__trial_timepointAvg>0)=true;
        var_LOG=false(size(spd_cNeuron__trial_timepoints,1),1); var_LOG(spd_cNeuron__trial_timepointVar>0)=true;
        trialSelect_LOG=logical(floor((avg_LOG+var_LOG)/2));
        
        % keep track of number of (valid) trials for each stimulus
        for stimulusI=1:nStimuli
            trialIs=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimulusI);
            trialSelect_LOG_cStim=trialSelect_LOG(trialIs);           
            nTrials__neuron_stim(neuronI,stimulusI)=numel(trialIs);
            nValidTrials__neuron_stim(neuronI,stimulusI)=sum(trialSelect_LOG_cStim);
        end % stimulusI
        
        % select trials and assemble data for saving
        trialSelect_LOG__neuron{neuronI}=trialSelect_LOG;                
        spd_4RSA__neuron__trial_timepoints{neuronI}=spd_cNeuron__trial_timepoints(trialSelect_LOG,:);
        trialInfo_4RSA__neuron__trial_stimCat{neuronI}=trialInfo_cNeuron__trial_stimCat(trialSelect_LOG,:);                        
    end % neuronI  
    
    % save selected data    
    save(fullfile(resultsPath,['dataSelectionInfo_4RSA_',subjStr{subjectI}]),'neuronSelect_LOG','trialSelect_LOG__neuron','nTrials__neuron_stim','nValidTrials__neuron_stim'); 
    save(fullfile(resultsPath,['trialInfo_4RSA_',subjStr{subjectI}]),'trialInfo_4RSA__neuron__trial_stimCat');
    save(fullfile(resultsPath,['respProp_4RSA_',subjStr{subjectI}]),'APindex_4RSA','responseProp_4RSA__neuron_prop','responseProp_str');
    save(fullfile(resultsPath,['SPDsingleTrial_4RSA_',subjStr{subjectI}]),'spd_4RSA__neuron__trial_timepoints','-v7.3');
        
    % average spd functions across repetitions of the same stimulus (use selected data)
    for neuronI=1:sum(neuronSelect_LOG)
        spd_cNeuron__trial_timepoints=spd_4RSA__neuron__trial_timepoints{neuronI};
        trialInfo_cNeuron__trial_stimCat=trialInfo_4RSA__neuron__trial_stimCat{neuronI};
        
        for stimulusI=1:nStimuli
            % select single-trial spike-density functions for current stimulus
            trialIs=find(trialInfo_cNeuron__trial_stimCat(:,1)==stimulusI);
            if isempty(trialIs)
                disp('no data for current stimulus');
            end
            spd_cNeuron_cStim__trial_timepoints=spd_cNeuron__trial_timepoints(trialIs,:);
            % average spike-density functions across trials
            spdTrialAvg_cNeuron_cStim__timepoints=mean(spd_cNeuron_cStim__trial_timepoints,1);
            spdTrialAvg_4RSA__neuron_stim_timepoints(neuronI,stimulusI,:)=spdTrialAvg_cNeuron_cStim__timepoints;
        end % stimulusI
    
    end % neuronI    
    save(fullfile(resultsPath,['SPDexemplar_4RSA_',subjStr{subjectI}]),'spdTrialAvg_4RSA__neuron_stim_timepoints','-v7.3');
            
    % clear variables
    clear spd__neuron__trial_timepoints trialInfo__neuron__trial_stimCat APindex responseProp__neuron_prop files fileIs_resp fileIs_graph
    clear spdTrialAvg__neuron_stim_timepoints nTrials__neuron_stim nValidTrials__neuron_stim
    clear visResp_LOG trialInfo_LOG APindex_LOG neuronSelect_LOG nTrials__neuron_stim nValidTrials__neuron_stim            
    clear trialSelect_LOG__neuron spd_4RSA__neuron__trial_timepoints trialInfo_4RSA__neuron__trial_stimCat
    clear spdTrialAvg_4RSA__neuron_stim_timepoints
    
end % subjectI



%% get files (authors: Russell Thompson, Ian Charest, MRC CBU)
function files = get_files(direc, filt)
% =========================================================================
% return a list of files
% filt = filter string
% direc = cell array of directory names
if nargin~=2, error('get_files:missing inputs, Please input folder(s) and file filter.'); end%if
files = [];
if ischar(direc) % if direc is already a character array
    currDir = direc;
    tmp = dir(fullfile(currDir,filt)); % find all files matching f*.nii
    tmp = [repmat([currDir filesep],size(tmp,1),1) char(tmp.name)]; % build the full path name for these files
    files = char(files,tmp);
else % if direc is a cell array
    if size(direc,1)>size(direc,2)
        nRuns=size(direc,1);
    else
        nRuns=size(direc,2);
    end
    for runI=1:nRuns % loop through each EPI session
        currDir = char(direc{runI});
        tmp = dir(fullfile(currDir,filt)); % find all files matching f*.nii
        tmp = [repmat([currDir filesep],size(tmp,1),1) char(tmp.name)]; % build the full path name for these files
        files = char(files,tmp);
    end
end
files = files(~all(files'==' ')',:);


