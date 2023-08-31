function plx500(files);
%%%%%%%%%%%%%%%%%%
% plx500(files); %
%%%%%%%%%%%%%%%%%%
% written by AHB, Sept2007,
% substantially updated April 2008 to include more rigorous methods of
% classifcation.
% Analyzes data for RSVP500 (rsvp500f.tim, rsvp500s.tim) task
% Incoming files must be run through plx_readnexfile and plx_makespikemat
% before analysis is possible
% files = optional argument, list files as strings.  Otherwise, program
% will load files listed in default analyze500.txt

%%% SETUP DEFAULTS
syncexcel=0;
warning off;
close all
hmiconfig=generate_hmi_configplex; % generates and loads config file
parnumlist=[500]; % list of paradigm numbers
xscale=-200:400; % default time window
metric_col=2; % column metric (1=spk counts, 2=mean spden, 3=peak spden, 4=mean peak, 5=area under curve)
anova=1; % 1=perform anovas, ~1=skip anova analysis
mpwin=[-30 30]; % mean peak window
minlatency=50; % minimum latencies

%%% CURRENT METRICS (used to define "response")
%% to add more epochs, make changes to RESPSTRUCT
baseline=[-200 50]; % window over which baseline response is calculated
epoch1=[50 300]; % early cue response window
epoch2=[100 400]; % middle cue response window
epoch3=[50 200]; % late cue response window
area_baseline=[-200 50];
area_epoch1=[50 300];

%% spot checked 152 neurons using 5 methods (Stew01/08-02/08) shows that
%% most neurons don't care but the best method is mean spden (50-250ms/50-300).
%% Peak measures don't deal with inhibited responses well at all


%%%  LOAD FILE LIST
if nargin==0,
    error('You must specify an individual filename or monkey initial (''S''/''W'').')
elseif strcmp(files,'S')==1
    disp('Analyzing all RSVP500 files for Stewie...')
    % Pulls files from HMI_PhysiologyNotes
    include=xlsread(hmiconfig.excelfile,'RSVP500','A9:A1000'); % alphanumeric, Gridlocation
    [crap,filest]=xlsread(hmiconfig.excelfile,'RSVP500','B9:B1000');
    filesx=filest(find(include==1)); clear include; clear files
    for ff=1:size(filesx,1),
        temp=char(filesx(ff)); files(ff)=cellstr(temp(1:12));
    end
elseif strcmp(files,'W')==1
    disp('Analyzing all RSVP500 files for Wiggum...')
    % Pulls files from HMI_PhysiologyNotes
    include=xlsread(hmiconfig.excelfile,'RSVP500','G9:G1000'); % alphanumeric, Gridlocation
    [crap,filest]=xlsread(hmiconfig.excelfile,'RSVP500','H9:H1000');
    filesx=filest(find(include==1)); clear include; clear files
    for ff=1:size(filesx,1),
        temp=char(filesx(ff));
        files(ff)=cellstr(temp(1:12));
    end
end
%%% ANALYZE INDIVIDUAL FILES
disp('*********************************************************************')
disp('plx500.m - Analysis program for RSVP500-series datafiles (April 2008)')
disp('*********************************************************************')
for f=1:length(files), % perform following operations on each nex file listed
    close all
    filename=char(files(f));
    %%% identify sheet name
    if filename(1)=='S', sheetname='RSVP Cells_S';
    elseif filename(1)=='W', sheetname='RSVP Cells_W';
    end
    disp('Removing previous files...')
    % remove previous files
    %killfiles=dir([hmiconfig.rsvp500spks,filename,'*-500*data.mat']); % graphstructs
    %for kf=1:size(killfiles,1),
    %    disp(['...deleting ',killfiles(kf).name])
    %    delete([hmiconfig.rsvp500spks,killfiles(kf).name]);
    %end
    killfilesfig=dir([hmiconfig.figure_dir,'rsvp500',filesep,filename,'-*.*']); % figures
    for kf=1:size(killfilesfig,1),
        disp(['...deleting ',killfilesfig(kf).name])
        delete([hmiconfig.figure_dir,'rsvp500',filesep,killfilesfig(kf).name]);
    end
    
    %plx_processnexfile({filename},0);
    
    %%% Test to see if plx_processnexfile has been run...
    disp('Finding spike matrix...')
    if exist([hmiconfig.spikedir,filename,'_spkmat.mat'])==2,
        disp('...Found! Continuing analysis...')
    else
        disp('...File Not Found! Attempting to process nexfile...')
        try
            plx_processnexfile({filename},0);
        catch
            error('...Unable to processnexfile.  Please mark this file to continue...')
        end
    end
    
    disp(['Analyzing spike activity from ',filename])
    tempstruct=load([hmiconfig.spikedir,filename,'_spkmat.mat']);
    tempbehav=tempstruct.behav_matrix(:,[1 3 4 30 40 44]); % load behavioural data
    tempbehav(:,7)=tempbehav(:,6)-tempbehav(:,5); % solve for cue onset time (aligned to the beginning of each trial, in ms?)
    tempspike=tempstruct.spikesig;
    clear tempstruct
    foundunits=size(tempspike,2);
    if length(find(ismember(tempbehav(:,2),parnumlist)))<1,
        disp(['..No RSVP500 trials found!!  Skipping this file.'])
    else
        disp(['..found ',num2str(size(tempbehav,1)),' trials...'])
        disp(['..found ',num2str(foundunits),' unit(s)...'])
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%% Unit loop
        for un=1:foundunits, % performed for each unit
            disp(['....analyzing ',char(tempspike(un).labels)])
            %%% setup structures
            spikestructsingle=struct('label',[],'faces_spk',[],'fruit_spk',[],'bodyp_spk',[],'places_spk',[],'objct_spk',[],...
                'faces_ts',[],'fruit_ts',[],'bodyp_ts',[],'places_ts',[],'objct_ts',[]);
            %%% check to see if graphstructsingle and respstructsingle
            %%% already exist for this unit
            unitname=char(tempspike(un).labels);
            if exist([hmiconfig.rsvp500spks,unitname(1:end-4),'-500graphdata.mat'])~=2, % doesn't exist
                graphstructsingle=struct('label',[],'faces_avg',[],'faces_sem',[],'fruit_avg',[],'fruit_sem',[],'bodyp_avg',[],...
                    'bodyp_sem',[],'places_avg',[],'places_sem',[],'objct_avg',[],'objct_sem',[],'allconds',[],'allconds_avg',[],'allconds_sem',[],...
                    'bestconds',[],'worstconds',[],'baseline',[],'cueresponse',[],'spden_trial',[]);
            else
                disp('......graphstructsingle already exists.  Loading previously saved data...')
                load([hmiconfig.rsvp500spks,unitname(1:end-4),'-500graphdata.mat'])
            end
            if exist([hmiconfig.rsvp500spks,unitname(1:end-4),'-500responsedata.mat'])~=2, % doesn't exist
                respstructsingle=struct('label',[],'spk_baseline',[],'m_baseline',[],'p_baseline',[],'trial_m_baseline',[],'trial_p_baseline',[],...
                    'anova_latency',[],'anova_baseline',[],'anova_epoch',[],'anova_within_group',[],...
                    'latency',[],'validrsp',[],'cat_avg',[],'cat_bst',[],'cat_sensory',[],...
                    'raw_si',[],'face_trad',[],'pairwise',[],'preferred_category',[],'preferred_sensory',[]);
            else
                disp('......respstructsingle already exists.  Loading previously saved data...')
                load([hmiconfig.rsvp500spks,unitname(1:end-4),'-500responsedata.mat'])
            end
            graphstructsingle.label=tempspike(un).labels; % paste label into GRAPH structure
            respstructsingle.label=tempspike(un).labels; % paste label into RESPONSE structure
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% GENERATE SPIKE DENSITY FUNCTIONS %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            disp('......generating spike density functions and graphstruct...')
            %%% paste individual trial info into SPIKESTRUCT structure
            pointer=find(ismember(tempbehav(:,2),parnumlist)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),hmiconfig.faces500)==1); % select only correct 500 series trials - faces
            spikestructsingle.faces_spk=tempspike(un).spikes(pointer,:);
            spikestructsingle.faces_ts=ceil(tempbehav(pointer,7)*1000); % round cue onset timestamps to nearest ms
            graphstructsingle.faces_rast=prep_raster(spikestructsingle.faces_spk,spikestructsingle.faces_ts,xscale);
            pointer=find(ismember(tempbehav(:,2),parnumlist)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),hmiconfig.fruit500)==1); % select only correct 500 series trials -
            spikestructsingle.fruit_spk=tempspike(un).spikes(pointer,:);
            spikestructsingle.fruit_ts=ceil(tempbehav(pointer,7)*1000);
            graphstructsingle.fruit_rast=prep_raster(spikestructsingle.fruit_spk,spikestructsingle.fruit_ts,xscale);
            pointer=find(ismember(tempbehav(:,2),parnumlist)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),hmiconfig.bodyp500)==1); % select only correct 500 series trials
            spikestructsingle.bodyp_spk=tempspike(un).spikes(pointer,:);
            spikestructsingle.bodyp_ts=ceil(tempbehav(pointer,7)*1000);
            graphstructsingle.bodyp_rast=prep_raster(spikestructsingle.bodyp_spk,spikestructsingle.bodyp_ts,xscale);
            pointer=find(ismember(tempbehav(:,2),parnumlist)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),hmiconfig.places500)==1); % select only correct 500 series trials
            spikestructsingle.places_spk=tempspike(un).spikes(pointer,:);
            spikestructsingle.places_ts=ceil(tempbehav(pointer,7)*1000);
            graphstructsingle.places_rast=prep_raster(spikestructsingle.places_spk,spikestructsingle.places_ts,xscale);
            pointer=find(ismember(tempbehav(:,2),parnumlist)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),hmiconfig.objct500)==1); % select only correct 500 series trials
            spikestructsingle.objct_spk=tempspike(un).spikes(pointer,:);
            spikestructsingle.objct_ts=ceil(tempbehav(pointer,7)*1000);
            graphstructsingle.objct_rast=prep_raster(spikestructsingle.objct_spk,spikestructsingle.objct_ts,xscale);

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Prepare GRAPHSTRUCT structure (used for storing average spden functions)
            graphstructsingle.allconds=unique(tempbehav(ismember(tempbehav(:,2),parnumlist)==1,3));
            graphstructsingle.allconds_avg=zeros(length(unique(graphstructsingle.allconds)),5000);
            graphstructsingle.allconds_sem=zeros(length(unique(graphstructsingle.allconds)),5000);
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Generate average spike density functions for each CONDITION (listed in CAT_avg/CAT_sem)
            [graphstructsingle.faces_avg,graphstructsingle.faces_sem]=plx_avgspden(spikestructsingle.faces_spk,spikestructsingle.faces_ts,5000,1,hmiconfig.gausskernel);
            [graphstructsingle.fruit_avg,graphstructsingle.fruit_sem]=plx_avgspden(spikestructsingle.fruit_spk,spikestructsingle.fruit_ts,5000,1,hmiconfig.gausskernel);
            [graphstructsingle.bodyp_avg,graphstructsingle.bodyp_sem]=plx_avgspden(spikestructsingle.bodyp_spk,spikestructsingle.bodyp_ts,5000,1,hmiconfig.gausskernel);
            [graphstructsingle.places_avg,graphstructsingle.places_sem]=plx_avgspden(spikestructsingle.places_spk,spikestructsingle.places_ts,5000,1,hmiconfig.gausskernel);
            [graphstructsingle.objct_avg,graphstructsingle.objct_sem]=plx_avgspden(spikestructsingle.objct_spk,spikestructsingle.objct_ts,5000,1,hmiconfig.gausskernel);
            graphstructsingle.allconds_avg=zeros(100,5000)*.1;

            %%%%%%%%%%%%%%%%%%
            %%% Trial loop %%%
            %%%%%%%%%%%%%%%%%%
            disp('......analyzing each TRIAL...')
            for tr=1:size(tempspike(un).spikes,1),
                %%% Generate spike density function for each TRIAL
                tempspikes=tempspike(un).spikes(tr,:);
                temp_ts=ceil(tempbehav(tr,7)*1000); % round cue onset timestamps to nearest ms
                [graphstructsingle.spden_trial(tr,:),junk]=plx_avgspden(tempspikes,temp_ts,5000,1,hmiconfig.gausskernel);
                %respstructsingle.trial_spk_baseline(tr)=length(find(tempspikes>temp_ts+baseline(1)&tempspikes<temp_ts+baseline(2)));
                %respstructsingle.trial_spk_epoch1(tr)=length(find(tempspikes>temp_ts+epoch1(1)&tempspikes<temp_ts+epoch1(2)));
                %respstructsingle.trial_spk_epoch2(tr)=length(find(tempspikes>temp_ts+epoch2(1)&tempspikes<temp_ts+epoch2(2)));
                %respstructsingle.trial_spk_epoch3(tr)=length(find(tempspikes>temp_ts+epoch3(1)&tempspikes<temp_ts+epoch3(2)));
                respstructsingle.trial_m_baseline(tr)=mean(graphstructsingle.spden_trial(tr,baseline(1)+1000:baseline(2)+1000)');
                respstructsingle.trial_m_epoch1(tr)=mean(graphstructsingle.spden_trial(tr,epoch1(1)+1000:epoch1(2)+1000)');
                respstructsingle.trial_m_epoch2(tr)=mean(graphstructsingle.spden_trial(tr,epoch2(1)+1000:epoch2(2)+1000)');
                respstructsingle.trial_m_epoch3(tr)=mean(graphstructsingle.spden_trial(tr,epoch3(1)+1000:epoch3(2)+1000)');
                %[respstructsingle.trial_p_baseline(tr),ind]=max(graphstructsingle.spden_trial(tr,baseline(1)+1000:baseline(2)+1000)');
                %[respstructsingle.trial_p_epoch1(tr),ind]=max(graphstructsingle.spden_trial(tr,epoch1(1)+1000:epoch1(2)+1000)');
                %[respstructsingle.trial_p_epoch2(tr)]=max(graphstructsingle.spden_trial(tr,epoch2(1)+1000:epoch2(2)+1000)');
                %[respstructsingle.trial_p_epoch3(tr),ind]=max(graphstructsingle.spden_trial(tr,epoch3(1)+1000:epoch3(2)+1000)');
                %respstructsingle.trial_mp_baseline(tr)=mean(graphstructsingle.spden_trial(tr,1000+ind+baseline(1)+mpwin(1):1000+ind+baseline(2)+mpwin(2))');
                %respstructsingle.trial_mp_epoch1(tr)=mean(graphstructsingle.spden_trial(tr,1000+ind+epoch1(1)+mpwin(1):1000+ind+epoch1(1)+mpwin(2))');
                %respstructsingle.trial_mp_epoch2(tr)=mean(graphstructsingle.spden_trial(tr,1000+ind+epoch2(1)+mpwin(1):1000+ind+epoch2(1)+mpwin(2))');
                %respstructsingle.trial_mp_epoch3(tr)=mean(graphstructsingle.spden_trial(tr,1000+ind+epoch3(1)+mpwin(1):1000+ind+epoch3(1)+mpwin(2))');
                %respstructsingle.trial_area_baseline(tr)=sum(graphstructsingle.spden_trial(tr,1000+area_baseline(1):1000+area_baseline(2))');
                %respstructsingle.trial_area_epoch1(tr)=sum(graphstructsingle.spden_trial(tr,1000+area_epoch1(1):1000+area_epoch1(2))');
                respstructsingle.trial_id(tr,1)=tempbehav(tr,3); % paste stimulus number
                switch tempbehav(tr,3) %% assign category number
                    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}, respstructsingle.trial_id(tr,2)=1;
                    case {21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40}, respstructsingle.trial_id(tr,2)=2;
                    case {41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60}, respstructsingle.trial_id(tr,2)=3;
                    case {61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80}, respstructsingle.trial_id(tr,2)=4;
                    case {81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100}, respstructsingle.trial_id(tr,2)=5;
                end
            end
            %%%%%%%%%%%%%%%%
            %%% ANALYSIS %%%
            %%%%%%%%%%%%%%%%
            %%% Condition loop
            disp('......analyzing each CONDITION...')
            for cnd=1:100, % first loop creates average spike density function for each condition
                %%% Generate average spike density functions for each STIMULUS (graphstruct.allconds_avg(cnd,:) and graphstruct.allconds_sem(cnd,:))
                pointer=find(ismember(tempbehav(:,2),parnumlist)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),cnd)==1); % locate all correct trials that match condition number
                if isempty(pointer)==1, % if no trials are found, paste zero values
                    graphstructsingle.allconds_avg(cnd,:)=zeros(1,5000);
                    graphstructsingle.allconds_sem(cnd,:)=zeros(1,5000);
                else
                    tempspikes=tempspike(un).spikes(pointer,:);
                    temp_ts=ceil(tempbehav(pointer,7)*1000); % round cue onset timestamps to nearest ms
                    [graphstructsingle.allconds_avg(cnd,:),graphstructsingle.allconds_sem(cnd,:)]=plx_avgspden(tempspikes,temp_ts,5000,1,hmiconfig.gausskernel);
                    if isnan(graphstructsingle.allconds_avg(cnd,1))==1, % fills in empty rows
                        graphstructsingle.allconds_avg(cnd,:)=zeros(1,5000);
                    end
                end
            end
            for cnd=1:100, % first loop creates average spike density function for each condition
                pointer=find(ismember(tempbehav(:,2),parnumlist)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),cnd)==1); % select correct trials matching cnd#
                tempspikes=tempspike(un).spikes(pointer,:);
                temp_ts=ceil(tempbehav(pointer,7)*1000); % round cue onset timestamps to nearest ms
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                %%% Paste spike counts, etc into RESPSTRUCT structure
                tempbaseline=0; temp1=0; temp2=0; temp3=0;
                for tr=1:size(tempspikes,1), % scroll through each trial for each condition
                    tempbaseline=tempbaseline+(length(find(tempspikes(tr,:)>=temp_ts(tr)+baseline(1) & tempspikes(tr,:)<=temp_ts(tr)+baseline(2))));
                    temp1=temp1+(length(find(tempspikes(tr,:)>=temp_ts(tr)+epoch1(1) & tempspikes(tr,:)<=temp_ts(tr)+epoch1(2))));
                    temp2=temp2+(length(find(tempspikes(tr,:)>=temp_ts(tr)+epoch2(1) & tempspikes(tr,:)<=temp_ts(tr)+epoch2(2))));
                    temp3=temp3+(length(find(tempspikes(tr,:)>=temp_ts(tr)+epoch3(1) & tempspikes(tr,:)<=temp_ts(tr)+epoch3(2))));
                end
                %%% calculate mean baseline measures
                %respstructsingle.spk_baseline(cnd)=tempbaseline/size(tempspikes,1); % normalize to spikes/trial
                respstructsingle.m_baseline(cnd)=mean(mean(graphstructsingle.allconds_avg(:,baseline(1)+1000:baseline(2)+1000)')); % average baseline rate
                %[respstructsingle.p_baseline(cnd),ind1]=max(graphstructsingle.allconds_avg(cnd,baseline(1)+1000:baseline(2)+1000)');
                %respstructsingle.p_baseline(cnd)=mean(mean(graphstructsingle.allconds_avg(:,baseline(1)+1000:baseline(2)+1000)')); % average baseline rate
                %respstructsingle.mp_baseline(cnd)=mean(graphstructsingle.allconds_avg(cnd,baseline(1)+1000+ind1+mpwin(1):baseline(2)+1000+ind1+mpwin(2))');
                %respstructsingle.mp_baseline(cnd)=mean(mean(graphstructsingle.allconds_avg(:,baseline(1)+1000:baseline(2)+1000)')); % average baseline rate
                %respstructsingle.spk_epoch1(cnd)=temp1/size(tempspikes,1); % normalize to spikes/trial
                %respstructsingle.spk_epoch2(cnd)=temp2/size(tempspikes,1); % normalize to spikes/trial
                %respstructsingle.spk_epoch3(cnd)=temp3/size(tempspikes,1); % normalize to spikes/trial
                respstructsingle.m_epoch1(cnd)=mean(graphstructsingle.allconds_avg(cnd,epoch1(1)+1000:epoch1(2)+1000)');
                respstructsingle.m_epoch2(cnd)=mean(graphstructsingle.allconds_avg(cnd,epoch2(1)+1000:epoch2(2)+1000)');
                respstructsingle.m_epoch3(cnd)=mean(graphstructsingle.allconds_avg(cnd,epoch3(1)+1000:epoch3(2)+1000)');
                %[respstructsingle.p_epoch1(cnd),ind2]=max(graphstructsingle.allconds_avg(cnd,epoch1(1)+1000:epoch1(2)+1000)');
                %[respstructsingle.p_epoch2(cnd),ind3]=max(graphstructsingle.allconds_avg(cnd,epoch2(1)+1000:epoch2(2)+1000)');
                %[respstructsingle.p_epoch3(cnd),ind4]=max(graphstructsingle.allconds_avg(cnd,epoch3(1)+1000:epoch3(2)+1000)');
                %respstructsingle.mp_epoch1(cnd)=mean(graphstructsingle.allconds_avg(cnd,epoch1(1)+1000+ind2+mpwin(1):epoch1(2)+1000+ind2+mpwin(2))');
                %respstructsingle.mp_epoch2(cnd)=mean(graphstructsingle.allconds_avg(cnd,epoch2(1)+1000+ind3+mpwin(1):epoch2(2)+1000+ind3+mpwin(2))');
                %respstructsingle.mp_epoch3(cnd)=mean(graphstructsingle.allconds_avg(cnd,epoch3(1)+1000+ind4+mpwin(1):epoch3(2)+1000+ind4+mpwin(2))');
                %respstructsingle.area_baseline(cnd)=sum(graphstructsingle.allconds_avg(cnd,1000+area_baseline(1):1000+area_baseline(2))');
                %respstructsingle.area_epoch1(cnd)=sum(graphstructsingle.allconds_avg(cnd,1000+area_epoch1(1):1000+area_epoch1(2))');
                %%% calculate condition latency
                respstructsingle.latency(cnd)=plx_calclatency(graphstructsingle.allconds_avg(cnd,:),...
                    mean(graphstructsingle.allconds_sem(cnd,:)),1000,respstructsingle.m_baseline(cnd),size(tempspikes,1));
                if respstructsingle.latency(cnd)<minlatency, respstructsingle.latency(cnd)=0; end % remove any latencies less than 50ms
                %%% Determine validity of individual responses (currently uses only epoch1 to classify individual stimulus response)
                %if respstructsingle.spk_epoch1(cnd) < (2*respstructsingle.spk_baseline(cnd)),
                %    respstructsingle.validrsp(cnd,1)=0; % classifies as valid/invalid (>2X baseline)
                %else respstructsingle.validrsp(cnd,1)=1; end
                if respstructsingle.m_epoch1(cnd) < (2*respstructsingle.m_baseline(cnd)),
                    respstructsingle.validrsp(cnd,2)=0; % classifies as valid/invalid (>2X baseline)
                else respstructsingle.validrsp(cnd,2)=1; end
                %if respstructsingle.p_epoch1(cnd) < (2*respstructsingle.p_baseline(cnd)),
                %    respstructsingle.validrsp(cnd,3)=0; % classifies as valid/invalid (>2X baseline)
                %else respstructsingle.validrsp(cnd,3)=1; end
                %if respstructsingle.mp_epoch1(cnd) < (2*respstructsingle.mp_baseline(cnd)),
                %    respstructsingle.validrsp(cnd,4)=0; % classifies as valid/invalid (>2X baseline)
                %else respstructsingle.validrsp(cnd,4)=1; end
                %if respstructsingle.area_epoch1(cnd) < (2*respstructsingle.area_baseline(cnd)),
                %    respstructsingle.validrsp(cnd,5)=0; % classifies as valid/invalid (>2X baseline)
                %else respstructsingle.validrsp(cnd,5)=1; end
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Subtract Baseline and repeat epoch quantification
            %spk_avgbaseline=mean(respstructsingle.spk_baseline);
            m_avgbaseline=mean(respstructsingle.m_baseline);
            %p_avgbaseline=mean(respstructsingle.p_baseline);
            %mp_avgbaseline=mean(respstructsingle.mp_baseline);
            %area_avgbaseline=mean(respstructsingle.area_baseline);
            for cnd=1:100, % third condition loop subtracts average baseline and recalculates analysis parameters
                %respstructsingle.spk_epoch1_nobase(cnd)=respstructsingle.spk_epoch1(cnd)-spk_avgbaseline;
                %respstructsingle.spk_epoch2_nobase(cnd)=respstructsingle.spk_epoch2(cnd)-spk_avgbaseline;
                %respstructsingle.spk_epoch3_nobase(cnd)=respstructsingle.spk_epoch3(cnd)-spk_avgbaseline;
                respstructsingle.m_epoch1_nobase(cnd)=respstructsingle.m_epoch1(cnd)-m_avgbaseline;
                respstructsingle.m_epoch2_nobase(cnd)=respstructsingle.m_epoch2(cnd)-m_avgbaseline;
                respstructsingle.m_epoch3_nobase(cnd)=respstructsingle.m_epoch3(cnd)-m_avgbaseline;
                %respstructsingle.p_epoch1_nobase(cnd)=respstructsingle.p_epoch1(cnd)-p_avgbaseline;
                %respstructsingle.p_epoch2_nobase(cnd)=respstructsingle.p_epoch2(cnd)-p_avgbaseline;
                %respstructsingle.p_epoch3_nobase(cnd)=respstructsingle.p_epoch3(cnd)-p_avgbaseline;
                %respstructsingle.mp_epoch1_nobase(cnd)=respstructsingle.mp_epoch1(cnd)-mp_avgbaseline;
                %respstructsingle.mp_epoch2_nobase(cnd)=respstructsingle.mp_epoch2(cnd)-mp_avgbaseline;
                %respstructsingle.mp_epoch3_nobase(cnd)=respstructsingle.mp_epoch3(cnd)-mp_avgbaseline;
                %respstructsingle.area_epoch1_nobase(cnd)=respstructsingle.area_epoch1(cnd)-area_avgbaseline;
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Determine average/best categorical responses (1 col/metric)
            [respstructsingle.cat_avg(1,:),respstructsingle.cat_sem(1,:),respstructsingle.cat_bst(1,:)]=rsp_avgbest(respstructsingle,hmiconfig.faces500);
            [respstructsingle.cat_avg(2,:),respstructsingle.cat_sem(2,:),respstructsingle.cat_bst(2,:)]=rsp_avgbest(respstructsingle,hmiconfig.fruit500);
            [respstructsingle.cat_avg(3,:),respstructsingle.cat_sem(3,:),respstructsingle.cat_bst(3,:)]=rsp_avgbest(respstructsingle,hmiconfig.places500);
            [respstructsingle.cat_avg(4,:),respstructsingle.cat_sem(4,:),respstructsingle.cat_bst(4,:)]=rsp_avgbest(respstructsingle,hmiconfig.bodyp500);
            [respstructsingle.cat_avg(5,:),respstructsingle.cat_sem(5,:),respstructsingle.cat_bst(5,:)]=rsp_avgbest(respstructsingle,hmiconfig.objct500);
            %%% No baseline
            [respstructsingle.cat_avg_nobase(1,:),respstructsingle.cat_sem_nobase(1,:),respstructsingle.cat_bst_nobase(1,:)]=rsp_avgbest_nobase(respstructsingle,hmiconfig.faces500);
            [respstructsingle.cat_avg_nobase(2,:),respstructsingle.cat_sem_nobase(2,:),respstructsingle.cat_bst_nobase(2,:)]=rsp_avgbest_nobase(respstructsingle,hmiconfig.fruit500);
            [respstructsingle.cat_avg_nobase(3,:),respstructsingle.cat_sem_nobase(3,:),respstructsingle.cat_bst_nobase(3,:)]=rsp_avgbest_nobase(respstructsingle,hmiconfig.places500);
            [respstructsingle.cat_avg_nobase(4,:),respstructsingle.cat_sem_nobase(4,:),respstructsingle.cat_bst_nobase(4,:)]=rsp_avgbest_nobase(respstructsingle,hmiconfig.bodyp500);
            [respstructsingle.cat_avg_nobase(5,:),respstructsingle.cat_sem_nobase(5,:),respstructsingle.cat_bst_nobase(5,:)]=rsp_avgbest_nobase(respstructsingle,hmiconfig.objct500);
            %%% Calculate mean categorical latencies
            [respstructsingle.cat_latency(1,1),respstructsingle.cat_latency(1,2)]=mean_sem(nonzeros(respstructsingle.latency(hmiconfig.faces500)));
            [respstructsingle.cat_latency(2,1),respstructsingle.cat_latency(2,2)]=mean_sem(nonzeros(respstructsingle.latency(hmiconfig.fruit500)));
            [respstructsingle.cat_latency(3,1),respstructsingle.cat_latency(3,2)]=mean_sem(nonzeros(respstructsingle.latency(hmiconfig.places500)));
            [respstructsingle.cat_latency(4,1),respstructsingle.cat_latency(4,2)]=mean_sem(nonzeros(respstructsingle.latency(hmiconfig.bodyp500)));
            [respstructsingle.cat_latency(5,1),respstructsingle.cat_latency(5,2)]=mean_sem(nonzeros(respstructsingle.latency(hmiconfig.objct500)));
            %%% Sensory/Non-sensory for each category (compares the mean baseline to the mean response in epoch1)
            %respstructsingle.cat_sensory(1,1)=signrank(respstructsingle.spk_baseline(hmiconfig.faces500),respstructsingle.spk_epoch1(hmiconfig.faces500));
            %respstructsingle.cat_sensory(2,1)=signrank(respstructsingle.spk_baseline(hmiconfig.fruit500),respstructsingle.spk_epoch1(hmiconfig.fruit500));
            %respstructsingle.cat_sensory(3,1)=signrank(respstructsingle.spk_baseline(hmiconfig.places500),respstructsingle.spk_epoch1(hmiconfig.places500));
            %respstructsingle.cat_sensory(4,1)=signrank(respstructsingle.spk_baseline(hmiconfig.bodyp500),respstructsingle.spk_epoch1(hmiconfig.bodyp500));
            %respstructsingle.cat_sensory(5,1)=signrank(respstructsingle.spk_baseline(hmiconfig.objct500),respstructsingle.spk_epoch1(hmiconfig.objct500));
            respstructsingle.cat_sensory(1,2)=signrank(respstructsingle.m_baseline(hmiconfig.faces500),respstructsingle.m_epoch1(hmiconfig.faces500));
            respstructsingle.cat_sensory(2,2)=signrank(respstructsingle.m_baseline(hmiconfig.fruit500),respstructsingle.m_epoch1(hmiconfig.fruit500));
            respstructsingle.cat_sensory(3,2)=signrank(respstructsingle.m_baseline(hmiconfig.places500),respstructsingle.m_epoch1(hmiconfig.places500));
            respstructsingle.cat_sensory(4,2)=signrank(respstructsingle.m_baseline(hmiconfig.bodyp500),respstructsingle.m_epoch1(hmiconfig.bodyp500));
            respstructsingle.cat_sensory(5,2)=signrank(respstructsingle.m_baseline(hmiconfig.objct500),respstructsingle.m_epoch1(hmiconfig.objct500));
            %respstructsingle.cat_sensory(1,3)=signrank(respstructsingle.p_baseline(hmiconfig.faces500),respstructsingle.p_epoch1(hmiconfig.faces500));
            %respstructsingle.cat_sensory(2,3)=signrank(respstructsingle.p_baseline(hmiconfig.fruit500),respstructsingle.p_epoch1(hmiconfig.fruit500));
            %respstructsingle.cat_sensory(3,3)=signrank(respstructsingle.p_baseline(hmiconfig.places500),respstructsingle.p_epoch1(hmiconfig.places500));
            %respstructsingle.cat_sensory(4,3)=signrank(respstructsingle.p_baseline(hmiconfig.bodyp500),respstructsingle.p_epoch1(hmiconfig.bodyp500));
            %respstructsingle.cat_sensory(5,3)=signrank(respstructsingle.p_baseline(hmiconfig.objct500),respstructsingle.p_epoch1(hmiconfig.objct500));
            %respstructsingle.cat_sensory(1,4)=signrank(respstructsingle.mp_baseline(hmiconfig.faces500),respstructsingle.mp_epoch1(hmiconfig.faces500));
            %respstructsingle.cat_sensory(2,4)=signrank(respstructsingle.mp_baseline(hmiconfig.fruit500),respstructsingle.mp_epoch1(hmiconfig.fruit500));
            %respstructsingle.cat_sensory(3,4)=signrank(respstructsingle.mp_baseline(hmiconfig.places500),respstructsingle.mp_epoch1(hmiconfig.places500));
            %respstructsingle.cat_sensory(4,4)=signrank(respstructsingle.mp_baseline(hmiconfig.bodyp500),respstructsingle.mp_epoch1(hmiconfig.bodyp500));
            %respstructsingle.cat_sensory(5,4)=signrank(respstructsingle.mp_baseline(hmiconfig.objct500),respstructsingle.mp_epoch1(hmiconfig.objct500));
            %respstructsingle.cat_sensory(1,5)=signrank(respstructsingle.area_baseline(hmiconfig.faces500),respstructsingle.area_epoch1(hmiconfig.faces500));
            %respstructsingle.cat_sensory(2,5)=signrank(respstructsingle.area_baseline(hmiconfig.fruit500),respstructsingle.area_epoch1(hmiconfig.fruit500));
            %respstructsingle.cat_sensory(3,5)=signrank(respstructsingle.area_baseline(hmiconfig.places500),respstructsingle.area_epoch1(hmiconfig.places500));
            %respstructsingle.cat_sensory(4,5)=signrank(respstructsingle.area_baseline(hmiconfig.bodyp500),respstructsingle.area_epoch1(hmiconfig.bodyp500));
            %respstructsingle.cat_sensory(5,5)=signrank(respstructsingle.area_baseline(hmiconfig.objct500),respstructsingle.area_epoch1(hmiconfig.objct500));

            %%% Determine EXCITE/INHIBIT/BOTH
            % algorithm based on significant difference between baseline
            % and epoch1.
            respstructsingle.excite_inhibit=zeros(1,5);
            if respstructsingle.cat_sensory(1,metric_col)<0.06, % faces
                if mean(respstructsingle.m_baseline(hmiconfig.faces500))<mean(respstructsingle.m_epoch1(hmiconfig.faces500)),
                    respstructsingle.excite_inhibit(1)=1; else respstructsingle.excite_inhibit(1)=-1;
                end
            end
            if respstructsingle.cat_sensory(2,metric_col)<0.06, % fruit
                if mean(respstructsingle.m_baseline(hmiconfig.fruit500))<mean(respstructsingle.m_epoch1(hmiconfig.fruit500)),
                    respstructsingle.excite_inhibit(2)=1; else respstructsingle.excite_inhibit(2)=-1;
                end
            end
            if respstructsingle.cat_sensory(3,metric_col)<0.06, % places
                if mean(respstructsingle.m_baseline(hmiconfig.places500))<mean(respstructsingle.m_epoch1(hmiconfig.places500)),
                    respstructsingle.excite_inhibit(3)=1; else respstructsingle.excite_inhibit(3)=-1;
                end
            end
            if respstructsingle.cat_sensory(4,metric_col)<0.06, % bodyparts
                if mean(respstructsingle.m_baseline(hmiconfig.bodyp500))<mean(respstructsingle.m_epoch1(hmiconfig.bodyp500)),
                    respstructsingle.excite_inhibit(4)=1; else respstructsingle.excite_inhibit(4)=-1;
                end
            end
            if respstructsingle.cat_sensory(5,metric_col)<0.06, % objects
                if mean(respstructsingle.m_baseline(hmiconfig.objct500))<mean(respstructsingle.m_epoch1(hmiconfig.objct500)),
                    respstructsingle.excite_inhibit(5)=1; else respstructsingle.excite_inhibit(5)=-1;
                end
            end
            excitemarkers=find(respstructsingle.excite_inhibit==1);
            inhibitmarkers=find(respstructsingle.excite_inhibit==-1);
            if isempty(excitemarkers)~=1 & isempty(inhibitmarkers)~=1, respstructsingle.excitetype='Both';
            elseif isempty(excitemarkers)~=1 & isempty(inhibitmarkers)==1, respstructsingle.excitetype='Excite';
            elseif isempty(excitemarkers)==1 & isempty(inhibitmarkers)~=1, respstructsingle.excitetype='Inhibit';
            elseif isempty(excitemarkers)==1 & isempty(inhibitmarkers)==1, respstructsingle.excitetype='Non-Responsive';
            end
            catnames={'Faces','Fruit','Places','BodyParts','Objects'};
            if strcmp(respstructsingle.excitetype,'Excite')==1
                [m1,m2]=max(respstructsingle.cat_avg(:,2));
                respstructsingle.pref_excite=catnames(m2);
                respstructsingle.pref_inhibit='None';
            elseif strcmp(respstructsingle.excitetype,'Inhibit')==1,
                [m1,m2]=min(respstructsingle.cat_avg(:,2));
                respstructsingle.pref_excite='None';
                respstructsingle.pref_inhibit=catnames(m2);
            elseif strcmp(respstructsingle.excitetype,'Both')==1,
                [m1,m2]=max(respstructsingle.cat_avg(:,2));
                respstructsingle.pref_excite=catnames(m2);
                [m1,m2]=min(respstructsingle.cat_avg(:,2));
                respstructsingle.pref_inhibit=catnames(m2);
            elseif strcmp(respstructsingle.excitetype,'Non-Responsive')==1,
                respstructsingle.pref_excite='None';
                respstructsingle.pref_inhibit='None';
            end

            %%% adding label to preferred category
            [junk,ind]=max(respstructsingle.cat_avg(:,metric_col));
            switch ind
                case 1, respstructsingle.preferred_category='Faces'; respstructsingle.preferred_sensory=respstructsingle.cat_sensory(1,metric_col);
                case 2, respstructsingle.preferred_category='Fruit'; respstructsingle.preferred_sensory=respstructsingle.cat_sensory(2,metric_col);
                case 3, respstructsingle.preferred_category='Places'; respstructsingle.preferred_sensory=respstructsingle.cat_sensory(3,metric_col);
                case 4, respstructsingle.preferred_category='BodyParts'; respstructsingle.preferred_sensory=respstructsingle.cat_sensory(4,metric_col);
                case 5, respstructsingle.preferred_category='Objects'; respstructsingle.preferred_sensory=respstructsingle.cat_sensory(5,metric_col);
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Determine best and worst stimulus (based on peak responses over cue response window (defined above)
            tempmax=max(graphstructsingle.allconds_avg(:,epoch1(1)+1000:epoch1(2)+1000)'); % calculate peak response over cue response
            tempmin=min(graphstructsingle.allconds_avg(:,epoch1(1)+1000:epoch1(2)+1000)'); % calculate peak response over cue response
            catind=1:20:101;
            for category=1:length(catind)-1,
                [crap,graphstructsingle.bestconds(category)]=max(tempmax(catind(category):catind(category+1)-1));
                [crap,graphstructsingle.worstconds(category)]=min(tempmin(catind(category):catind(category+1)-1));
            end
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Calculate selectivity indices (here is the place to make
            %%% changes to the algorithms !!!
            % each column refers to category
            % each row refers to metric
            %respstructsingle.cat_si(1,1)=calc_si(respstructsingle.cat_avg,1,1); % spk, faces
            %respstructsingle.cat_si(2,1)=calc_si(respstructsingle.cat_avg,2,1); % spk, fruit
            %respstructsingle.cat_si(3,1)=calc_si(respstructsingle.cat_avg,3,1); % spk, places
            %respstructsingle.cat_si(4,1)=calc_si(respstructsingle.cat_avg,4,1); % spk, bodyparts
            %respstructsingle.cat_si(5,1)=calc_si(respstructsingle.cat_avg,5,1); % spk, objects
            respstructsingle.cat_si(1,2)=calc_si(respstructsingle.cat_avg,1,2); % mean, faces
            respstructsingle.cat_si(2,2)=calc_si(respstructsingle.cat_avg,2,2); % mean, fruit
            respstructsingle.cat_si(3,2)=calc_si(respstructsingle.cat_avg,3,2); % mean, places
            respstructsingle.cat_si(4,2)=calc_si(respstructsingle.cat_avg,4,2); % mean, bodyparts
            respstructsingle.cat_si(5,2)=calc_si(respstructsingle.cat_avg,5,2); % mean, objects
            %respstructsingle.cat_si(1,3)=calc_si(respstructsingle.cat_avg,1,3); % peak, faces
            %respstructsingle.cat_si(2,3)=calc_si(respstructsingle.cat_avg,2,3); % peak, fruit
            %respstructsingle.cat_si(3,3)=calc_si(respstructsingle.cat_avg,3,3); % peak, places
            %respstructsingle.cat_si(4,3)=calc_si(respstructsingle.cat_avg,4,3); % peak, bodyparts
            %respstructsingle.cat_si(5,3)=calc_si(respstructsingle.cat_avg,5,3); % peak, objects
            %respstructsingle.cat_si(1,4)=calc_si(respstructsingle.cat_avg,1,4); % meanpeak, faces
            %respstructsingle.cat_si(2,4)=calc_si(respstructsingle.cat_avg,2,4); % meanpeak, fruit
            %respstructsingle.cat_si(3,4)=calc_si(respstructsingle.cat_avg,3,4); % meanpeak, places
            %respstructsingle.cat_si(4,4)=calc_si(respstructsingle.cat_avg,4,4); % meanpeak, bodyparts
            %respstructsingle.cat_si(5,4)=calc_si(respstructsingle.cat_avg,5,4); % meanpeak, objects
            %respstructsingle.cat_si(1,5)=calc_si(respstructsingle.cat_avg,1,5); % area, faces
            %respstructsingle.cat_si(2,5)=calc_si(respstructsingle.cat_avg,2,5); % area, fruit
            %respstructsingle.cat_si(3,5)=calc_si(respstructsingle.cat_avg,3,5); % area, places
            %respstructsingle.cat_si(4,5)=calc_si(respstructsingle.cat_avg,4,5); % area, bodyparts
            %respstructsingle.cat_si(5,5)=calc_si(respstructsingle.cat_avg,5,5); % area, objects
            %%% face selectivity without FRUIT
            %respstructsingle.cat_si(6,1)=calc_si_nofruit(respstructsingle.cat_avg,1,1);
            respstructsingle.cat_si(6,2)=calc_si_nofruit(respstructsingle.cat_avg,1,2);
            %respstructsingle.cat_si(6,3)=calc_si_nofruit(respstructsingle.cat_avg,1,3);
            %respstructsingle.cat_si(6,4)=calc_si_nofruit(respstructsingle.cat_avg,1,4);
            %respstructsingle.cat_si(6,5)=calc_si_nofruit(respstructsingle.cat_avg,1,5);
            %%% calculate non-specific selectivity
            %respstructsingle.raw_si(1)=calc_rawsi(respstructsingle.cat_avg,1);
            respstructsingle.raw_si(2)=calc_rawsi(respstructsingle.cat_avg,2);
            %respstructsingle.raw_si(3)=calc_rawsi(respstructsingle.cat_avg,3);
            %respstructsingle.raw_si(4)=calc_rawsi(respstructsingle.cat_avg,4);
            %respstructsingle.raw_si(5)=calc_rawsi(respstructsingle.cat_avg,5);
            nonface=mean(respstructsingle.cat_avg(2:4,metric_col));
            %%% Calculate selectivity indices AFTER subtracting baseline
            %respstructsingle.cat_si_nobase(1,1)=calc_si(respstructsingle.cat_avg_nobase,1,1);
            %respstructsingle.cat_si_nobase(2,1)=calc_si(respstructsingle.cat_avg_nobase,2,1);
            %respstructsingle.cat_si_nobase(3,1)=calc_si(respstructsingle.cat_avg_nobase,3,1);
            %respstructsingle.cat_si_nobase(4,1)=calc_si(respstructsingle.cat_avg_nobase,4,1);
            %respstructsingle.cat_si_nobase(5,1)=calc_si(respstructsingle.cat_avg_nobase,5,1);
            respstructsingle.cat_si_nobase(1,2)=calc_si(respstructsingle.cat_avg_nobase,1,2);
            respstructsingle.cat_si_nobase(2,2)=calc_si(respstructsingle.cat_avg_nobase,2,2);
            respstructsingle.cat_si_nobase(3,2)=calc_si(respstructsingle.cat_avg_nobase,3,2);
            respstructsingle.cat_si_nobase(4,2)=calc_si(respstructsingle.cat_avg_nobase,4,2);
            respstructsingle.cat_si_nobase(5,2)=calc_si(respstructsingle.cat_avg_nobase,5,2);
            %respstructsingle.cat_si_nobase(1,3)=calc_si(respstructsingle.cat_avg_nobase,1,3);
            %respstructsingle.cat_si_nobase(2,3)=calc_si(respstructsingle.cat_avg_nobase,2,3);
            %respstructsingle.cat_si_nobase(3,3)=calc_si(respstructsingle.cat_avg_nobase,3,3);
            %respstructsingle.cat_si_nobase(4,3)=calc_si(respstructsingle.cat_avg_nobase,4,3);
            %respstructsingle.cat_si_nobase(5,3)=calc_si(respstructsingle.cat_avg_nobase,5,3);
            %respstructsingle.cat_si_nobase(1,4)=calc_si(respstructsingle.cat_avg_nobase,1,4);
            %respstructsingle.cat_si_nobase(2,4)=calc_si(respstructsingle.cat_avg_nobase,2,4);
            %respstructsingle.cat_si_nobase(3,4)=calc_si(respstructsingle.cat_avg_nobase,3,4);
            %respstructsingle.cat_si_nobase(4,4)=calc_si(respstructsingle.cat_avg_nobase,4,4);
            %respstructsingle.cat_si_nobase(5,4)=calc_si(respstructsingle.cat_avg_nobase,5,4);
            %respstructsingle.cat_si_nobase(1,4)=calc_si(respstructsingle.cat_avg_nobase,1,5);
            %respstructsingle.cat_si_nobase(2,4)=calc_si(respstructsingle.cat_avg_nobase,2,5);
            %respstructsingle.cat_si_nobase(3,4)=calc_si(respstructsingle.cat_avg_nobase,3,5);
            %respstructsingle.cat_si_nobase(4,4)=calc_si(respstructsingle.cat_avg_nobase,4,5);
            %respstructsingle.cat_si_nobase(5,4)=calc_si(respstructsingle.cat_avg_nobase,5,5);
            %%% face selectivity without FRUIT
            %respstructsingle.cat_si_nobase(6,1)=calc_si_nofruit(respstructsingle.cat_avg_nobase,1,1);
            respstructsingle.cat_si_nobase(6,2)=calc_si_nofruit(respstructsingle.cat_avg_nobase,1,2);
            %respstructsingle.cat_si_nobase(6,3)=calc_si_nofruit(respstructsingle.cat_avg_nobase,1,3);
            %respstructsingle.cat_si_nobase(6,4)=calc_si_nofruit(respstructsingle.cat_avg_nobase,1,4);
            %respstructsingle.cat_si_nobase(6,5)=calc_si_nofruit(respstructsingle.cat_avg_nobase,1,5);
            %%% calculate non-specific selectivity
            %respstructsingle.raw_si_nobase(1)=calc_rawsi(respstructsingle.cat_avg_nobase,1);
            respstructsingle.raw_si_nobase(2)=calc_rawsi(respstructsingle.cat_avg_nobase,2);
            %respstructsingle.raw_si_nobase(3)=calc_rawsi(respstructsingle.cat_avg_nobase,3);
            %respstructsingle.raw_si_nobase(4)=calc_rawsi(respstructsingle.cat_avg_nobase,4);
            %respstructsingle.raw_si_nobase(5)=calc_rawsi(respstructsingle.cat_avg_nobase,5);

            %%% calculate pure selectivity
            respstructsingle.pure_si(1,:)=calc_puresi(respstructsingle.cat_avg,metric_col,1);
            respstructsingle.pure_si(2,:)=calc_puresi(respstructsingle.cat_avg,metric_col,2);
            respstructsingle.pure_si(3,:)=calc_puresi(respstructsingle.cat_avg,metric_col,3);
            respstructsingle.pure_si(4,:)=calc_puresi(respstructsingle.cat_avg,metric_col,4);
            respstructsingle.pure_si(5,:)=calc_puresi(respstructsingle.cat_avg,metric_col,5);

            %%% Solve for traditional face-selectivity
            nonface=mean(respstructsingle.cat_avg(2:4,metric_col));
            if respstructsingle.cat_avg(1,metric_col)>(2*nonface), respstructsingle.face_trad=1;
            else respstructsingle.face_trad=0; end
            %respstructsingle.pairwise=pairwisematrix(respstructsingle,1,metric_col);

            %%% Solve for waveform parameters
            signame=char(respstructsingle.label);
            % Following code is a correction for 32 channel data
            if signame(14)=='A',
                wavedata=load([hmiconfig.wave_raw,signame(1:19),'_raw.mat']);
            else
                wavedata=load([hmiconfig.wave_raw,signame(1:20),'_raw.mat']);
            end
            wf_data=mean(wavedata.waverawdata');
            [respstructsingle.wf_params(1) respstructsingle.wf_params(2)]=min(wf_data);
            [respstructsingle.wf_params(3) respstructsingle.wf_params(4)]=max(wf_data);
            respstructsingle.wf_params(5)=(respstructsingle.wf_params(4)-respstructsingle.wf_params(2))*25;

            %%% Conduct ROC analyses (as per Bell et al 2009)
            disp('..Performing ROC analysis...')
            respstructsingle.roc_analysis=zeros(5,5);
%             pointer1=find(respstructsingle.trial_id(:,2)==1);
%             pointer2=find(respstructsingle.trial_id(:,2)==2);
%             pointer3=find(respstructsingle.trial_id(:,2)==3);
%             pointer4=find(respstructsingle.trial_id(:,2)==4);
%             pointer5=find(respstructsingle.trial_id(:,2)==5);
%             for rr=1:5, % once per condition
%                 pointerR=find(respstructsingle.trial_id(:,2)==rr);
%                 respstructsingle.roc_analysis(rr,1)=plx500_calcROC(respstructsingle.trial_m_epoch1(pointerR)',respstructsingle.trial_m_epoch1(pointer1)');
%                 respstructsingle.roc_analysis(rr,2)=plx500_calcROC(respstructsingle.trial_m_epoch1(pointerR)',respstructsingle.trial_m_epoch1(pointer2)');
%                 respstructsingle.roc_analysis(rr,3)=plx500_calcROC(respstructsingle.trial_m_epoch1(pointerR)',respstructsingle.trial_m_epoch1(pointer3)');
%                 respstructsingle.roc_analysis(rr,4)=plx500_calcROC(respstructsingle.trial_m_epoch1(pointerR)',respstructsingle.trial_m_epoch1(pointer4)');
%                 respstructsingle.roc_analysis(rr,5)=plx500_calcROC(respstructsingle.trial_m_epoch1(pointerR)',respstructsingle.trial_m_epoch1(pointer5)');
%             end

            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Prepare and complete ANOVA analysis %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% This is the only difference between quick and analyze
            if anova==1, % allow user to shutoff anova analysis to save time/stats not available
                % prepare stimulus and category vectors
                %stimulus_id=1:100';
                category_id=[ones(1,20) ones(1,20)*2 ones(1,20)*3 ones(1,20)*4 ones(1,20)*5]';
                % anova on latencies
                respstructsingle.anova_latency=anova1(respstructsingle.latency,category_id,'off');
                % anova on category responses
                %respstructsingle.anova_baseline(1)=anova1(respstructsingle.spk_baseline,category_id,'off');
                respstructsingle.anova_baseline(2)=anova1(respstructsingle.m_baseline,category_id,'off');
                %respstructsingle.anova_baseline(3)=anova1(respstructsingle.p_baseline,category_id,'off');
                %respstructsingle.anova_baseline(4)=anova1(respstructsingle.mp_baseline,category_id,'off');
                %respstructsingle.anova_baseline(5)=anova1(respstructsingle.area_baseline,category_id,'off');
                %respstructsingle.anova_epoch(1,1)=anova1(respstructsingle.spk_epoch1,category_id,'off');
                respstructsingle.anova_epoch(1,2)=anova1(respstructsingle.m_epoch1,category_id,'off');
                %respstructsingle.anova_epoch(1,3)=anova1(respstructsingle.p_epoch1,category_id,'off');
                %respstructsingle.anova_epoch(1,4)=anova1(respstructsingle.mp_epoch1,category_id,'off');
                %respstructsingle.anova_epoch(1,5)=anova1(respstructsingle.area_epoch1,category_id,'off');
                %respstructsingle.anova_epoch(2,1)=anova1(respstructsingle.spk_epoch2,category_id,'off');
                %respstructsingle.anova_epoch(2,2)=anova1(respstructsingle.m_epoch2,category_id,'off');
                %respstructsingle.anova_epoch(2,3)=anova1(respstructsingle.p_epoch2,category_id,'off');
                %respstructsingle.anova_epoch(2,4)=anova1(respstructsingle.mp_epoch2,category_id,'off');
                %respstructsingle.anova_epoch(2,5)=anova1(respstructsingle.area_epoch1,category_id,'off');
                %respstructsingle.anova_epoch(3,1)=anova1(respstructsingle.spk_epoch3,category_id,'off');
                %respstructsingle.anova_epoch(3,2)=anova1(respstructsingle.m_epoch3,category_id,'off');
                %respstructsingle.anova_epoch(3,3)=anova1(respstructsingle.p_epoch3,category_id,'off');
                %respstructsingle.anova_epoch(3,4)=anova1(respstructsingle.mp_epoch3,category_id,'off');
                %respstructsingle.anova_epoch(3,5)=anova1(respstructsingle.area_epoch1,category_id,'off');
                %%% No baseline
                %respstructsingle.anova_epoch_nobase(1,1)=anova1(respstructsingle.spk_epoch1_nobase,category_id,'off');
                respstructsingle.anova_epoch_nobase(1,2)=anova1(respstructsingle.m_epoch1_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(1,3)=anova1(respstructsingle.p_epoch1_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(1,4)=anova1(respstructsingle.mp_epoch1_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(1,5)=anova1(respstructsingle.area_epoch1_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(2,1)=anova1(respstructsingle.spk_epoch2_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(2,2)=anova1(respstructsingle.m_epoch2_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(2,3)=anova1(respstructsingle.p_epoch2_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(2,4)=anova1(respstructsingle.mp_epoch2_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(2,5)=anova1(respstructsingle.area_epoch1_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(3,1)=anova1(respstructsingle.spk_epoch3_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(3,2)=anova1(respstructsingle.m_epoch3_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(3,3)=anova1(respstructsingle.p_epoch3_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(3,4)=anova1(respstructsingle.mp_epoch3_nobase,category_id,'off');
                %respstructsingle.anova_epoch_nobase(3,5)=anova1(respstructsingle.area_epoch1_nobase,category_id,'off');
                % within category analysis
                for ca=1:5,
                    pointer=find(respstructsingle.trial_id(:,2)==ca);
                    respstructsingle.anova_within_group(ca,1,2)=anova1(respstructsingle.trial_m_epoch1(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,1,3)=anova1(respstructsingle.trial_p_epoch1(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,1,4)=anova1(respstructsingle.trial_mp_epoch1(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,1,5)=anova1(respstructsingle.trial_area_epoch1(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,2,1)=anova1(respstructsingle.trial_m_epoch2(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,2,2)=anova1(respstructsingle.trial_p_epoch2(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,2,3)=anova1(respstructsingle.trial_mp_epoch2(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,3,1)=anova1(respstructsingle.trial_m_epoch3(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,3,2)=anova1(respstructsingle.trial_p_epoch3(pointer),respstructsingle.trial_id(pointer,1),'off');
                    %respstructsingle.anova_within_group(ca,3,2)=anova1(respstructsingle.trial_mp_epoch3(pointer),respstructsingle.trial_id(pointer,1),'off');
                end % end category loop
            end % end anova analysis loop
            
            %%% Appended Fields (Added to plx500.m March19, 2009)
            %%% RAW SI FOR BOTH EXCITATORY AND INHIBITED RESPONSES
            rsp=respstructsingle.cat_avg(:,2)';
            cols=1:5;
            % excite_rawsi
            [val ind]=max(rsp);
            othercols=find(cols~=ind);
            non_ind=mean(rsp(othercols));
            maincol=rsp(ind);
            respstructsingle.excite_rawsi=(maincol-non_ind)/(maincol+non_ind);
            % inhibit_rawsi
            [val ind]=min(rsp);
            othercols=find(cols~=ind);
            non_ind=mean(rsp(othercols));
            maincol=rsp(ind);
            respstructsingle.inhibit_rawsi=(maincol-non_ind)/(maincol+non_ind);
            %%%%%%%%%%% MARCH 12, 2009 %%%%%%%%%%
            %%% WITHOUT Fruit Measures
            % Excite/Inhibit/Both - No Fruit
            tmpexciteinhibit=respstructsingle.excite_inhibit([1 3 4 5]);
            excitemarkers=find(tmpexciteinhibit==1);
            inhibitmarkers=find(tmpexciteinhibit==-1);
            if isempty(excitemarkers)~=1 & isempty(inhibitmarkers)~=1, respstructsingle.excitetype_nofruit='Both';
            elseif isempty(excitemarkers)~=1 & isempty(inhibitmarkers)==1, respstructsingle.excitetype_nofruit='Excite';
            elseif isempty(excitemarkers)==1 & isempty(inhibitmarkers)~=1, respstructsingle.excitetype_nofruit='Inhibit';
            elseif isempty(excitemarkers)==1 & isempty(inhibitmarkers)==1, respstructsingle.excitetype_nofruit='Non-Responsive';
            end
            % CategorySelectivity (ANOVA)
            category_id=[ones(1,20) ones(1,20)*2 ones(1,20)*3 ones(1,20)*4 ones(1,20)*5]'; pointer=find(category_id~=2); % eliminate fruit
            respstructsingle.catanova_nofruit=anova1(respstructsingle.m_epoch1(pointer),category_id(pointer),'off');
            if respstructsingle.catanova_nofruit<=0.05 & strcmp(respstructsingle.excitetype_nofruit,'Non-Responsive')~=1
                respstructsingle.selective_nofruit='Selective';
            elseif respstructsingle.catanova_nofruit>0.05 & strcmp(respstructsingle.excitetype_nofruit,'Non-Responsive')~=1
                respstructsingle.selective_nofruit='Not Selective';
            else
                respstructsingle.selective_nofruit='Non-Responsive';
            end
            %%% Preferred Category - No Fruit
            [junk,ind]=max(respstructsingle.cat_avg([1 3 4 5],2));
            switch ind
                case 1, respstructsingle.prefcat_excite_nofruit='Faces';
                case 2, respstructsingle.prefcat_excite_nofruit='Places';
                case 3, respstructsingle.prefcat_excite_nofruit='BodyParts';
                case 4, respstructsingle.prefcat_excite_nofruit='Objects';
            end
            %%% Preferred Category - No Fruit
            [junk,ind]=min(respstructsingle.cat_avg([1 3 4 5],2));
            switch ind
                case 1, respstructsingle.prefcat_inhibit_nofruit='Faces';
                case 2, respstructsingle.prefcat_inhibit_nofruit='Places';
                case 3, respstructsingle.prefcat_inhibit_nofruit='BodyParts';
                case 4, respstructsingle.prefcat_inhibit_nofruit='Objects';
            end
            %%% Raw Selectivity Indices
            % RawSIs
            rsp=respstructsingle.cat_avg([1 3 4 5],2)'; cols=1:4;
            % excite_rawsi
            [val ind]=max(rsp);
            othercols=find(cols~=ind);
            non_ind=mean(rsp(othercols));
            maincol=rsp(ind);
            respstructsingle.excite_rawsi_nofruit=(maincol-non_ind)/(maincol+non_ind);
            % inhibit_rawsi
            [val ind]=min(rsp);
            othercols=find(cols~=ind);
            non_ind=mean(rsp(othercols));
            maincol=rsp(ind);
            respstructsingle.inhibit_rawsi_nofruit=(maincol-non_ind)/(maincol+non_ind);
            % CatSpecificSI
            for c=1:4,
                othercols=find(cols~=c);
                non_ind=mean(rsp(othercols));
                maincol=rsp(c);
                respstructsingle.cat_si_nofruit(c)=(maincol-non_ind)/(maincol+non_ind);
            end
            %%% Solve for traditional face-selectivity
            nonface=mean(respstructsingle.cat_avg([3 4 5],2));
            if respstructsingle.cat_avg(1,2)>(2*nonface), respstructsingle.face_trad_nofruit=1;
            else respstructsingle.face_trad_nofruit=0; end

            %%% Statistical comparison of preferred categories (excited/suppressed)
            %%% to average response to remaining categories - excluding FRUIT
            % Added March 19, 2009
            respstructsingle.cat_id=[ones(20,1);ones(20,1)*2;ones(20,1)*3;ones(20,1)*4;ones(20,1)*5];
            % comparison matrix
            for c1=1:5,
                for c2=1:5,
                    % use mean cat responses
                    pointer1=find(respstructsingle.cat_id==c1);
                    pointer2=find(respstructsingle.cat_id==c2);
                    respstructsingle.stats_rsp_matrix_avg(c1,c2)=ranksum(respstructsingle.m_epoch1(pointer1),respstructsingle.m_epoch1(pointer2));
                    % use trial responses
                    pointer1=find(respstructsingle.trial_id(:,2)==c1);
                    pointer2=find(respstructsingle.trial_id(:,2)==c2);
                    respstructsingle.stats_rsp_matrix_trial(c1,c2)=ranksum(respstructsingle.trial_m_epoch1(pointer1),respstructsingle.trial_m_epoch1(pointer2));
                end
            end
            cats={'Faces','Fruit','Places','BodyParts','Objects'}; catcols=[1 3 4 5]; % eliminate fruits
            % compare preferred excitatory category to remaining categories
            catid=find(strcmp(cats,respstructsingle.prefcat_excite_nofruit)==1);
            otherid=catcols(find(catcols~=catid));
            pointer1=find(respstructsingle.cat_id==catid);
            pointer2=find(ismember(respstructsingle.cat_id,otherid)==1);
            respstructsingle.stats_prefexcite_v_others_nofruit=...
                ranksum(respstructsingle.m_epoch1(pointer1),respstructsingle.m_epoch1(pointer2));

            % compare preferred suppressed category to remaining categories
            catid=find(strcmp(cats,respstructsingle.prefcat_inhibit_nofruit)==1);
            otherid=catcols(find(catcols~=catid));
            pointer1=find(respstructsingle.cat_id==catid);
            pointer2=find(ismember(respstructsingle.cat_id,otherid)==1);
            respstructsingle.stats_prefinhibit_v_others_nofruit=...
                ranksum(respstructsingle.m_epoch1(pointer1),respstructsingle.m_epoch1(pointer2));


            %%%%% SAVE DATA
            respstructsingle.datemodified=date;
            outputfname = [hmiconfig.rsvp500spks,unitname(1:end-4),'-500responsedata.mat'];
            save(outputfname,'respstructsingle')
            outputfname = [hmiconfig.rsvp500spks,unitname(1:end-4),'-500graphdata.mat'];
            disp('Saving average spike density functions...')
            save(outputfname,'graphstructsingle');
            
            %%% graph the file
            disp('Graphing data...')
            plotneuron(hmiconfig,xscale,graphstructsingle,respstructsingle,char(files(f)),metric_col)
            plotneuron_figure(hmiconfig,xscale,graphstructsingle,respstructsingle,char(files(f)),metric_col)
            %%% sync with excel
            if syncexcel==1,
                try
                    respstruct=syncexcel(respstructsingle,hmiconfig,sheetname);
                catch
                    disp('**** Unable to sync with Excel ****');
                end
            end
        end % end loop for each unit

        %         %%% Cross-Correlation Analysis
        %         disp('..conducting cross-correlation analysis...')
        %         range=[-100 100]; % 300 ms on either side
        %         binsize=5; % 5 ms binsize
        %         winsize=length(range(1):binsize:range(2));
        %         numunits=size(spikestruct,2);
        %         for d1=1:numunits, % unit#1
        %             disp(['....correlating ',char(spikestruct(d1).label),' with']);
        %             for d2=1:numunits, % unit#2
        %                 disp(['......',char(spikestruct(d2).label)])
        %                 [xcorr.histograms(d1,d2,1,1:winsize),xcorr.ci(d1,d2,1,:)]=plx_xcorr(spikestruct(d1).faces_spk,spikestruct(d2).faces_spk,range,binsize);
        %                 [xcorr.histograms(d1,d2,2,1:winsize),xcorr.ci(d1,d2,2,:)]=plx_xcorr(spikestruct(d1).fruit_spk,spikestruct(d2).fruit_spk,range,binsize);
        %                 [xcorr.histograms(d1,d2,3,1:winsize),xcorr.ci(d1,d2,3,:)]=plx_xcorr(spikestruct(d1).places_spk,spikestruct(d2).places_spk,range,binsize);
        %                 [xcorr.histograms(d1,d2,4,1:winsize),xcorr.ci(d1,d2,4,:)]=plx_xcorr(spikestruct(d1).objct_spk,spikestruct(d2).objct_spk,range,binsize);
        %                 [xcorr.histograms(d1,d2,5,1:winsize),xcorr.ci(d1,d2,5,:)]=plx_xcorr(spikestruct(d1).bodyp_spk,spikestruct(d2).bodyp_spk,range,binsize);
        %                 alldata1=[spikestruct(d1).faces_spk;spikestruct(d1).fruit_spk;spikestruct(d1).places_spk;...
        %                     spikestruct(d1).objct_spk;spikestruct(d1).bodyp_spk];
        %                 alldata2=[spikestruct(d2).faces_spk;spikestruct(d2).fruit_spk;spikestruct(d2).places_spk;...
        %                     spikestruct(d2).objct_spk;spikestruct(d2).bodyp_spk];
        %                 [xcorr.histograms(d1,d2,6,1:winsize),xcorr.ci(d1,d2,6,:)]=plx_xcorr(alldata1,alldata2,range,binsize);
        %             end
        %         end

        clear tempbehav tempspike spk_avgbaseline m_avgbaseline p_avgbaseline
        %%% Sync Info from Excel %%%


        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

        %plot_xcorr(hmiconfig,xcorr,numunits,range(1):binsize:range(2),spikestruct)
    end % end contains RSVP trials
end % ends file loop
return

function plotneuron(hmiconfig,xscale,graphstruct,respstruct,fname,metric_col)
fontsize_sml=7; fontsize_med=8; fontsize_lrg=9;
%%% determining baseline %%%
switch metric_col
    case 1, avg_baseline=mean(respstruct.spk_baseline); avg_baseline1=mean(respstruct.spk_baseline);
    case 2, avg_baseline=mean(respstruct.m_baseline); avg_baseline1=mean(respstruct.m_baseline);
    case 3, avg_baseline=mean(respstruct.p_baseline); avg_baseline1=mean(respstruct.p_baseline);
    case 4, avg_baseline=mean(respstruct.mp_baseline); avg_baseline1=mean(respstruct.mp_baseline);
    case 5, avg_baseline=mean(respstruct.mp_baseline); avg_baseline1=mean(respstruct.area_baseline);
end
metric_col_list={'SpikeCounts','Mean Spden','Peak Spden','MeanPeak','Area'};
figure
clf; cla; set(gcf,'Units','Normalized'); set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
xrange=(1000+xscale(1)):(1000+xscale(end));
subplot(4,4,[1 5 9 13]) % colour plot
pcolor(xscale,1:100,graphstruct.allconds_avg(:,xrange))
shading flat
%caxis([10 90])
hold on
plot([xscale(1) xscale(end)],[21 21],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[41 41],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[61 61],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[81 81],'w-','LineWidth',1)
plot([0 0],[0 100],'w:','LineWidth',1)
colorbar('SouthOutside')
text(401,graphstruct.bestconds(1)+0.5,['\leftarrow',num2str(graphstruct.bestconds(1))],'FontSize',fontsize_lrg,'FontWeight','Bold')
text(401,graphstruct.bestconds(2)+20.5,['\leftarrow',num2str(graphstruct.bestconds(2))],'FontSize',fontsize_lrg,'FontWeight','Bold')
text(401,graphstruct.bestconds(3)+40.5,['\leftarrow',num2str(graphstruct.bestconds(3))],'FontSize',fontsize_lrg,'FontWeight','Bold')
text(401,graphstruct.bestconds(4)+60.5,['\leftarrow',num2str(graphstruct.bestconds(4))],'FontSize',fontsize_lrg,'FontWeight','Bold')
text(401,graphstruct.bestconds(5)+80.5,['\leftarrow',num2str(graphstruct.bestconds(5))],'FontSize',fontsize_lrg,'FontWeight','Bold')
text(0,101.5,'0','FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(1),101.5,num2str(xscale(1)),'FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(end),101.5,num2str(xscale(end)),'FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(1)-(abs(xscale(1))*.5),10,['Faces (n=',num2str(length(find(respstruct.trial_id(:,2)==1))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),30,['Fruit (n=',num2str(length(find(respstruct.trial_id(:,2)==2))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),50,['Places (n=',num2str(length(find(respstruct.trial_id(:,2)==3))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),70,['Bodyparts (n=',num2str(length(find(respstruct.trial_id(:,2)==4))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),90,['Objects (n=',num2str(length(find(respstruct.trial_id(:,2)==5))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
set(gca,'FontSize',7); xlim([xscale(1) xscale(end)]); box off; axis off; axis ij; ylim([0 100]);
signame=char(graphstruct.label);

subplot(4,4,[2]) % average spike density functions
hold on
plot(xscale,graphstruct.faces_avg(xrange),'r-','LineWidth',2)
plot(xscale,graphstruct.fruit_avg(xrange),'m-','LineWidth',2)
plot(xscale,graphstruct.places_avg(xrange),'b-','LineWidth',2)
plot(xscale,graphstruct.bodyp_avg(xrange),'y-','LineWidth',2)
plot(xscale,graphstruct.objct_avg(xrange),'g-','LineWidth',2)
plot([xscale(1) xscale(end)],[avg_baseline avg_baseline],'k--','LineWidth',0.25)
h=axis;
plot([respstruct.cat_latency(1) respstruct.cat_latency(1)],[0 h(4)],'r-','LineWidth',0.25)
plot([respstruct.cat_latency(2) respstruct.cat_latency(2)],[0 h(4)],'m-','LineWidth',0.25)
plot([respstruct.cat_latency(3) respstruct.cat_latency(3)],[0 h(4)],'b-','LineWidth',0.25)
plot([respstruct.cat_latency(4) respstruct.cat_latency(4)],[0 h(4)],'y-','LineWidth',0.25)
plot([respstruct.cat_latency(5) respstruct.cat_latency(5)],[0 h(4)],'g-','LineWidth',0.25)
plot([0 0],[0 h(4)],'k:','LineWidth',0.5)
plot([xscale(1) xscale(end)],[0 0],'k-')
ylabel('sp/s','FontSize',fontsize_med); set(gca,'FontSize',fontsize_med); xlim([xscale(1) xscale(end)]); box off;
try
    title([signame(1:end-4),': RSVP500 Task (',char(metric_col_list(metric_col)),') ',...
        char(respstruct.gridlocation),' (',char(respstruct.APIndex),') - ',num2str(respstruct.depth),'um (',date,')'],'FontSize',10,'FontWeight','Bold');
catch
    title([signame(1:end-4),': RSVP500 Task (',char(metric_col_list(metric_col)),') (',date,')'],'FontSize',10,'FontWeight','Bold'); % if unable to sync
end

subplot(4,4,3) % raster plots
hold on
rast=plx_makerasters(graphstruct.faces_rast,[0 20],0);
plot(rast(:,1),rast(:,2),'ro','MarkerSize',2,'MarkerFaceColor',[1 0 0])
rast=plx_makerasters(graphstruct.fruit_rast,[20 40],0);
plot(rast(:,1),rast(:,2),'mo','MarkerSize',2,'MarkerFaceColor',[1 0 1])
rast=plx_makerasters(graphstruct.places_rast,[40 60],0);
plot(rast(:,1),rast(:,2),'bo','MarkerSize',2,'MarkerFaceColor',[0 0 1])
rast=plx_makerasters(graphstruct.bodyp_rast,[60 80],0);
plot(rast(:,1),rast(:,2),'yo','MarkerSize',2,'MarkerFaceColor',[1 1 0])
rast=plx_makerasters(graphstruct.objct_rast,[80 100],0);
plot(rast(:,1),rast(:,2),'go','MarkerSize',2,'MarkerFaceColor',[0 1 0])
plot([0 0],[0 100],'k-','LineWidth',0.5)
plot([-200 400],[20 20],'k-','LineWidth',0.5)
plot([-200 400],[40 40],'k-','LineWidth',0.5)
plot([-200 400],[60 60],'k-','LineWidth',0.5)
plot([-200 400],[80 80],'k-','LineWidth',0.5)
set(gca,'YDir','reverse','FontSize',7); xlim([xscale(1) xscale(end)]); ylim([0 100]);

subplot(4,4,6) % average responses
hold on
errorbar(1:5,respstruct.cat_avg(:,metric_col),respstruct.cat_sem(:,metric_col));
bar(1:5,respstruct.cat_avg(:,metric_col))
plot([0.25 5.75],[avg_baseline1 avg_baseline1],'k--','LineWidth',0.25)
ylabel('sp/s','FontSize',fontsize_med); set(gca,'FontSize',fontsize_med);
set(gca,'XTick',1:5); set(gca,'XTickLabels',{'F','Ft','Pl','Bp','Ob'})
title('Average Category Response','FontSize',fontsize_lrg)
h=axis; ylim([0 h(4)]); xlim([0.5 5.5]);
text(0.8,h(4)*0.9,['anova: p=',num2str(respstruct.anova_epoch(1,metric_col),'%1.2g')],'FontSize',fontsize_med)

subplot(4,4,7) % mean latencies
hold on
errorbar(1:5,respstruct.cat_latency(:,1),respstruct.cat_latency(:,2));
bar(1:5,respstruct.cat_latency(:,1))
ylabel('ms','FontSize',fontsize_med); set(gca,'FontSize',fontsize_med);
set(gca,'XTick',1:5); set(gca,'XTickLabels',{'F','Ft','Pl','Bp','Ob'}); ylim([0 200]);
title('Average Response Latency','FontSize',fontsize_lrg); xlim([0.5 5.5]);
h=axis;
text(0.8,h(4)*0.9,['anova: p=',num2str(respstruct.anova_latency,'%1.2g')],'FontSize',fontsize_med)

%%% text details
%%% Within Category Anovas(determines within group selectivity)
text(6.5,190,'Within Category ANOVAs','FontSize',fontsize_lrg,'FontWeight','Bold')
if respstruct.anova_within_group(1,1,2)<0.051,
    text(6.5,150,['Faces:   ',num2str(respstruct.anova_within_group(1,1,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,150,['Faces:   ',num2str(respstruct.anova_within_group(1,1,metric_col),'%1.2g')],'FontSize',fontsize_med); end
if respstruct.anova_within_group(2,1,2)<0.051,
    text(6.5,120,['Fruit:   ',num2str(respstruct.anova_within_group(2,1,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,120,['Fruit:   ',num2str(respstruct.anova_within_group(2,1,metric_col),'%1.2g')],'FontSize',fontsize_med); end
if respstruct.anova_within_group(3,1,2)<0.051,
    text(6.5,90,['Places:  ',num2str(respstruct.anova_within_group(3,1,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,90,['Places:  ',num2str(respstruct.anova_within_group(3,1,metric_col),'%1.2g')],'FontSize',fontsize_med); end
if respstruct.anova_within_group(4,1,2)<0.051,
    text(6.5,60,['BodyP:   ',num2str(respstruct.anova_within_group(4,1,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,60,['BodyP:   ',num2str(respstruct.anova_within_group(4,1,metric_col),'%1.2g')],'FontSize',fontsize_med); end
if respstruct.anova_within_group(5,1,2)<0.051,
    text(6.5,30,['Objects: ',num2str(respstruct.anova_within_group(5,1,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,30,['Objects: ',num2str(respstruct.anova_within_group(5,1,metric_col),'%1.2g')],'FontSize',fontsize_med); end
%%% Valid Responses
text(6.5,-10,'Valid responses according to category?','FontSize',fontsize_lrg,'FontWeight','Bold')
if respstruct.cat_sensory(1,metric_col)<=0.05,
    text(6.5,-50,['Faces:   ',num2str(respstruct.cat_sensory(1,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,-50,['Faces:   ',num2str(respstruct.cat_sensory(1,metric_col),'%1.2g')],'FontSize',fontsize_med); end
if respstruct.cat_sensory(2,metric_col)<=0.05,
    text(6.5,-80,['Fruit:   ',num2str(respstruct.cat_sensory(2,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,-80,['Fruit:   ',num2str(respstruct.cat_sensory(2,metric_col),'%1.2g')],'FontSize',fontsize_med); end
if respstruct.cat_sensory(3,metric_col)<=0.05,
    text(6.5,-110,['Places:  ',num2str(respstruct.cat_sensory(3,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,-110,['Places:  ',num2str(respstruct.cat_sensory(3,metric_col),'%1.2g')],'FontSize',fontsize_med); end
if respstruct.cat_sensory(4,metric_col)<=0.05,
    text(6.5,-140,['BodyP:   ',num2str(respstruct.cat_sensory(4,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,-140,['BodyP:   ',num2str(respstruct.cat_sensory(4,metric_col),'%1.2g')],'FontSize',fontsize_med); end
if respstruct.cat_sensory(5,metric_col)<=0.05,
    text(6.5,-170,['Objects: ',num2str(respstruct.cat_sensory(5,metric_col),'%1.2g')],'FontSize',fontsize_med,'Color','r')
else text(6.5,-170,['Objects: ',num2str(respstruct.cat_sensory(5,metric_col),'%1.2g')],'FontSize',fontsize_med); end

%%% Summary
text(6.5,-210,'Summary (Automated)','FontSize',10,'FontWeight','Bold')
if isempty(find(respstruct.cat_sensory(:,metric_col)<=0.05))~=1, % sensory or not
    if max(respstruct.cat_avg_nobase(:,metric_col))<=0 & respstruct.anova_epoch(1,metric_col)<=0.05,
        text(6.5,-250,'Inhibited (category-selective)','FontSize',fontsize_lrg,'Color','r','FontWeight','Bold');
    elseif max(respstruct.cat_avg_nobase(:,metric_col))<=0 & respstruct.anova_epoch(1,metric_col)>0.05,
        text(6.5,-250,'Inhibited (not selective)','FontSize',fontsize_lrg,'Color','r','FontWeight','Bold');
    elseif max(respstruct.cat_avg_nobase(:,metric_col))>0 & respstruct.anova_epoch(1,metric_col)<=0.05, % selective or not
        text(6.5,-250,'Sensory (category-selective)','FontSize',fontsize_lrg,'Color','g','FontWeight','Bold');
    elseif max(respstruct.cat_avg_nobase(:,metric_col))>0 & respstruct.anova_epoch(1,metric_col)>0.05,
        text(6.5,-250,'Sensory (not selective)','FontSize',fontsize_lrg,'Color','g','FontWeight','Bold');
    end
else
    text(6.5,-250,'Non-responsive','FontSize',9,'Color','r');
end
text(6.5,-290,['Preferred category:     ',respstruct.preferred_category],'FontSize',fontsize_med,'FontWeight','Bold')
text(6.5,-320,['Sel pref: ',num2str(respstruct.raw_si(metric_col),'%1.2g'),' /(no baseline): ',num2str(respstruct.raw_si_nobase(metric_col),'%1.2g')],'FontSize',fontsize_med)
text(6.5,-350,['Type: ',char(respstruct.excitetype),],'FontSize',fontsize_med,'FontWeight','Bold')
text(6.5,-380,['Faces: ',num2str(respstruct.excite_inhibit(1)),],'FontSize',fontsize_med)
text(8.5,-380,['Fruit: ',num2str(respstruct.excite_inhibit(2)),],'FontSize',fontsize_med)
text(10.5,-380,['Places: ',num2str(respstruct.excite_inhibit(3)),],'FontSize',fontsize_med)
text(6.5,-410,['Bodyparts: ',num2str(respstruct.excite_inhibit(4)),],'FontSize',fontsize_med)
text(10.5,-410,['Objects: ',num2str(respstruct.excite_inhibit(5)),],'FontSize',fontsize_med)
try
    text(6.5,-450,'Summary (Confirmed)','FontSize',10,'FontWeight','Bold')
    if strcmp(respstruct.conf_neurtype,'Sensory')==1, % sensory or not
        if strcmp(respstruct.conf_excite,'Inhibit')~=1 & strcmp(respstruct.conf_selective,'Selective')==1,
            text(6.5,-490,'Sensory (category-selective)','FontSize',fontsize_lrg,'Color','g','FontWeight','Bold');
        elseif strcmp(respstruct.conf_excite,'Inhibit')~=1 & strcmp(respstruct.conf_selective,'Not Selective')==1,
            text(6.5,-490,'Sensory (not selective)','FontSize',fontsize_lrg,'Color','g','FontWeight','Bold');
        elseif strcmp(respstruct.conf_excite,'Inhibit')==1 & strcmp(respstruct.conf_selective,'Selective')==1,
            text(6.5,-490,'Inhibited (category-selective)','FontSize',fontsize_lrg,'Color','r','FontWeight','Bold');
        elseif strcmp(respstruct.conf_excite,'Inhibit')==1 & strcmp(respstruct.conf_selective,'Not Selective')==1,
            text(6.5,-490,'Inhibited (not selective)','FontSize',fontsize_lrg,'Color','r','FontWeight','Bold');
        end
    else
        text(6.5,-490,'Non-responsive','FontSize',9,'Color','r');
    end
    text(6.5,-520,['Preferred category:     ',char(respstruct.conf_preferred_cat)],'FontSize',fontsize_med,'FontWeight','Bold')
    text(6.5,-550,['Type:     ',char(respstruct.conf_excite)],'FontSize',fontsize_med,'FontWeight','Bold')
    text(6.5,-580,['Quality:     ',num2str(respstruct.Quality)],'FontSize',fontsize_med,'FontWeight','Bold')
end

subplot(4,4,[10 11]) % selectivity
hold on
bar(1:7,[[respstruct.cat_si(:,metric_col);respstruct.raw_si(metric_col)],[respstruct.cat_si_nobase(:,metric_col);respstruct.raw_si_nobase(metric_col)]],'group')
plot([0.5 7.5],[0.1 0.1],'b:','LineWidth',0.25)
plot([0.5 7.5],[-0.1 -0.1],'b:','LineWidth',0.25)
plot([0.5 7.5],[0.33 0.33],'r:','LineWidth',0.25)
plot([0.5 7.5],[-0.33 -0.33],'r:','LineWidth',0.25)
ylabel('SI','FontSize',fontsize_med); set(gca,'FontSize',fontsize_med);
set(gca,'XTick',1:7); set(gca,'XTickLabels',{'Face','Fruit','Place','Bodyp','Object','FnoFt','RAW'}); ylim([-0.7 0.7]);
if respstruct.face_trad==1, text(0.4,0.4,'Face-Selective (2x)','FontSize',fontsize_lrg,'FontWeight','Bold','Color','r'); end
title('Selectivity Indices (vs. average of all other categories)','FontSize',fontsize_lrg)

subplot(4,4,[14 15]) % pure selectivity
hold on
bar(1:5,respstruct.pure_si,'group')
ylabel('SI','FontSize',fontsize_med); set(gca,'FontSize',fontsize_med);
set(gca,'XTick',1:5); set(gca,'XTickLabels',{'Face','Fruit','Place','Bodyp','Object'}); ylim([-0.7 0.7]);
legend('vs.F','vs.Ft','vs.Pl','vs.Bp','vs.Ob','Orientation','Horizontal','Location','SouthOutside')
title('Pure Selectivity Indices (vs. Single Categories)','FontSize',fontsize_lrg)

%subplot(3,3,3) % pairwise matrix
%matrix=respstruct.pairwise;
%matrix(:,6)=0; matrix(6,:)=0;
%pcolor(0.5:1:5.5,0.5:1:5.5,matrix)
%caxis([0 0.11]); colorbar; axis square; set(gca,'Ydir','reverse')
%set(gca,'YTick',1:5); set(gca,'XTick',1:5)
%set(gca,'YTickLabel',{'F','Ft','Pl','Bp','Ob'});
%set(gca,'XTickLabel',{'F','Ft','Pl','Bp','Ob'});
%set(gca,'FontSize',fontsize_med)
%title('Pairwise Statistical Comparisons','FontSize',fontsize_lrg)

%New panel - shows Spikes
subplot(4,4,4)
if signame(14)=='A',
    wavedata=load([hmiconfig.wave_raw,signame(1:19),'_raw.mat']);
else
    wavedata=load([hmiconfig.wave_raw,signame(1:20),'_raw.mat']);
end
hold on
try plot(-200:25:575,wavedata.waverawdata(:,1:end)','-','Color',[0.5 0.5 0.5],'LineWidth',0.01); end
minval=min(min(wavedata.waverawdata));
maxval=max(max(wavedata.waverawdata));
rangeadj=abs(minval-maxval)*0.1;
plot([(respstruct.wf_params(2)*25)-200 (respstruct.wf_params(2)*25)-200],[minval maxval],'g-')
plot([(respstruct.wf_params(4)*25)-200 (respstruct.wf_params(4)*25)-200],[minval maxval],'g-')
text(200,-1.6,['Duration: ',num2str(respstruct.wf_params(5)),' us'],'FontSize',7)
plot([-200 600],[0 0],'k:'); xlim([-200 600]);
plot(-200:25:575,mean(wavedata.waverawdata'),'r-','LineWidth',2); ylim([minval-rangeadj maxval+rangeadj]);
xlabel('Time (us)'); ylabel('Amplitude (mV)'); set(gca,'FontSize',fontsize_med)
title('Unit Waveforms','FontSize',fontsize_lrg)

%matfigname=[hmiconfig.figure_dir,'rsvp500',filesep,signame(1:end-4),'_rsvp500_',char(metric_col_list(metric_col)),'.fig'];
jpgfigname=[hmiconfig.figure_dir,'rsvp500',filesep,signame(1:end-4),'_rsvp500_',char(metric_col_list(metric_col)),'.jpg'];
%illfigname=[hmiconfig.rootdir,filesep,signame(1:end-4),'_rsvp500_',char(metric_col_list(metric_col)),'.ai'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
%print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
%hgsave(matfigname);
if hmiconfig.printer==1, % prints the figure to the default printer (if printer==1)
    print
end

return




function plotneuron_figure(hmiconfig,xscale,graphstruct,respstruct,fname,metric_col)
% Generate figure suitable for paper/presentations
fontsize_sml=7; fontsize_med=8; fontsize_lrg=9;
figure
clf; cla; set(gcf,'Units','Normalized'); set(gcf,'Position',[0.1 0.1 0.8 0.8]); set(gca,'FontName','Arial')
xrange=(1000+xscale(1)):(1000+xscale(end));
subplot(3,2,1) % average spike density functions
hold on
plot(xscale,graphstruct.faces_avg(xrange),'r-','LineWidth',2)
plot(xscale,graphstruct.places_avg(xrange),'b-','LineWidth',2)
plot(xscale,graphstruct.bodyp_avg(xrange),'y-','LineWidth',2)
plot(xscale,graphstruct.objct_avg(xrange),'g-','LineWidth',2)
h=axis;
plot([0 0],[0 h(4)],'k:','LineWidth',0.5)
plot([xscale(1) xscale(end)],[0 0],'k-')
ylabel('sp/s','FontSize',fontsize_med); set(gca,'FontSize',fontsize_med); xlim([xscale(1) xscale(end)]); box off;
subplot(3,2,[3 5]) % colour plot
% reorganize graphstruct.allconds_avg
fc=graphstruct.allconds_avg(1:20,xrange);
bp=graphstruct.allconds_avg(61:80,xrange);
ob=graphstruct.allconds_avg(81:100,xrange);
pl=graphstruct.allconds_avg(41:60,xrange);
grph=[fc;bp;ob;pl];
pcolor(xscale,1:80,grph)
shading flat
%caxis([10 90])
hold on; %colormap(hot)
plot([xscale(1) xscale(end)],[21 21],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[41 41],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[61 61],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[81 81],'w-','LineWidth',1)
plot([0 0],[0 80],'w:','LineWidth',1)
colorbar('SouthOutside')
text(0,81.5,'0','FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(1),81.5,num2str(xscale(1)),'FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(end),81.5,num2str(xscale(end)),'FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(1)-(abs(xscale(1))*.5),10,['Faces (n=',num2str(length(find(respstruct.trial_id(:,2)==1))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),70,['Places (n=',num2str(length(find(respstruct.trial_id(:,2)==3))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),30,['Bodyparts (n=',num2str(length(find(respstruct.trial_id(:,2)==4))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),50,['Objects (n=',num2str(length(find(respstruct.trial_id(:,2)==5))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
set(gca,'FontSize',7); xlim([xscale(1) xscale(end)]); box off; axis off; axis ij; ylim([0 80]);
signame=char(graphstruct.label);

subplot(3,2,[4 6]) % colour plot
% reorganize graphstruct.allconds_avg
normalizer=max([max(max(fc)) max(max(bp)) max(max(ob)) max(max(pl))]);
normalizer=mean(respstruct.trial_m_baseline);
baselinenorm=mean(respstruct.trial_m_baseline)/normalizer;
fc_norm=(fc/normalizer);
bp_norm=(bp/normalizer);
ob_norm=(ob/normalizer);
pl_norm=(pl/normalizer);
grph=[fc_norm;bp_norm;ob_norm;pl_norm];
pcolor(xscale,1:80,grph)
shading flat
caxis([0 10]);
hold on
plot([xscale(1) xscale(end)],[21 21],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[41 41],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[61 61],'w-','LineWidth',1)
plot([xscale(1) xscale(end)],[81 81],'w-','LineWidth',1)
plot([0 0],[0 80],'w:','LineWidth',1)
colorbar('SouthOutside')
text(0,81.5,'0','FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(1),81.5,num2str(xscale(1)),'FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(end),81.5,num2str(xscale(end)),'FontSize',fontsize_sml,'HorizontalAlignment','Center')
text(xscale(1)-(abs(xscale(1))*.5),10,['Faces (n=',num2str(length(find(respstruct.trial_id(:,2)==1))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),70,['Places (n=',num2str(length(find(respstruct.trial_id(:,2)==3))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),30,['Bodyparts (n=',num2str(length(find(respstruct.trial_id(:,2)==4))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
text(xscale(1)-(abs(xscale(1))*.5),50,['Objects (n=',num2str(length(find(respstruct.trial_id(:,2)==5))),')'],'FontSize',fontsize_med,'HorizontalAlignment','Center','Rotation',90)
set(gca,'FontSize',7); xlim([xscale(1) xscale(end)]); box off; axis off; axis ij; ylim([0 80]);
signame=char(graphstruct.label);
jpgfigname=[hmiconfig.rootdir,filesep,signame(1:end-4),'.jpg'];
illfigname=[hmiconfig.rootdir,filesep,signame(1:end-4),'.ai'];
print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
return






function [av,sm,bst]=rsp_avgbest(data,stim_no);
%%% one column for each metric (currently 3)
%[av(1),sm(1)]=mean_sem(data.spk_epoch1(stim_no));
[av(2),sm(2)]=mean_sem(data.m_epoch1(stim_no));
%[av(3),sm(3)]=mean_sem(data.p_epoch1(stim_no));
%[av(4),sm(4)]=mean_sem(data.mp_epoch1(stim_no));
%[av(5),sm(5)]=mean_sem(data.area_epoch1(stim_no));
%bst(1)=max(data.spk_epoch1(stim_no));
bst(2)=max(data.m_epoch1(stim_no));
%bst(3)=max(data.p_epoch1(stim_no));
%bst(4)=max(data.mp_epoch1(stim_no));
%bst(5)=max(data.area_epoch1(stim_no));
return

function [av,sm,bst]=rsp_avgbest_nobase(data,stim_no);
%%% one column for each metric (currently 3)
%[av(1),sm(1)]=mean_sem(data.spk_epoch1_nobase(stim_no));
[av(2),sm(2)]=mean_sem(data.m_epoch1_nobase(stim_no));
%[av(3),sm(3)]=mean_sem(data.p_epoch1_nobase(stim_no));
%[av(4),sm(4)]=mean_sem(data.mp_epoch1_nobase(stim_no));
%[av(5),sm(5)]=mean_sem(data.area_epoch1_nobase(stim_no));
%bst(1)=max(data.spk_epoch1_nobase(stim_no));
bst(2)=max(data.m_epoch1_nobase(stim_no));
%bst(3)=max(data.p_epoch1_nobase(stim_no));
%bst(4)=max(data.mp_epoch1_nobase(stim_no));
%bst(5)=max(data.area_epoch1_nobase(stim_no));
return

function si=calc_si(data,ind_col,metric_col)
[numcol,junk]=size(data);
cols=1:numcol;
othercols=find(cols~=ind_col);
non_ind=mean(data(othercols,metric_col));
maincol=data(ind_col,metric_col);
si=(maincol-non_ind)/(maincol+non_ind);
return

function si=calc_si_nofruit(data,ind_col,metric_col)
[numcol,junk]=size(data);
cols=1:numcol;
non_ind=mean(data([2 3 4],metric_col));
maincol=data(ind_col,metric_col);
si=(maincol-non_ind)/(maincol+non_ind);
return

function si=calc_rawsi(data,metric_col)
[numcol,junk]=size(data);
cols=1:numcol;
[val,ind]=max(data(:,metric_col));
othercols=find(cols~=ind);
non_ind=mean(data(othercols,metric_col));
maincol=data(ind,metric_col);
si=(maincol-non_ind)/(maincol+non_ind);
return

function output=calc_puresi(data,metric_col,refcol);
maincol=data(refcol,metric_col);
for cc=1:5,
    noncol=data(cc,metric_col);
    output(cc)=(maincol-noncol)/(maincol+noncol);
end
return

function matrix=pairwisematrix(datastruct,epoch,metric_col);
catbin=1:20:101;
matrix=zeros(5,5);
switch epoch
    case 1
        dataspk=datastruct.spk_epoch1;
        datam=datastruct.m_epoch1;
        datap=datastruct.p_epoch1;
        datamp=datastruct.mp_epoch1;
    case 2
        dataspk=datastruct.spk_epoch2;
        datam=datastruct.m_epoch2;
        datap=datastruct.p_epoch2;
        datamp=datastruct.mp_epoch2;
    case 3
        dataspk=datastruct.spk_epoch3;
        datam=datastruct.m_epoch3;
        datap=datastruct.p_epoch3;
        datamp=datastruct.mp_epoch3;
end
switch metric_col
    case 1
        for cat=1:5,
            matrix(cat,1)=ranksum(dataspk(catbin(cat):catbin(cat)+19),dataspk(1:20));
            matrix(cat,2)=ranksum(dataspk(catbin(cat):catbin(cat)+19),dataspk(21:40));
            matrix(cat,3)=ranksum(dataspk(catbin(cat):catbin(cat)+19),dataspk(41:60));
            matrix(cat,4)=ranksum(dataspk(catbin(cat):catbin(cat)+19),dataspk(61:80));
            matrix(cat,5)=ranksum(dataspk(catbin(cat):catbin(cat)+19),dataspk(81:100));
        end
    case 2
        for cat=1:5,
            matrix(cat,1)=ranksum(datam(catbin(cat):catbin(cat)+19),datam(1:20));
            matrix(cat,2)=ranksum(datam(catbin(cat):catbin(cat)+19),datam(21:40));
            matrix(cat,3)=ranksum(datam(catbin(cat):catbin(cat)+19),datam(41:60));
            matrix(cat,4)=ranksum(datam(catbin(cat):catbin(cat)+19),datam(61:80));
            matrix(cat,5)=ranksum(datam(catbin(cat):catbin(cat)+19),datam(81:100));
        end
    case 3
        for cat=1:5,
            matrix(cat,1)=ranksum(datap(catbin(cat):catbin(cat)+19),datap(1:20));
            matrix(cat,2)=ranksum(datap(catbin(cat):catbin(cat)+19),datap(21:40));
            matrix(cat,3)=ranksum(datap(catbin(cat):catbin(cat)+19),datap(41:60));
            matrix(cat,4)=ranksum(datap(catbin(cat):catbin(cat)+19),datap(61:80));
            matrix(cat,5)=ranksum(datap(catbin(cat):catbin(cat)+19),datap(81:100));
        end
    case 4
        for cat=1:5,
            matrix(cat,1)=ranksum(datamp(catbin(cat):catbin(cat)+19),datamp(1:20));
            matrix(cat,2)=ranksum(datamp(catbin(cat):catbin(cat)+19),datamp(21:40));
            matrix(cat,3)=ranksum(datamp(catbin(cat):catbin(cat)+19),datamp(41:60));
            matrix(cat,4)=ranksum(datamp(catbin(cat):catbin(cat)+19),datamp(61:80));
            matrix(cat,5)=ranksum(datamp(catbin(cat):catbin(cat)+19),datamp(81:100));
        end
end
return

function output=syncexcel(respstruct,hmiconfig,sheetname);
output=respstruct;
[crap,PlxFile]=xlsread(hmiconfig.excelfile,sheetname,'B1:B1000'); % alpha, PlexonFilename
[crap,UnitName]=xlsread(hmiconfig.excelfile,sheetname,'C1:C1000'); % alpha, Unitname
disp('..Syncing with Excel...')
%% find matching units
tempname=char(respstruct(1).label);
tempname=[tempname(1:12),'.plx'];
matches=find(strcmp(tempname,PlxFile)==1);
[crap,PlxFile]=xlsread(hmiconfig.excelfile,sheetname,['B',num2str(matches(1)),':B',num2str(matches(end))]); % alpha, PlexonFilename
[crap,UnitName]=xlsread(hmiconfig.excelfile,sheetname,['C',num2str(matches(1)),':C',num2str(matches(end))]); % alpha, Unitname
[crap,GridLoc]=xlsread(hmiconfig.excelfile,sheetname,['E',num2str(matches(1)),':E',num2str(matches(end))]); % alphanumeric, Gridlocation
Depth=xlsread(hmiconfig.excelfile,sheetname,['F',num2str(matches(1)),':F',num2str(matches(end))]); % numeric, Depth
[crap,APIndex]=xlsread(hmiconfig.excelfile,sheetname,['G',num2str(matches(1)),':G',num2str(matches(end))]); % alphanumeric, APindex
[crap,NeurType]=xlsread(hmiconfig.excelfile,sheetname,['J',num2str(matches(1)),':J',num2str(matches(end))]); % alphanumeric, neuron type
[crap,ConfPref]=xlsread(hmiconfig.excelfile,sheetname,['L',num2str(matches(1)),':L',num2str(matches(end))]); % alphanumeric, pref category
[crap,ConfSele]=xlsread(hmiconfig.excelfile,sheetname,['N',num2str(matches(1)),':N',num2str(matches(end))]); % alphanumeric, conf selective
[crap,ConfExci]=xlsread(hmiconfig.excelfile,sheetname,['P',num2str(matches(1)),':P',num2str(matches(end))]); % alphanumeric, conf excite
Quality=xlsread(hmiconfig.excelfile,sheetname,['Q',num2str(matches(1)),':Q',num2str(matches(end))]); % alphanumeric, conf excite

for mt=1:length(matches),
    tempname=char(PlxFile(mt));
    for un=1:size(output,2),
        unitn=char(UnitName(mt));
        if strcmp([tempname(1:12),'-',unitn(1:7),'.mat'],output(un).label)==1
            output(un).gridlocation=GridLoc(mt);
            output(un).depth=Depth(mt);
            output(un).APIndex=APIndex(mt);
            try
                output(un).conf_neurtype=NeurType(mt);
                output(un).conf_preferred_cat=ConfPref(mt);
                output(un).conf_selective=ConfSele(mt);
                output(un).conf_excite=ConfExci(mt);
                output(un).Quality=Quality(mt);
            catch
                disp('*** Unable to sync ''confirmed'' classifications')
            end
            %%% Paste into Excel automated values
            if isempty(find(output(un).cat_sensory(:,2)<=0.05))~=1, % sensory or not
                if max(output(un).cat_avg_nobase(:,2))<=0,
                    neurontype='Sensory'; % neurontype='Inhibited';
                elseif max(output(un).cat_avg_nobase(:,2))>0,
                    neurontype='Sensory';
                else
                    neurontype='Non-Responsive';
                end
            else
                neurontype='Non-Responsive';
            end
            excite=output(un).excitetype;
            prefcat=output(un).preferred_category;
            if output(un).anova_epoch(1,2)<=0.05, % determine if ANOVA for category responses is significant
                selective='Selective';
            else
                selective='Not Selective';
            end
            xlswrite(hmiconfig.excelfile,{neurontype},sheetname,['I',num2str(matches(mt)),':I',num2str(matches(mt))])
            xlswrite(hmiconfig.excelfile,{prefcat},sheetname,['K',num2str(matches(mt)),':K',num2str(matches(mt))])
            xlswrite(hmiconfig.excelfile,{selective},sheetname,['M',num2str(matches(mt)),':M',num2str(matches(mt))])
            xlswrite(hmiconfig.excelfile,{excite},sheetname,['O',num2str(matches(mt)),':O',num2str(matches(mt))])
        end
    end
end
return

function output=prep_raster(spks,ts,xscale);
output=spks*0; bleh=[];
for x=1:length(ts),
    bleh(x,:)=spks(x,:)-ts(x);
    pointer=find(ismember(bleh(x,:),xscale)==1);
    output(x,1:length(pointer))=bleh(x,pointer);
end
output=output(isnan(output(:,1))==0,:);
return
