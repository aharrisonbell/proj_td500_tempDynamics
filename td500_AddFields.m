function td500_AddFields(initial,criterion)
% f570_AddFields(files); %
% written by AHB, May 2014
global lsnconfig 

%%%  LOAD FILE LIST
if nargin==0,
    error('You must specify an individual filename or monkey initial (''S''/''W'').')
elseif strcmp(initial,'S')==1
    disp('Analyzing all RSVP500 files for Stewie...')
    try
        % Pulls files from HMI_PhysiologyNotes
        include=xlsread(lsnconfig.excelfile,'RSVP500','A10:A1000'); % alphanumeric, Gridlocation
        [~,filest]=xlsread(lsnconfig.excelfile,'RSVP500','B10:B1000');
    catch
        disp('*** Unable to load Excel Spreadsheet - Using matlab file instead ***')
        load([lsnconfig.datadir,'Stewie_FileList.mat']);
    end
    filesx=filest(ismember(include,criterion)); %%% CHANGE CRITERION HERE
    if isempty(filesx)==1,
        disp('* No files marked for processing')
        return
    end
    for ff=1:size(filesx,1),
        temp=filesx{ff};
        files(ff)=cellstr(temp(1:12));
    end
    save([lsnconfig.datadir,'Stewie_FileList.mat'],'filest','include');
elseif strcmp(initial,'W')==1
    disp('Analyzing all RSVP500 files for Wiggum...')
    try
        % Pulls files from HMI_PhysiologyNotes
        include=xlsread(lsnconfig.excelfile,'RSVP500','D10:D1000'); % alphanumeric, Gridlocation
        [~,filest]=xlsread(lsnconfig.excelfile,'RSVP500','E10:E1000');
    catch
        disp('*** Unable to load Excel Spreadsheet - Using matlab file instead ***')
        load([config.datadir,'Wiggum_FileList.mat']);
    end
    filesx=filest(ismember(include,criterion)); %%% CHANGE CRITERION HERE
    if isempty(filesx)==1,
        disp('* No files marked for processing')
        return
    end
    for ff=1:size(filesx,1),
        temp=filesx{ff};
        files(ff)=cellstr(temp(1:12));
    end
    save([lsnconfig.datadir,'Wiggum_FileList.mat'],'filest','include');
end
clear include initial

disp('**************************************************************************************')
disp('td500_AddFields.m - Add fields to previously created RSVP500 respstructsingle matrices')
disp('**************************************************************************************')
for f=1:length(files), % perform following operations on each nex file listed
    close all % close all figure windows
    filename=char(files(f));
    disp(['Modifying files from ',filename,' (',num2str(f),' of ',num2str(length(files)),')'])
    
    % Begin analysis
    
    
    try
        tempstruct=load([lsnconfig.spikedir,filename,'_spkmat.mat']);
    catch
        disp('Unable to find SPKMAT.MAT - rerunning PLX_PROCESSNEXFILE')
        plx_processnexfile({filename},0);
        plx500({filename})
        close all;
        tempstruct=load([lsnconfig.spikedir,filename,'_spkmat.mat']);
    end
    
    tempbehav=tempstruct.behav_matrix(:,[1 3 4 30 40 44]); % load behavioural data
    tempbehav(:,7)=tempbehav(:,6)-tempbehav(:,5); % solve for cue onset time (aligned to the beginning of each trial, in ms?)
    tempspike=tempstruct.spikesig;
    clear tempstruct
    foundunits=size(tempspike,2);
    if length(find(tempbehav(:,2)==500))<1,
        disp(['.No RSVP500 trials found!!  Skipping this file.'])
    else
        disp(['.found ',num2str(size(tempbehav(tempbehav(:,2)==500,:),1)),' trials...'])
        disp(['.found ',num2str(foundunits),' unit(s)...'])
 
        % Unit loop
        for un=1:foundunits, % performed for each unit
            disp(['..analyzing ',char(tempspike(un).labels)])
            newname=char(tempspike(un).labels); fullunitname=[filename,newname(13:20)];
            load([lsnconfig.rsvp500spks,fullunitname,'-500graphdata.mat']);
            load([lsnconfig.rsvp500spks,fullunitname,'-500responsedata.mat']);
            
            % Add Fields
            % EPSP and shorter Gaussian spike density functions
            % disp('...Creating new spike density functions...')
            
            %%% Setup Spike Density Structures
            spikestruct=struct('label',[],'faces_spk',[],'fruit_spk',[],'bodyp_spk',[],'places_spk',[],'objct_spk',[],...
                'faces_ts',[],'fruit_ts',[],'bodyp_ts',[],'places_ts',[],'objct_ts',[]);
            
            graphstructsingle.label=fullunitname; % paste label into GRAPH structure
            respstructsingle.label=fullunitname; % paste label into RESPONSE DATA structure
            spikestruct.label=fullunitname; % paste label into SPIKE DENSITY FUNCTION structure
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% GENERATE SPIKE DENSITY FUNCTIONS %%%
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            graphstructsingle.spden_trial=[];
            for tr=1:size(tempbehav,1),
                respstructsingle.trial_id(tr,1)=tempbehav(tr,3); % paste stimulus number
                switch tempbehav(tr,3) %% assign category number
                    case {1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20}, respstructsingle.trial_id(tr,2)=1;
                    case {21,22,23,24,25,26,27,28,29,30,31,32,33,34,35,36,37,38,39,40}, respstructsingle.trial_id(tr,2)=2;
                    case {41,42,43,44,45,46,47,48,49,50,51,52,53,54,55,56,57,58,59,60}, respstructsingle.trial_id(tr,2)=3;
                    case {61,62,63,64,65,66,67,68,69,70,71,72,73,74,75,76,77,78,79,80}, respstructsingle.trial_id(tr,2)=4;
                    case {81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98,99,100}, respstructsingle.trial_id(tr,2)=5;
                end
                %%% Generate spike density function for each TRIAL
                tempspikes=tempspike(un).spikes(tr,:);
                temp_ts=ceil(tempbehav(tr,7)*1000); % round cue onset timestamps to nearest ms
                
                [graphstructsingle.spden_trial(tr,:),~]=plx_avgspden(tempspikes,temp_ts,601,1,lsnconfig.gausskernel,lsnconfig.xscale);
            end
            
            
            pointer=find(ismember(tempbehav(:,2),500)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),lsnconfig.faces500)==1); % select only correct 500 series trials - faces
            spikestruct.faces_spk=tempspike(un).spikes(pointer,:);
            spikestruct.faces_ts=ceil(tempbehav(pointer,7)*1000); % round cue onset timestamps to nearest ms
            spikestruct.faces_rast=plx_prepRaster(spikestruct.faces_spk,spikestruct.faces_ts,lsnconfig.xscale);
            
            pointer=find(ismember(tempbehav(:,2),500)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),lsnconfig.fruit500)==1); % select only correct 500 series trials -
            spikestruct.fruit_spk=tempspike(un).spikes(pointer,:);
            spikestruct.fruit_ts=ceil(tempbehav(pointer,7)*1000);
            spikestruct.fruit_rast=plx_prepRaster(spikestruct.fruit_spk,spikestruct.fruit_ts,lsnconfig.xscale);
            
            pointer=find(ismember(tempbehav(:,2),500)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),lsnconfig.bodyp500)==1); % select only correct 500 series trials
            spikestruct.bodyp_spk=tempspike(un).spikes(pointer,:);
            spikestruct.bodyp_ts=ceil(tempbehav(pointer,7)*1000);
            spikestruct.bodyp_rast=plx_prepRaster(spikestruct.bodyp_spk,spikestruct.bodyp_ts,lsnconfig.xscale);
            
            pointer=find(ismember(tempbehav(:,2),500)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),lsnconfig.places500)==1); % select only correct 500 series trials
            spikestruct.places_spk=tempspike(un).spikes(pointer,:);
            spikestruct.places_ts=ceil(tempbehav(pointer,7)*1000);
            spikestruct.places_rast=plx_prepRaster(spikestruct.places_spk,spikestruct.places_ts,lsnconfig.xscale);
            
            pointer=find(ismember(tempbehav(:,2),500)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),lsnconfig.objct500)==1); % select only correct 500 series trials
            spikestruct.objct_spk=tempspike(un).spikes(pointer,:);
            spikestruct.objct_ts=ceil(tempbehav(pointer,7)*1000);
            spikestruct.objct_rast=plx_prepRaster(spikestruct.objct_spk,spikestruct.objct_ts,lsnconfig.xscale);
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %%% Generate average spike density functions for each CONDITION (listed in CAT_avg/CAT_sem)
            % Gaussian kernel: 5ms
            [spikestruct.faces_avg_smlK,spikestruct.faces_sem_smlK]  =plx_avgspden(spikestruct.faces_spk,spikestruct.faces_ts,601,1,lsnconfig.gausskernelsml,lsnconfig.xscale);
            [spikestruct.fruit_avg_smlK,spikestruct.fruit_sem_smlK]  =plx_avgspden(spikestruct.fruit_spk,spikestruct.fruit_ts,601,1,lsnconfig.gausskernelsml,lsnconfig.xscale);
            [spikestruct.bodyp_avg_smlK,spikestruct.bodyp_sem_smlK]  =plx_avgspden(spikestruct.bodyp_spk,spikestruct.bodyp_ts,601,1,lsnconfig.gausskernelsml,lsnconfig.xscale);
            [spikestruct.places_avg_smlK,spikestruct.places_sem_smlK]=plx_avgspden(spikestruct.places_spk,spikestruct.places_ts,601,1,lsnconfig.gausskernelsml,lsnconfig.xscale);
            [spikestruct.objct_avg_smlK,spikestruct.objct_sem_smlK]  =plx_avgspden(spikestruct.objct_spk,spikestruct.objct_ts,601,1,lsnconfig.gausskernelsml,lsnconfig.xscale);
            
            % Gaussian kernel: 10ms
            [spikestruct.faces_avg,spikestruct.faces_sem]  =plx_avgspden(spikestruct.faces_spk,spikestruct.faces_ts,601,1,lsnconfig.gausskernel,lsnconfig.xscale);
            [spikestruct.fruit_avg,spikestruct.fruit_sem]  =plx_avgspden(spikestruct.fruit_spk,spikestruct.fruit_ts,601,1,lsnconfig.gausskernel,lsnconfig.xscale);
            [spikestruct.bodyp_avg,spikestruct.bodyp_sem]  =plx_avgspden(spikestruct.bodyp_spk,spikestruct.bodyp_ts,601,1,lsnconfig.gausskernel,lsnconfig.xscale);
            [spikestruct.places_avg,spikestruct.places_sem]=plx_avgspden(spikestruct.places_spk,spikestruct.places_ts,601,1,lsnconfig.gausskernel,lsnconfig.xscale);
            [spikestruct.objct_avg,spikestruct.objct_sem]  =plx_avgspden(spikestruct.objct_spk,spikestruct.objct_ts,601,1,lsnconfig.gausskernel,lsnconfig.xscale);
            
            % EPSP kernel: 1/20ms
            [spikestruct.faces_avg_epsp,spikestruct.faces_sem_epsp]  =plx_avgspden(spikestruct.faces_spk,spikestruct.faces_ts,601,2,lsnconfig.gausskernelsml,lsnconfig.xscale);
            [spikestruct.fruit_avg_epsp,spikestruct.fruit_sem_epsp]  =plx_avgspden(spikestruct.fruit_spk,spikestruct.fruit_ts,601,2,lsnconfig.gausskernelsml,lsnconfig.xscale);
            [spikestruct.bodyp_avg_epsp,spikestruct.bodyp_sem_epsp]  =plx_avgspden(spikestruct.bodyp_spk,spikestruct.bodyp_ts,601,2,lsnconfig.gausskernelsml,lsnconfig.xscale);
            [spikestruct.places_avg_epsp,spikestruct.places_sem_epsp]=plx_avgspden(spikestruct.places_spk,spikestruct.places_ts,601,2,lsnconfig.gausskernelsml,lsnconfig.xscale);
            [spikestruct.objct_avg_epsp,spikestruct.objct_sem_epsp]  =plx_avgspden(spikestruct.objct_spk,spikestruct.objct_ts,601,2,lsnconfig.gausskernelsml,lsnconfig.xscale);
            
            % Condition Loop
            % disp('.....Scrolling through each condition...')
            for cnd=1:100, % first loop creates average spike density function for each condition
                pointer=find(ismember(tempbehav(:,2),500)==1 & tempbehav(:,4)==6 & ismember(tempbehav(:,3),cnd)==1);
                if isempty(pointer)==1, % if no trials are found, paste zero values
                    spikestruct.allconds_avg_smlK(cnd,:)=zeros(1,601);
                    spikestruct.allconds_sem_smlK(cnd,:)=zeros(1,601);
                    spikestruct.allconds_avg(cnd,:)=zeros(1,601);
                    spikestruct.allconds_sem(cnd,:)=zeros(1,601);
                    spikestruct.allconds_avg_epsp(cnd,:)=zeros(1,601);
                    spikestruct.allconds_sem_epsp(cnd,:)=zeros(1,601);
                else
                    tempspikes=tempspike(un).spikes(pointer,:);
                    temp_ts=ceil(tempbehav(pointer,7)*1000); % round cue onset timestamps to nearest ms
                    [spikestruct.allconds_avg(cnd,:),spikestruct.allconds_sem(cnd,:)]=...
                        plx_avgspden(tempspikes,temp_ts,601,1,lsnconfig.gausskernel,lsnconfig.xscale);
                    [spikestruct.allconds_avg_smlK(cnd,:),spikestruct.allconds_sem_smlK(cnd,:)]=...
                        plx_avgspden(tempspikes,temp_ts,601,1,lsnconfig.gausskernelsml,lsnconfig.xscale);
                    [spikestruct.allconds_avg_epsp(cnd,:),spikestruct.allconds_sem_epsp(cnd,:)]=...
                        plx_avgspden(tempspikes,temp_ts,601,2,lsnconfig.gausskernel,lsnconfig.xscale);
                end
            end
            
            % Resave output file for population analysis
            respstructsingle.datemodified=date;
            graphstructsingle.datemodified=date;
            spikestruct.datemodified=date;
            
            unitname=char(respstructsingle.label);
            outputfname = [lsnconfig.rsvp500spks,unitname,'-500_NeuronData.mat'];
            save(outputfname,'respstructsingle','graphstructsingle','spikestruct');
            
            clear respstructsingle graphstructsingle spikestruct
        end % end loop for each unit
        clear tempbehav tempspike
    end
end
return
