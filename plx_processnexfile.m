function plx_processnexfile(sheetname,manual)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% plx_processnexfile(filez); %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% based on readstevenexfilegg.m, written by Steve Gotts (NIMH)
% edited by AHB, Sept2007, Modified Nov2007
% Modified February 12, 2016
% Creates matrices for analog, spike, and strobed channels from nex file
% NEX file must first be created from PLX file using OFS
% filez = optional argument, list filez as strings.  Otherwise, program
% will load files listed in default nexlist.txt
%
% This program uses a nested program written by PLEXON (readNexFile.m) and
% then arranges the output into the desired files

%%% SETUP DEFAULTS
hmiconfig = generate_hmi_configplex; % generates and loads config file
% channels = {'LFP1','LFP2','LFP3','LFP4','AD05','AD06','AD07','AD08','AD09','AD10','AD11','AD12','AD13','AD14'};
channels = {'LFP1','LFP2','LFP3','LFP4','AD14'}; % channel defaults
behav_only=0;

%%%  LOAD FILE LIST
if nargin==0,
    include=xlsread(hmiconfig.plxlist,hmiconfig.sheetnex,'A2:A1000'); % alphanumeric, Gridlocation
    manual=xlsread(hmiconfig.plxlist,hmiconfig.sheetnex,'B2:B1000'); % alphanumeric, Gridlocation
    [~,filez]=xlsread(hmiconfig.plxlist,hmiconfig.sheetnex,'C2:C1000'); % numeric, Depth
    lfpopt=xlsread(hmiconfig.plxlist,hmiconfig.sheetnex,'D2:D1000'); % numeric, lfpopt
    filez=filez(include==1); manual=manual(include==1); lfpopt=lfpopt(include==1);
    clear include
    skiplfpeye=0; %#ok<*NASGU>
elseif nargin==1,
    include=xlsread(hmiconfig.plxlist,sheetname,'A2:A800'); % alphanumeric, Gridlocation
    manual=xlsread(hmiconfig.plxlist,sheetname,'B2:B800'); % alphanumeric, Gridlocation
    [~,filez]=xlsread(hmiconfig.plxlist,sheetname,'C2:C800'); % numeric, Depth
    [~,filez]=xlsread(hmiconfig.plxlist,hmiconfig.sheetnex,'D2:D1000'); % numeric, lfpopt
    filez=filez(include==1); manual=manual(include==1); lfpopt=lfpopt(include==1); %#ok<NODEF>
    clear include
    skiplfpeye=0;
elseif nargin==2,
    filez=sheetname;
    lfpopt=2;
end

for f=1:length(filez), % perform following operations on each nex file listed
    %%% LOAD INDIVIDUAL FILE %%%
    close all
    disp(' ')
    cprintf('*green',['Processing ', char(filez(f)),'\n']);
    disp('*****************************************')
    tempname=char(filez(f));
    nexfname=strcat(hmiconfig.indir,filez(f),'.nex');
    outputmatname=[hmiconfig.outdir,tempname,'_NexFile.mat'];
    nexFile = []; % create empty matrix
    disp('Loading nexfile')
    fid = fopen(char(nexfname), 'r');
    if(fid == -1)
        error 'Unable to open file'
        return %#ok<UNRCH>
    end
    magic = fread(fid, 1, 'int32');
    if magic ~= 827868494
        error 'The file is not a valid .nex file'
    end
    nexFile.version = fread(fid, 1, 'int32');
    nexFile.comment = deblank(char(fread(fid, 256, 'char')'));
    nexFile.freq    = fread(fid, 1, 'double');
    nexFile.tbeg    = fread(fid, 1, 'int32')./nexFile.freq;
    nexFile.tend    = fread(fid, 1, 'int32')./nexFile.freq;
    nvar = fread(fid, 1, 'int32');
    % Display comments
    disp(['Experiment length: ',num2str(round(nexFile.tend)/60), ' min '])
    disp(['Comments contained in ',char(filez(f)),':'])
    disp('-----------------------------------------------------')
    disp(nexFile.comment)
    disp('-----------------------------------------------------')
    % skip location of next header and padding
    fseek(fid, 260, 'cof');
    neuronCount = 0;
    eventCount = 0;
    intervalCount = 0;
    waveCount = 0;
    popCount = 0;
    contCount = 0;
    markerCount = 0;

    % read all variables
    for i=1:nvar
        type = fread(fid, 1, 'int32');
        varVersion = fread(fid, 1, 'int32');
        name = deblank(char(fread(fid, 64, 'char')'));
        offset = fread(fid, 1, 'int32');
        n = fread(fid, 1, 'int32');
        WireNumber = fread(fid, 1, 'int32');
        UnitNumber = fread(fid, 1, 'int32');
        Gain = fread(fid, 1, 'int32');
        Filter = fread(fid, 1, 'int32');
        XPos = fread(fid, 1, 'double');
        YPos = fread(fid, 1, 'double');
        WFrequency = fread(fid, 1, 'double');   % wf sampling fr.
        ADtoMV  = fread(fid, 1, 'double');      % coeff to convert from AD values to Millivolts.
        NPointsWave = fread(fid, 1, 'int32');   % number of points in each wave
        NMarkers = fread(fid, 1, 'int32');      % how many values are associated with each marker
        MarkerLength = fread(fid, 1, 'int32');  % how many characters are in each marker value
        MVOfffset = fread(fid, 1, 'double');    % coeff to shift AD values in Millivolts: mv = raw*ADtoMV+MVOfffset
        filePosition = ftell(fid);
        switch type
            case 0 % neuron
                neuronCount = neuronCount+1;
                nexFile.neurons{neuronCount,1}.name = name;
                fseek(fid, offset, 'bof');
                nexFile.neurons{neuronCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
                fseek(fid, filePosition, 'bof');
            case 1 % event
                eventCount = eventCount+1;
                nexFile.events{eventCount,1}.name = name;
                fseek(fid, offset, 'bof');
                nexFile.events{eventCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
                fseek(fid, filePosition, 'bof');
            case 2 % interval
                intervalCount = intervalCount+1;
                nexFile.intervals{intervalCount,1}.name = name;
                fseek(fid, offset, 'bof');
                nexFile.intervals{intervalCount,1}.intStarts = fread(fid, [n 1], 'int32')./nexFile.freq;
                nexFile.intervals{intervalCount,1}.intEnds = fread(fid, [n 1], 'int32')./nexFile.freq;
                fseek(fid, filePosition, 'bof');
            case 3 % waveform
                waveCount = waveCount+1;
                nexFile.waves{waveCount,1}.name = name;
                nexFile.waves{waveCount,1}.NPointsWave = NPointsWave;
                nexFile.waves{waveCount,1}.WFrequency = WFrequency;
                fseek(fid, offset, 'bof');
                nexFile.waves{waveCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
                nexFile.waves{waveCount,1}.waveforms = fread(fid, [NPointsWave n], 'int16').*ADtoMV + MVOfffset;
                fseek(fid, filePosition, 'bof');
            case 4 % population vector
                popCount = popCount+1;
                nexFile.popvectors{popCount,1}.name = name;
                fseek(fid, offset, 'bof');
                nexFile.popvectors{popCount,1}.weights = fread(fid, [n 1], 'double');
                fseek(fid, filePosition, 'bof');
            case 5 % continuous variable
                contCount = contCount+1;
                nexFile.contvars{contCount,1}.name = name;
                nexFile.contvars{contCount,1}.ADFrequency = WFrequency;
                fseek(fid, offset, 'bof');
                nexFile.contvars{contCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
                nexFile.contvars{contCount,1}.fragmentStarts = fread(fid, [n 1], 'int32') + 1;
                nexFile.contvars{contCount,1}.data = fread(fid, [NPointsWave 1], 'int16').*ADtoMV + MVOfffset;
                fseek(fid, filePosition, 'bof');
            case 6 % marker
                markerCount = markerCount+1;
                nexFile.markers{markerCount,1}.name = name;
                fseek(fid, offset, 'bof');
                nexFile.markers{markerCount,1}.timestamps = fread(fid, [n 1], 'int32')./nexFile.freq;
                for i=1:NMarkers %#ok<FXSET>
                    nexFile.markers{markerCount,1}.values{i,1}.name = deblank(char(fread(fid, 64, 'char')'));
                    for p = 1:n
                        nexFile.markers{markerCount,1}.values{i,1}.strings{p, 1} = deblank(char(fread(fid, MarkerLength, 'char')'));
                    end
                end
                fseek(fid, filePosition, 'bof');
            otherwise
                disp (['unknown variable type ' num2str(type)]);
        end
        dummy = fread(fid, 60, 'char');
    end
    fclose(fid);

    % clear excess variables
    clear eventCount intervalCount waveCount popCount contCount markerCount
    clear type varVersion name offset n  WireNumber UnitNumber Gain Filter XPos YPos WFrequency
    clear ADtoMV NPointsWave NMarkers MarkerLength MVOfffset filePosition nvar p i magic

    % remove previous spike channel timestamp files
    disp('Removing previous spike channel timestamp files...')
    dirlist=dir([tempname,'*.mat']);
    try delete(dirlist.name); catch; end;
   
    % export spike channel timestamps
    disp('Exporting spike channel timestamps...');
    if isfield(nexFile,'neurons')==1,
        for nn=1:size(nexFile.neurons,1),
            disp(['...Loading spike channel (',nexFile.neurons{nn}.name,') data from ',tempname,'...']);
            neuronfname=[tempname,'-',nexFile.neurons{nn}.name];
            ts=nexFile.neurons{nn}.timestamps';
            freq=nexFile.freq;
            fullname=[hmiconfig.neurondir,neuronfname,'.mat'];
            %%% SAVE NEURON FILE
            save(fullname,'ts','freq','outputmatname')
        end
    else
        disp('...No spike channels found in this file.  Skipping this step')
        behav_only=1;
    end
    
    % export waveform templates and raw waveforms
    disp(' ')
    try
        figure; clf; cla; set(gcf,'Units','Normalized'); set(gcf,'Position',[0.1 0.3 0.8 0.5]); set(gca,'FontName','Arial')
        disp('Exporting waveform data ...');
        
        for nn=1:neuronCount, % scroll through each unit
            waverawname=[hmiconfig.wave_raw,tempname,'-',nexFile.waves{nn+neuronCount,1}.name,'_raw.mat'];
            wavetempname=[hmiconfig.wavetemp,tempname,'-',nexFile.waves{nn,1}.name,'.mat'];
            delete(waverawname);
            delete(wavetempname);
            disp(['...waveform data for unit ',nexFile.waves{nn+neuronCount,1}.name,'...'])
            try
                waverawdata=nexFile.waves{nn,1}.waveforms;
                wavetempdata=nexFile.waves{nn+neuronCount,1}.waveforms;
                subplot(1,neuronCount,nn)
                hold on
                for wf=1:20:size(waverawdata,2), plot(waverawdata(:,wf),'-','Color',[0.5 0.5 0.5],'LineWidth',0.01); end
                plot(wavetempdata,'r:','LineWidth',2)
                plot(mean(waverawdata'),'b-','LineWidth',4); %#ok<UDIM>
                ylabel('mV','FontSize',7); xlabel([nexFile.waves{nn+neuronCount,1}.name]); set(gca,'FontSize',7); axis square;
                if nn==1, title([tempname,' (',date,')'],'FontSize',12','FontWeight','Bold'); end
                save(waverawname,'waverawdata')
                save(wavetempname,'wavetempdata')
            catch
                disp('*** No waveform data found this unit ***')
            end
        end
        jpgfigname=[hmiconfig.wave_raw,filesep,'Figures',filesep,tempname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
    catch
    end
    % clear excess variables
    clear neuronCount nexFile neuronfname dummy ts freq nn matname outputmatname

    % export event codes (strobed bits)
    disp('Exporting event codes...');
    n = 0; nm = 0; nl = 0; ts = 0; m = 0; names = 0;
    nexfilename = char(filez(f));
    matfilename = nexfilename;
    fid = fopen(char(nexfname), 'r');
    if(fid == -1)
        disp('Cannot open file');
        return
    end
    varname='Strobed';
    magic = fread(fid, 1, 'int32');
    if magic ~= 827868494, error 'The file is not a valid .nex file'; end
    nexFile.version = fread(fid, 1, 'int32');
    nexFile.comment = deblank(char(fread(fid, 256, 'char')'));
    nexFile.freq = fread(fid, 1, 'double'); % frequency of collection of timestamped data
    nexFile.tbeg = fread(fid, 1, 'int32')./nexFile.freq; % trial beginning
    nexFile.tend = fread(fid, 1, 'int32')./nexFile.freq; % trial end -- duration = (tend-tbeg)/freq
    nvar = fread(fid, 1, 'int32'); % number of variables in the file
    fseek(fid, 260, 'cof');
    neuronCount=0; eventCount=0; intervalCount=0; waveCount=0; popCount=0; contCount=0; markerCount=0; name=zeros(1, 64); found=0; %#ok<*PREALL>
    for i=1:nvar
        type = fread(fid, 1, 'int32');      % Interpretation of type values: 0-neuron, 1-event, 2-interval, 3-waveform, 4-population vector, 5-continuous variable, 6-marker
        var_version = fread(fid, 1, 'int32');
        name = fread(fid, [1 64], 'char');  % [nvar 64] array of variable names
        offset = fread(fid, 1, 'int32');
        n = fread(fid, 1, 'int32');         % length of variable
        dummy = fread(fid, 32, 'char');
        adfreq = fread(fid, 1, 'double');
        adtomv = fread(fid, 1, 'double');
        npw = fread(fid, 1, 'int32');
        nm = fread(fid, 1, 'int32');        % number of fields in each marker
        nl = fread(fid, 1, 'int32');        % number of characters in each marker field
        dummy = fread(fid, 68, 'char');
        if type == 6 % marker
            found = 1;
            fseek(fid, offset, 'bof');
            ts = fread(fid, [1 n], 'int32');
            names = zeros(1,64);
            m = zeros(n,nl);
            names = fread(fid, [1 64], 'char');
            for j=1:n
                m(j,1:nl) = fread(fid, [1 nl], 'char');
            end
            Strobed=zeros(n,2);
            ts = ts/nexFile.freq;
            for k=1:n
                Strobed(k,1)=ts(k);
                Strobed(k,2)=(m(k,1)-48)*10000+(m(k,2)-48)*1000+(m(k,3)-48)*100+(m(k,4)-48)*10+m(k,5)-48;
            end
            break
        end
    end
    fclose(fid);
    clear ts m found dummy adfreq adtomv npw nm nl nexFile

    disp('Generate behaviour matrix...')
    %%% Identify indices for trial start and stops
    % start by saving whatever is pulled from the actual Strobed channel
    startpoints=find(Strobed(:,2)==hmiconfig.TRIAL_START); % identify all start points
    endpoints=find(Strobed(:,2)==hmiconfig.TRIAL_END); % identify all end points
    disp(' ')
    disp(['...Number of startpoints from NEX file: ',num2str(length(startpoints))])
    disp(['...Number of endpoints from NEX file: ',num2str(length(endpoints))])
    save([hmiconfig.startend,filesep,'RawValues',filesep,char(filez(f)),'_startpoints.txt'],'startpoints','-ascii','-tabs');
    save([hmiconfig.startend,filesep,'RawValues',filesep,char(filez(f)),'_endpoints.txt'],'endpoints','-ascii','-tabs');

    % either automatically generate startpoints or load manually edited ones
    if manual(f)==1,
        disp('...Loading manually generated start and stop points')
        if exist([hmiconfig.startend,char(filez(f)),'_startend.txt'],'file')==0,
            cd(hmiconfig.startend);
            [startendfname,startendpath]=uigetfile('*_startend.txt');
            startends=load([startendpath,startendfname]);
        else
            startends=load([hmiconfig.startend,char(filez(f)),'_startend.txt']);
        end
        startpoints=startends(:,1);
        endpoints=startends(:,2);
%     elseif length(startpoints)~=length(endpoints)
%         %autopts=inputdlg('The number of start and end points do not match.  Would you like to use the automated startpoints (Y/N)?')
%         autopts='y';
%         if ismember(autopts,{'Y','y','yes','Yes','1'})==1,
%             disp('...Automatically generating start and stop points')
%             startpoints=[1;endpoints+1]; startpoints(end)=[]; % remove last one
%         else
%             disp('...Attempting to correct...')
%             disp('......first pass...')
%             [startpoints,endpoints]=plx_checkpoints(startpoints,endpoints);
%             if min(endpoints-startpoints)<4,
%                 disp('...There may be a problem with the number of start and end points.  Attempting to fix...')
%                 for ii=1:length(startpoints),
%                     if endpoints(ii)-startpoints(ii)<4, % 4 is arbitrary minimum value of codes that should be in basecode
%                         endpoints(ii)=[];
%                         break
%                     end
%                 end
%             end
%             disp('......second pass...')
%             if length(startpoints)~=length(endpoints),
%                 disp('...The number of start and end points do not match.  Loading correction program...')
%                 [startpoints,endpoints]=plx_checkpoints(startpoints,endpoints);
%             end
%             if length(startpoints)~=length(endpoints),
%                 error('...Cannot reconcile number of start and end points.  Take a closer look.');
%             end
%             for ii=1:length(startpoints),
%                 if startpoints(ii)>endpoints(ii),
%                     error(['There is still a problem.  Check trial ',num2str(ii)])
%                 end
%             end
%         end
    elseif length(startpoints)~=length(endpoints)
        %% create manual list of start and end points based on endpoints
        startpoints=[];
        startpoints(1)=1;
        for st=2:length(endpoints),
            startpoints(st)=endpoints(st-1)+1; %#ok<AGROW>
        end
        startpoints=startpoints';
    end
    %% Figure for Startpoints
    figure % generate a figure to quickly confirm startpoints and endpoints
    clf; cla;
    set(gcf,'Units','Normalized');
    set(gcf,'Position',[0.1 0.3 0.8 0.5])
    set(gca,'FontName','Arial')
    plot(endpoints-startpoints); hold on; plot([0 length(startpoints)],[0 0],'r-','LineWidth',1.5);
    ylim([-10 50]); ylabel('Difference b/w Start/End Points','FontSize',7)
    title([matfilename,' - Trial Start and End points']); set(gca,'FontSize',7)
    text(10,-5,['Number of trials: ',num2str(length(startpoints))],'FontSize',9)

    %% Create behavioural matrix
    numtrials=length(startpoints); % total number of trials
    behav_matrix=nan(numtrials,51);
    % behav_matrix structure:
    %  1 Trial Number (assigned here, not by Cortex)            %  2 Start Code
    %  3 Paradigm Number                                        %  4 Condition Number
    %  5 ITImin                                                 %  6 ITImax
    %  7 ITIstep                                                %  8 penalty
    %  9 onFPmin                                                % 10 onFPmax
    % 11 onFPstep                                               % 12 FPmin
    % 13 FPmax                                                  % 14 FPstep
    % 15 CUEtime                                                % 16 CUEmin
    % 17 CUEmax                                                 % 18 CUEstep
    % 19 NOGOtime                                               % 20 GOtime
    % 21 REWint                                                 % 22 fixwin_x
    % 23 fixwin_y                                               % 30 trial outcome
    % 31 CTOA (choice on - cue on)                              % 32 SRT
    % 33 TRIAL DURATION (endtime-starttime)                     % 34 cue duration
    % 35 saccade duration                                       % 36 start position (x)
    % 37 end position (x)                                       % 38 net saccade
    % 39 REMARK (1=correct side, 2=incorrect side)              % 40 Start Trial
    % 41 End Trial                                              % 42 FPon
    % 43 FPoff                                                  % 44 CUEon
    % 45 CUEoff                                                 % 46 TARG1on
    % 47 TARG1off                                               % 48 TARG2on
    % 49 TARG2off                                               % 50 Eye Left Window
    % 51 Reward                                                 % 52 reward size
    % 53 Block Number (from Cortex)                             % 54 Trial
    % Number (from Cortex)
    h = waitbar(0,'Analyzing behavioural outcome of each trial...');  
    for t = 1:numtrials, % run on each trial
        waitbar(t/numtrials)
        trialevents = Strobed(startpoints(t):endpoints(t),:); % TRIAL_START
        %%% Since the number of hard-coded variables for each paradigm
        %%% (e.g., ITItime,etc.) changes from paradigm to paradigm, the
        %%% following code allows for differences.  However, it must be
        %%% customized for each paradigm

        %%% ADDING NEW PARADIGMS REQUIRES ADDITIONAL CODE HERE!
        parnum = trialevents(2,2);
        switch parnum
            case {400,401,402,403,411,412,413,421,422,423} % morph dms task
                try
                    trialevents_t = trialevents([1 23:end],:); % remove all early trial ID information
                    behav_matrix(t,2:23) = trialevents(1:22,2);
                catch
                    disp(['...trial #',num2str(t),' does not contain the proper number of codes.  Skipping this trial...'])
                    behav_matrix(t,1) = t; % paste trial number
                    behav_matrix(t,40) = min(trialevents(trialevents(:,2)==hmiconfig.TRIAL_START,1)); % Start Trial
                    behav_matrix(t,41) = trialevents(trialevents(:,2)==hmiconfig.TRIAL_END,1); % End Trial
                    continue
                end
            case {450,451,452,453,454,455,600,601,602,603,604,605,610} % morph dms+prior probability task (and associated training tasks)
                try
                    trialevents_t = trialevents([1 23:end],:); % remove all early trial ID information
                    behav_matrix(t,2:23) = trialevents(1:22,2);
                    behav_matrix(t,53)=trialevents(23,2); % paste block number
                    behav_matrix(t,54)=trialevents(24,2); % paste trial number
                catch
                    disp(['...trial #',num2str(t),' does not contain the proper number of codes.  Skipping this trial...'])
                    behav_matrix(t,1) = t; % paste trial number
                    behav_matrix(t,40) = min(trialevents(trialevents(:,2)==hmiconfig.TRIAL_START,1)); % Start Trial
                    behav_matrix(t,41) = trialevents(trialevents(:,2)==hmiconfig.TRIAL_END,1); % End Trial
                    behav_matrix(t,53)=0;
                    behav_matrix(t,54)=0; 
                    continue
                end
            case {500 530 550 570 561} % rsvp tasks
                try
                    trialevents_t = trialevents([1 10:end],:); % remove all early trial ID information
                    behav_matrix(t,2:10) = trialevents(1:9,2);
                catch
                    disp(['...trial #',num2str(t),' does not contain the proper number of codes.  Skipping this trial...'])
                    behav_matrix(t,1) = t; % paste trial number
                    behav_matrix(t,40) = min(trialevents(trialevents(:,2)==hmiconfig.TRIAL_START,1)); % Start Trial
                    behav_matrix(t,41) = trialevents(trialevents(:,2)==hmiconfig.TRIAL_END,1); % End Trial
                    continue
                end
            case {510 511 512} % morph reward rsvp task
                try
                    trialevents_t = trialevents([1 10:end],:); % remove all early trial ID information
                    behav_matrix(t,2:16) = trialevents(1:15,2);
                catch
                    disp(['...trial #',num2str(t),' does not contain the proper number of codes.  Skipping this trial...'])
                    behav_matrix(t,1) = t; % paste trial number
                    behav_matrix(t,40) = min(trialevents(trialevents(:,2)==hmiconfig.TRIAL_START,1)); % Start Trial
                    behav_matrix(t,41) = trialevents(trialevents(:,2)==hmiconfig.TRIAL_END,1); % End Trial
                    continue
                end
            case {40,050,50,999,2000} % rsvp search task
                try
                    trialevents_t = trialevents([1 12:end],:); % remove all early trial ID information
                    behav_matrix(t,[2:4 12:15 19]) = trialevents([1:3 5:9],2);
                catch
                    disp(['...trial #',num2str(t),' does not contain the proper number of codes.  Skipping this trial...'])
                    behav_matrix(t,1) = t; % paste trial number
                    behav_matrix(t,40) = min(trialevents(trialevents(:,2)==hmiconfig.TRIAL_START,1)); % Start Trial
                    behav_matrix(t,41) = trialevents(trialevents(:,2)==hmiconfig.TRIAL_END,1); % End Trial
                    continue
                end
            otherwise
                disp(['...trial #',num2str(t),': unknown trial type (parnum: ',num2str(parnum),')'])
                try
                    trialevents_t = trialevents([1 21:end],:); % remove all early trial ID information
                    behav_matrix(t,2:21) = trialevents(1:20,2);
                catch
                    disp(['...trial #',num2str(t),' does not contain the proper number of codes.  Skipping this trial...'])
                    behav_matrix(t,1) = t; % paste trial number
                    %behav_matrix(t,40) = min(trialevents(find(trialevents(:,2)==hmiconfig.TRIAL_START),1)); % Start Trial
                    %behav_matrix(t,41) = trialevents(find(trialevents(:,2)==hmiconfig.TRIAL_END),1); % End Trial
                    continue
                end
        end
        %%% include excess columns so that output matches with behaviour matrix generated in hmi400.m (cortex file analysis)
        behav_matrix(t,1) = t; % paste trial number
        behav_matrix(t,30) = hmi_sorttrialplx(trialevents_t(:,2),hmiconfig); % Trial outcome (see hmi_sorttrial320 for details)
        %%% paste indices for critical trial events
        for te = 1:length(trialevents_t(:,2)),
            switch trialevents_t(te,2) % if more than one instance of each code appears, takes last one
                case hmiconfig.TRIAL_START
                    behav_matrix(t,40) = min(trialevents_t(trialevents_t(:,2)==hmiconfig.TRIAL_START,1)); % Start Trial
                case hmiconfig.TRIAL_END
                    behav_matrix(t,41) = trialevents_t(trialevents_t(:,2)==hmiconfig.TRIAL_END,1); % End Trial
                case hmiconfig.TURN_FIXSPOT_ON
                    behav_matrix(t,42) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.TURN_FIXSPOT_ON, 1, 'last' ),1); % FPon
                case hmiconfig.TURN_FIXSPOT_OFF
                    behav_matrix(t,43) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.TURN_FIXSPOT_OFF, 1, 'last' ),1); % FPoff
                case hmiconfig.CUE_ON
                    behav_matrix(t,44) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.CUE_ON, 1, 'last' ),1); % CUEon
                case hmiconfig.CUE_OFF
                    behav_matrix(t,45) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.CUE_OFF, 1, 'last' ),1); % CUEoff
                case hmiconfig.TARG1_ON
                    behav_matrix(t,46) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.TARG1_ON, 1, 'last' ),1); % TARG1on
                case hmiconfig.TARG1_OFF
                    behav_matrix(t,47) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.TARG1_OFF, 1, 'last' ),1); % TARG1off
                case hmiconfig.TARG2_ON
                    behav_matrix(t,48) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.TARG2_ON, 1, 'last' ),1); % TARG2on
                case hmiconfig.TARG2_OFF
                    behav_matrix(t,49) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.TARG2_OFF, 1, 'last' ),1); % TARG2off
                case hmiconfig.EYE_LEFT_WINDOW
                    behav_matrix(t,50) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.EYE_LEFT_WINDOW, 1, 'last' ),1); % Eye Left Window
                case hmiconfig.REWARD
                    behav_matrix(t,51) = trialevents_t(find(trialevents_t(:,2)==hmiconfig.REWARD, 1, 'last' ),1); % Reward
                case hmiconfig.ENC_REWARDSML
                    behav_matrix(t,51) = 1; % Small Reward Given
                case hmiconfig.ENC_REWARDLRG
                    behav_matrix(t,51) = 2; % Large Reward Given
            end
        end
        %%% FILL IN BEHAVIOURAL MATRIX
        behav_matrix(t,33) = (behav_matrix(t,41) - behav_matrix(t,40)); % Trial Duration
        if behav_matrix(t,30)>1, % monkey acquired fixation
            behav_matrix(t,31) = (behav_matrix(t,46) - behav_matrix(t,44))*1000; % CTOA
        end
        if behav_matrix(t,30)>2, % monkey made a choice
            behav_matrix(t,32) =(behav_matrix(t,50) - behav_matrix(t,46))*1000; % SRT
            behav_matrix(t,34) = (behav_matrix(t,45) - behav_matrix(t,44))*1000; % Cue duration
            %behav_matrix(t,35) = x; % Saccade duration
        end
    end
    close(h);
    disp('Paradigm codes found: ');
    disp(unique(behav_matrix(:,3)))
    clear trialevents trialevents_t t startpoints endpoints parnum Strobed

    %% Export eye position signals
    disp('Export eye position channels...')
    try
        % AD15 = Heye, AD16 = Veye, these names must match the channel names in the PLX file
        [foundornot,adfreq, nHeye, Heye_ts, fn, Heye] = plx_readnexcont(char(nexfname),nexfilename,'eyeH'); %#ok<ASGLU>
        [foundornot,adfreq, nVeye, Veye_ts, fn, Veye] = plx_readnexcont(char(nexfname),nexfilename,'eyeV'); %#ok<ASGLU>
        Heye=Heye';
        Veye=Veye';
        %%% Create Eye Matrices (structure) (in ms)
        maxtrial=round(max(behav_matrix(:,33)))*1000; % determine maximum trial length
        eyesig=struct('Heye',nan(numtrials,maxtrial),'Veye',nan(numtrials,maxtrial)); % create empty structure for eye traces
        for t=1:numtrials, % run on each trial
            indtriallength=length(behav_matrix(t,40)*1000:behav_matrix(t,41)*1000);
            eyesig.Heye(t,1:indtriallength)=Heye(behav_matrix(t,40)*1000:behav_matrix(t,41)*1000);
            eyesig.Veye(t,1:indtriallength)=Veye(behav_matrix(t,40)*1000:behav_matrix(t,41)*1000);
        end
    catch
        disp('*** Unable to access eye trace information.  Using blank values. ***')
        Heye=zeros(1,1); Veye=zeros(1,1); Heye_ts=0; Veye_ts=0; adfreq=1000; eyesig.Heye=zeros(numtrials,1000); eyesig.Veye=zeros(numtrials,1000);
    end
    outputmatname=strcat(hmiconfig.outdir,tempname,'.mat');
    
    %% Remove previous matfile
    disp(' ')
    try
        disp(['Deleting ',outputmatname])
        delete(outputmatname);
    catch
        disp('No previous BEHAV file found')
    end
    disp('Saving compiled behavioural data...')
    save(outputmatname,'behav_matrix','Heye','Veye','Heye_ts','Veye_ts','adfreq','eyesig');
    clear Heye Veye foundornot adfreq nHeye nVeye Heye_ts Veye_ts fn eyesig behav_matrix
    disp(' ')
    % Remove previous LFP files
    disp('Removing previous LFP files...')
    killfiles=dir([hmiconfig.LFPdir,nexfilename,'-*.mat']);
    for kf=1:size(killfiles,1),
        disp(['...deleting ',killfiles(kf).name])
        delete([hmiconfig.LFPdir,killfiles(kf).name]);
    end
    disp('Exporting LFP channels...')
    for c = 1:length(channels),
        plx_readLFPdata(channels(c),nexfname,nexfilename,hmiconfig,outputmatname)
    end
    disp('Correcting LFP data...')
    if lfpopt(f)==1, plx_processLFPs(1,{nexfilename},1); % use old correction kernel
    else plx_processLFPs(1,{nexfilename},2); % use new correction kernel
    end
    if behav_only~=1,
        disp('Generating spike matrices...')
        plx_makespikemat(matfilename);
    end
    cprintf('*green',['Finished processing ', char(filez(f)),'.\n']);
    disp(' ')
end
return

function trial_outcome = hmi_sorttrialplx(trialdata,hmiconfig)
% This is the function that classifies trials as correct/incorrect/somewhat
% correct - this is where to change those definitions!!!!
% Potential Outcomes:
% 0 Default Value
% 1 Monkey never acquired the FP
% 2 Monkey broke fixation during trial
% 3 Incorrect choice (EASY)
% 4 Incorrect choice (HARD) (but rewarded)
% 5 Correct choice (EASY)
% 6 Correct choice (HARD)
trial_outcome = 0; % default to completely incorrect (safe bet :)
if isempty(find(trialdata==hmiconfig.ENC_NO_FIXATION, 1))==0, % never acquired FP
    trial_outcome = 1;
end
if isempty(find(trialdata==hmiconfig.ENC_BREAK_FIX_ERROR, 1))==0, % broke fixation
    trial_outcome = 2;
end
if isempty(find(trialdata==hmiconfig.ENC_RESP_WRONG_HARD, 1))==0||isempty(find(trialdata==hmiconfig.ENC_RESP_WRONG, 1))==0, % incorrect choice, hard
    trial_outcome = 3;
end
if isempty(find(trialdata==hmiconfig.ENC_RESP_WRONG_EASY, 1))==0, % incorrect choice, easy
    trial_outcome = 4;
end
if isempty(find(trialdata==hmiconfig.ENC_RESP_CORRECT_HARD, 1))==0, % correct choice, easy
    trial_outcome = 5;
end
if isempty(find(trialdata==hmiconfig.ENC_RESP_CORRECT_EASY, 1))==0||isempty(find(trialdata==hmiconfig.ENC_RESP_CORRECT, 1))==0, % correct choice, hard
    trial_outcome = 6;
end
return

function plx_readLFPdata(chan,nfnamelong,nfnameshort,hmiconfig,dataforLFP) %#ok<INUSD>
[foundornot,adfreq,n,initialtime,fn,LFPdata] = plx_readnexcont(char(nfnamelong),nfnameshort,char(chan)); %#ok<ASGLU>
if foundornot==1
    %   for i=1:length(LFP)
    %       LFP_times(i)=initialtime + (i-1)/adfreq;
    %   end
    LFPsignal=LFPdata';
    LFPfile=nfnameshort; %(1:length(nfnameshort)-4);
    LFPfilename=strcat(hmiconfig.LFPdir,LFPfile,'-',chan,'.mat');
    save(char(LFPfilename),'LFPsignal','initialtime','adfreq','fn','dataforLFP');
    clear LFPdata LFPsignal n fn adfreq
end
return