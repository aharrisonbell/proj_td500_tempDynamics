function rsvp_postNeuron_avgSpden(monkinitial,neuronlist);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% rsvp_postNeuron_avgSpden(files); %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% written by AHB, June 2010, copied verbatim from plx500_multispden
% Generates a single set of spike density functions for a list (neuronlist)
% of units.
% neuronlist (optional) = text file containing the units to include
% (default below)

%%% SETUP DEFAULTS
warning off;
hmiconfig=generate_hmi_configplex; % generates and loads config file
xrange=-200:400;
normwin=-100:0;

% March29,2009
if monkinitial=='S',
    monkeyname='Stewie'; sheetname='RSVP Cells_S';
    neuronlist(1).grids={'A7L2','A7L1','A6L3'}; % BodyPart Selective
    neuronlist(2).grids={'A6L2','A6L0','A5L2','A5L1','A5L0'}; % Face Selective
    neuronlist(3).grids={'A4L2','A4L1','A4R1'}; % No Category Selectivity
    neuronlist(4).grids={'A2L5','A0L6','A0L2','A0L0','P1L1','P2L3','P3L5','P3L4', 'P4L2','P4L4','P5L3','P6L3'}; % Object Selective
    neuronlist(5).grids={'P6L2','P6L1','P7L2'}; % Face Selective
elseif monkinitial=='W',
    monkeyname='Wiggum'; sheetname='RSVP Cells_W';
    neuronlist(1).grids={'A7R0','A7L1','A6L1'}; %  Places
    neuronlist(2).grids={'A0R0','A0R1','A1R0','A2R1','A3R0'}; % Faces
    neuronlist(3).grids={'P1R0','P1R3'}; % Bodyparts
    neuronlist(4).grids={'P3R0','P3R2'}; % Objects
    neuronlist(5).grids={'P7R0','P7R2'}; %  Face Selective
end

disp('*****************************************************************')
disp('* plx500_multispden.m - Analysis program for neuronal data from *')
disp('*   RSVP500 datafiles.  This program generates a single set of  *')
disp('*   spike density functions - the average of several units.     *')
disp('*****************************************************************')
disp(' ')
%%% LOAD FILES
disp('Loading unit names...')
[data,numgrids,counts_matrix,allunits,unit_index,unitdata]=plx500_prepproject1data(hmiconfig,sheetname);
for nl=1:size(neuronlist,2),
    disp('Finding units from the desired patch...')
    gridind=find(ismember(unit_index.GridLoc,neuronlist(nl).grids)==1 & strcmp(unit_index.SensoryConf,'Sensory')==1);
    disp(['..found ',num2str(length(gridind)),' units'])
    %%% LOAD GRAPH DATA AND PASTE INTO MASTER MATRIX
    disp('Loading spike density functions for each unit...')
    graphs_excite = struct('faces_raw',[],'fruit_raw',[],'places_raw',[],'bodyparts_raw',[],'objects_raw',[],...
        'faces_norm',[],'fruit_norm',[],'places_norm',[],'bodyparts_norm',[],'objects_norm',[],...
        'responses',[]);
    graphs_inhibit = struct('faces_raw',[],'fruit_raw',[],'places_raw',[],'bodyparts_raw',[],'objects_raw',[],...
        'faces_norm',[],'fruit_norm',[],'places_norm',[],'bodyparts_norm',[],'objects_norm',[],...
        'responses',[]);
    graphs_both = struct('faces_raw',[],'fruit_raw',[],'places_raw',[],'bodyparts_raw',[],'objects_raw',[],...
        'faces_norm',[],'fruit_norm',[],'places_norm',[],'bodyparts_norm',[],'objects_norm',[],...
        'responses',[]);
    graphs_all = struct('faces_raw',[],'fruit_raw',[],'places_raw',[],'bodyparts_raw',[],'objects_raw',[],...
        'faces_norm',[],'fruit_norm',[],'places_norm',[],'bodyparts_norm',[],'objects_norm',[],...
        'responses',[]);
    for un=1:length(gridind),
        temp1=char(unit_index.PlxFile(gridind(un)));
        unitname=[temp1(1:12),'-',char(unit_index.UnitName(gridind(un)))];
        disp(['...loading data from ',unitname]);
        load([hmiconfig.rsvp500spks,unitname,'-500responsedata.mat']); % load unit data

        %%% Determine response type
        if strcmp(respstructsingle.conf_excite,'Non-Responsive')==1,
            disp('... Neuron is non-responsive.  Skipping...')
            continue
        end
        if ismember(respstructsingle.conf_excite,{'Excite'})==1,
            load([hmiconfig.rsvp500spks,unitname,'-500graphdata.mat']); % load spike density data
            xwindow=1000+xrange(1):1000+xrange(end);
            xnormwin=1000+normwin(1):1000+normwin(end);
            graphs_excite.faces_raw=[graphs_excite.faces_raw;graphstructsingle.faces_avg(xwindow)];
            graphs_excite.fruit_raw=[graphs_excite.fruit_raw;graphstructsingle.fruit_avg(xwindow)];
            graphs_excite.places_raw=[graphs_excite.places_raw;graphstructsingle.places_avg(xwindow)];
            graphs_excite.bodyparts_raw=[graphs_excite.bodyparts_raw;graphstructsingle.bodyp_avg(xwindow)];
            graphs_excite.objects_raw=[graphs_excite.objects_raw;graphstructsingle.objct_avg(xwindow)];

            % normalize to max of each response
            clip=graphstructsingle.allconds_avg(:,xwindow);
            normalizer=max(max(clip));
            graphs_excite.faces_norm=[graphs_excite.faces_norm;graphstructsingle.faces_avg(xwindow)/normalizer];
            graphs_excite.fruit_norm=[graphs_excite.fruit_norm;graphstructsingle.fruit_avg(xwindow)/normalizer];
            graphs_excite.places_norm=[graphs_excite.places_norm;graphstructsingle.places_avg(xwindow)/normalizer];
            graphs_excite.bodyparts_norm=[graphs_excite.bodyparts_norm;graphstructsingle.bodyp_avg(xwindow)/normalizer];
            graphs_excite.objects_norm=[graphs_excite.objects_norm;graphstructsingle.objct_avg(xwindow)/normalizer];

            % population responses
            graphs_excite.responses(un,1:100)=respstructsingle.m_epoch1;
        end
        if ismember(respstructsingle.conf_excite,{'Inhibit'})==1,
            load([hmiconfig.rsvp500spks,unitname,'-500graphdata.mat']); % load spike density data
            xwindow=1000+xrange(1):1000+xrange(end);
            xnormwin=1000+normwin(1):1000+normwin(end);
            graphs_inhibit.faces_raw=[graphs_inhibit.faces_raw;graphstructsingle.faces_avg(xwindow)];
            graphs_inhibit.fruit_raw=[graphs_inhibit.fruit_raw;graphstructsingle.fruit_avg(xwindow)];
            graphs_inhibit.places_raw=[graphs_inhibit.places_raw;graphstructsingle.places_avg(xwindow)];
            graphs_inhibit.bodyparts_raw=[graphs_inhibit.bodyparts_raw;graphstructsingle.bodyp_avg(xwindow)];
            graphs_inhibit.objects_raw=[graphs_inhibit.objects_raw;graphstructsingle.objct_avg(xwindow)];

            % normalize to max of each response
            clip=graphstructsingle.allconds_avg(:,xwindow);
            normalizer=max(max(clip));
            graphs_inhibit.faces_norm=[graphs_inhibit.faces_norm;graphstructsingle.faces_avg(xwindow)/normalizer];
            graphs_inhibit.fruit_norm=[graphs_inhibit.fruit_norm;graphstructsingle.fruit_avg(xwindow)/normalizer];
            graphs_inhibit.places_norm=[graphs_inhibit.places_norm;graphstructsingle.places_avg(xwindow)/normalizer];
            graphs_inhibit.bodyparts_norm=[graphs_inhibit.bodyparts_norm;graphstructsingle.bodyp_avg(xwindow)/normalizer];
            graphs_inhibit.objects_norm=[graphs_inhibit.objects_norm;graphstructsingle.objct_avg(xwindow)/normalizer];

            % population responses
            graphs_inhibit.responses(un,1:100)=respstructsingle.m_epoch1;

        end
        if ismember(respstructsingle.conf_excite,{'Both'})==1,
            load([hmiconfig.rsvp500spks,unitname,'-500graphdata.mat']); % load spike density data
            xwindow=1000+xrange(1):1000+xrange(end);
            xnormwin=1000+normwin(1):1000+normwin(end);
            graphs_both.faces_raw=[graphs_both.faces_raw;graphstructsingle.faces_avg(xwindow)];
            graphs_both.fruit_raw=[graphs_both.fruit_raw;graphstructsingle.fruit_avg(xwindow)];
            graphs_both.places_raw=[graphs_both.places_raw;graphstructsingle.places_avg(xwindow)];
            graphs_both.bodyparts_raw=[graphs_both.bodyparts_raw;graphstructsingle.bodyp_avg(xwindow)];
            graphs_both.objects_raw=[graphs_both.objects_raw;graphstructsingle.objct_avg(xwindow)];

            % normalize to max of each response
            clip=graphstructsingle.allconds_avg(:,xwindow);
            normalizer=max(max(clip));
            graphs_both.faces_norm=[graphs_both.faces_norm;graphstructsingle.faces_avg(xwindow)/normalizer];
            graphs_both.fruit_norm=[graphs_both.fruit_norm;graphstructsingle.fruit_avg(xwindow)/normalizer];
            graphs_both.places_norm=[graphs_both.places_norm;graphstructsingle.places_avg(xwindow)/normalizer];
            graphs_both.bodyparts_norm=[graphs_both.bodyparts_norm;graphstructsingle.bodyp_avg(xwindow)/normalizer];
            graphs_both.objects_norm=[graphs_both.objects_norm;graphstructsingle.objct_avg(xwindow)/normalizer];

            % population responses
            graphs_both.responses(un,1:100)=respstructsingle.m_epoch1;

        end
        if ismember(respstructsingle.conf_excite,{'Excite','Inhibit','Both'})==1,
            load([hmiconfig.rsvp500spks,unitname,'-500graphdata.mat']); % load spike density data
            xwindow=1000+xrange(1):1000+xrange(end);
            xnormwin=1000+normwin(1):1000+normwin(end);
            graphs_all.faces_raw=[graphs_all.faces_raw;graphstructsingle.faces_avg(xwindow)];
            graphs_all.fruit_raw=[graphs_all.fruit_raw;graphstructsingle.fruit_avg(xwindow)];
            graphs_all.places_raw=[graphs_all.places_raw;graphstructsingle.places_avg(xwindow)];
            graphs_all.bodyparts_raw=[graphs_all.bodyparts_raw;graphstructsingle.bodyp_avg(xwindow)];
            graphs_all.objects_raw=[graphs_all.objects_raw;graphstructsingle.objct_avg(xwindow)];


            % normalize to max of each response
            clip=graphstructsingle.allconds_avg(:,xwindow);
            normalizer=max(max(clip));
            graphs_all.faces_norm=[graphs_all.faces_norm;graphstructsingle.faces_avg(xwindow)/normalizer];
            graphs_all.fruit_norm=[graphs_all.fruit_norm;graphstructsingle.fruit_avg(xwindow)/normalizer];
            graphs_all.places_norm=[graphs_all.places_norm;graphstructsingle.places_avg(xwindow)/normalizer];
            graphs_all.bodyparts_norm=[graphs_all.bodyparts_norm;graphstructsingle.bodyp_avg(xwindow)/normalizer];
            graphs_all.objects_norm=[graphs_all.objects_norm;graphstructsingle.objct_avg(xwindow)/normalizer];

            % population responses
            graphs_all.responses(un,1:100)=respstructsingle.m_epoch1;

        end
        clear respstructsingle graphstructsingle
    end

    %%% GENERATE THE FIGURE
    figure; clf; cla; % colour map showing category selectivity/proportion
    set(gcf,'Units','Normalized','Position',[0.05 0.1 0.9 0.8])
    set(gca,'FontName','Arial','FontSize',8)
    subplot(4,3,1) % raw
    hold on
    plot(xrange,mean(graphs_excite.faces_raw),'r-','LineWidth',1.5)
    plot(xrange,mean(graphs_excite.fruit_raw),'m-','LineWidth',1.5)
    plot(xrange,mean(graphs_excite.places_raw),'b-','LineWidth',1.5)
    plot(xrange,mean(graphs_excite.bodyparts_raw),'y-','LineWidth',1.5)
    plot(xrange,mean(graphs_excite.objects_raw),'g-','LineWidth',1.5)
    % error bars
    plot(xrange,mean(graphs_excite.faces_raw)+sem(graphs_excite.faces_raw),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.fruit_raw)+sem(graphs_excite.fruit_raw),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.places_raw)+sem(graphs_excite.places_raw),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.bodyparts_raw)+sem(graphs_excite.bodyparts_raw),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.objects_raw)+sem(graphs_excite.objects_raw),'g-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.faces_raw)-sem(graphs_excite.faces_raw),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.fruit_raw)-sem(graphs_excite.fruit_raw),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.places_raw)-sem(graphs_excite.places_raw),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.bodyparts_raw)-sem(graphs_excite.bodyparts_raw),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.objects_raw)-sem(graphs_excite.objects_raw),'g-','LineWidth',0.5)
    xlabel('Time from stimulus onset (ms)','FontSize',8); xlim([xrange(1) xrange(end)]);
    ylabel('Firing rate (sp/s)','FontSize',8); set(gca,'FontSize',7);
    title([monkeyname,' Patch#',num2str(nl),' Firing Rate (n=',num2str(size(graphs_excite.faces_raw,1)),')'],'FontSize',9,'FontWeight','Bold')

    subplot(4,3,2) % norm
    hold on
    plot(xrange,mean(graphs_excite.faces_norm),'r-','LineWidth',1.5)
    plot(xrange,mean(graphs_excite.fruit_norm),'m-','LineWidth',1.5)
    plot(xrange,mean(graphs_excite.places_norm),'b-','LineWidth',1.5)
    plot(xrange,mean(graphs_excite.bodyparts_norm),'y-','LineWidth',1.5)
    plot(xrange,mean(graphs_excite.objects_norm),'g-','LineWidth',1.5)
    % error bars
    plot(xrange,mean(graphs_excite.faces_norm)+sem(graphs_excite.faces_norm),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.fruit_norm)+sem(graphs_excite.fruit_norm),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.places_norm)+sem(graphs_excite.places_norm),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.bodyparts_norm)+sem(graphs_excite.bodyparts_norm),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.objects_norm)+sem(graphs_excite.objects_norm),'g-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.faces_norm)-sem(graphs_excite.faces_norm),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.fruit_norm)-sem(graphs_excite.fruit_norm),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.places_norm)-sem(graphs_excite.places_norm),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.bodyparts_norm)-sem(graphs_excite.bodyparts_norm),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_excite.objects_norm)-sem(graphs_excite.objects_norm),'g-','LineWidth',0.5)
    xlabel('Time from stimulus onset (ms)','FontSize',8); xlim([xrange(1) xrange(end)]);
    ylabel('Normalized firing rate (sp/s)','FontSize',8); set(gca,'FontSize',7);
    title('Normalized Firing Rate','FontSize',9,'FontWeight','Bold')

    subplot(4,3,3) % responses
    hold on
    pcolor(1:100,1:size(graphs_excite.responses,1),graphs_excite.responses)
    shading flat; % colorbar('SouthOutside');
    plot([20 20],[0 size(graphs_excite.responses,1)],'k-','LineWidth',1)
    plot([40 40],[0 size(graphs_excite.responses,1)],'k-','LineWidth',1)
    plot([60 60],[0 size(graphs_excite.responses,1)],'k-','LineWidth',1)
    plot([80 80],[0 size(graphs_excite.responses,1)],'k-','LineWidth',1)
    set(gca,'FontSize',7); box off; axis ij; ylim([0 size(graphs_excite.responses,1)]);
    xlim([1 100]); set(gca,'XTick',[10 30 50 70 90],'XTickLabel',{'Faces','Fruit','Places','Bodyparts','Objects'})
    set(gca,'YTick',[1 size(graphs_excite.responses,1)])
    ylabel('Neuron Number','FontSize',7)
    xlabel('Stimulus Identity','FontSize',1)

    subplot(4,3,4) % raw
    hold on
    plot(xrange,mean(graphs_inhibit.faces_raw),'r-','LineWidth',1.5)
    plot(xrange,mean(graphs_inhibit.fruit_raw),'m-','LineWidth',1.5)
    plot(xrange,mean(graphs_inhibit.places_raw),'b-','LineWidth',1.5)
    plot(xrange,mean(graphs_inhibit.bodyparts_raw),'y-','LineWidth',1.5)
    plot(xrange,mean(graphs_inhibit.objects_raw),'g-','LineWidth',1.5)
    % error bars
    plot(xrange,mean(graphs_inhibit.faces_raw)+sem(graphs_inhibit.faces_raw),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.fruit_raw)+sem(graphs_inhibit.fruit_raw),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.places_raw)+sem(graphs_inhibit.places_raw),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.bodyparts_raw)+sem(graphs_inhibit.bodyparts_raw),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.objects_raw)+sem(graphs_inhibit.objects_raw),'g-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.faces_raw)-sem(graphs_inhibit.faces_raw),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.fruit_raw)-sem(graphs_inhibit.fruit_raw),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.places_raw)-sem(graphs_inhibit.places_raw),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.bodyparts_raw)-sem(graphs_inhibit.bodyparts_raw),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.objects_raw)-sem(graphs_inhibit.objects_raw),'g-','LineWidth',0.5)
    xlabel('Time from stimulus onset (ms)','FontSize',8); xlim([xrange(1) xrange(end)]);
    ylabel('Firing rate (sp/s)','FontSize',8); set(gca,'FontSize',7);
    title([monkeyname,' Patch#',num2str(nl),' Firing Rate (n=',num2str(size(graphs_inhibit.faces_raw,1)),')'],'FontSize',9,'FontWeight','Bold')

    subplot(4,3,5) % norm
    hold on
    plot(xrange,mean(graphs_inhibit.faces_norm),'r-','LineWidth',1.5)
    plot(xrange,mean(graphs_inhibit.fruit_norm),'m-','LineWidth',1.5)
    plot(xrange,mean(graphs_inhibit.places_norm),'b-','LineWidth',1.5)
    plot(xrange,mean(graphs_inhibit.bodyparts_norm),'y-','LineWidth',1.5)
    plot(xrange,mean(graphs_inhibit.objects_norm),'g-','LineWidth',1.5)
    % error bars
    plot(xrange,mean(graphs_inhibit.faces_norm)+sem(graphs_inhibit.faces_norm),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.fruit_norm)+sem(graphs_inhibit.fruit_norm),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.places_norm)+sem(graphs_inhibit.places_norm),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.bodyparts_norm)+sem(graphs_inhibit.bodyparts_norm),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.objects_norm)+sem(graphs_inhibit.objects_norm),'g-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.faces_norm)-sem(graphs_inhibit.faces_norm),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.fruit_norm)-sem(graphs_inhibit.fruit_norm),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.places_norm)-sem(graphs_inhibit.places_norm),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.bodyparts_norm)-sem(graphs_inhibit.bodyparts_norm),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_inhibit.objects_norm)-sem(graphs_inhibit.objects_norm),'g-','LineWidth',0.5)
    xlabel('Time from stimulus onset (ms)','FontSize',8); xlim([xrange(1) xrange(end)]);
    ylabel('Normalized firing rate (sp/s)','FontSize',8); set(gca,'FontSize',7);
    title('Normalized Firing Rate','FontSize',9,'FontWeight','Bold')

    subplot(4,3,6) % responses
    hold on
    pcolor(1:100,1:size(graphs_inhibit.responses,1),graphs_inhibit.responses)
    shading flat; % colorbar('SouthOutside');
    plot([20 20],[0 size(graphs_inhibit.responses,1)],'k-','LineWidth',1)
    plot([40 40],[0 size(graphs_inhibit.responses,1)],'k-','LineWidth',1)
    plot([60 60],[0 size(graphs_inhibit.responses,1)],'k-','LineWidth',1)
    plot([80 80],[0 size(graphs_inhibit.responses,1)],'k-','LineWidth',1)
    set(gca,'FontSize',7); box off; axis ij; ylim([0 size(graphs_inhibit.responses,1)]);
    xlim([1 100]); set(gca,'XTick',[10 30 50 70 90],'XTickLabel',{'Faces','Fruit','Places','Bodyparts','Objects'})
    set(gca,'YTick',[1 size(graphs_inhibit.responses,1)])
    ylabel('Neuron Number','FontSize',7)
    xlabel('Stimulus Identity','FontSize',1)

    subplot(4,3,7) % raw
    hold on
    plot(xrange,mean(graphs_both.faces_raw),'r-','LineWidth',1.5)
    plot(xrange,mean(graphs_both.fruit_raw),'m-','LineWidth',1.5)
    plot(xrange,mean(graphs_both.places_raw),'b-','LineWidth',1.5)
    plot(xrange,mean(graphs_both.bodyparts_raw),'y-','LineWidth',1.5)
    plot(xrange,mean(graphs_both.objects_raw),'g-','LineWidth',1.5)
    % error bars
    plot(xrange,mean(graphs_both.faces_raw)+sem(graphs_both.faces_raw),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.fruit_raw)+sem(graphs_both.fruit_raw),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.places_raw)+sem(graphs_both.places_raw),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.bodyparts_raw)+sem(graphs_both.bodyparts_raw),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.objects_raw)+sem(graphs_both.objects_raw),'g-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.faces_raw)-sem(graphs_both.faces_raw),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.fruit_raw)-sem(graphs_both.fruit_raw),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.places_raw)-sem(graphs_both.places_raw),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.bodyparts_raw)-sem(graphs_both.bodyparts_raw),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.objects_raw)-sem(graphs_both.objects_raw),'g-','LineWidth',0.5)
    xlabel('Time from stimulus onset (ms)','FontSize',8); xlim([xrange(1) xrange(end)]);
    ylabel('Firing rate (sp/s)','FontSize',8); set(gca,'FontSize',7);
    title([monkeyname,' Patch#',num2str(nl),' Firing Rate (n=',num2str(size(graphs_both.faces_raw,1)),')'],'FontSize',9,'FontWeight','Bold')

    subplot(4,3,8) % norm
    hold on
    plot(xrange,mean(graphs_both.faces_norm),'r-','LineWidth',1.5)
    plot(xrange,mean(graphs_both.fruit_norm),'m-','LineWidth',1.5)
    plot(xrange,mean(graphs_both.places_norm),'b-','LineWidth',1.5)
    plot(xrange,mean(graphs_both.bodyparts_norm),'y-','LineWidth',1.5)
    plot(xrange,mean(graphs_both.objects_norm),'g-','LineWidth',1.5)
    % error bars
    plot(xrange,mean(graphs_both.faces_norm)+sem(graphs_both.faces_norm),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.fruit_norm)+sem(graphs_both.fruit_norm),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.places_norm)+sem(graphs_both.places_norm),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.bodyparts_norm)+sem(graphs_both.bodyparts_norm),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.objects_norm)+sem(graphs_both.objects_norm),'g-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.faces_norm)-sem(graphs_both.faces_norm),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.fruit_norm)-sem(graphs_both.fruit_norm),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.places_norm)-sem(graphs_both.places_norm),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.bodyparts_norm)-sem(graphs_both.bodyparts_norm),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_both.objects_norm)-sem(graphs_both.objects_norm),'g-','LineWidth',0.5)
    xlabel('Time from stimulus onset (ms)','FontSize',8); xlim([xrange(1) xrange(end)]);
    ylabel('Normalized firing rate (sp/s)','FontSize',8); set(gca,'FontSize',7);
    title('Normalized Firing Rate','FontSize',9,'FontWeight','Bold')

    subplot(4,3,9) % responses
    hold on
    pcolor(1:100,1:size(graphs_both.responses,1),graphs_both.responses)
    shading flat; % colorbar('SouthOutside');
    plot([20 20],[0 size(graphs_both.responses,1)],'k-','LineWidth',1)
    plot([40 40],[0 size(graphs_both.responses,1)],'k-','LineWidth',1)
    plot([60 60],[0 size(graphs_both.responses,1)],'k-','LineWidth',1)
    plot([80 80],[0 size(graphs_both.responses,1)],'k-','LineWidth',1)
    set(gca,'FontSize',7); box off; axis ij; ylim([0 size(graphs_both.responses,1)]);
    xlim([1 100]); set(gca,'XTick',[10 30 50 70 90],'XTickLabel',{'Faces','Fruit','Places','Bodyparts','Objects'})
    set(gca,'YTick',[1 size(graphs_both.responses,1)])
    ylabel('Neuron Number','FontSize',7)
    xlabel('Stimulus Identity','FontSize',1)


    subplot(4,3,10) % raw
    hold on
    plot(xrange,mean(graphs_all.faces_raw),'r-','LineWidth',1.5)
    plot(xrange,mean(graphs_all.fruit_raw),'m-','LineWidth',1.5)
    plot(xrange,mean(graphs_all.places_raw),'b-','LineWidth',1.5)
    plot(xrange,mean(graphs_all.bodyparts_raw),'y-','LineWidth',1.5)
    plot(xrange,mean(graphs_all.objects_raw),'g-','LineWidth',1.5)
    % error bars
    plot(xrange,mean(graphs_all.faces_raw)+sem(graphs_all.faces_raw),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.fruit_raw)+sem(graphs_all.fruit_raw),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.places_raw)+sem(graphs_all.places_raw),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.bodyparts_raw)+sem(graphs_all.bodyparts_raw),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.objects_raw)+sem(graphs_all.objects_raw),'g-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.faces_raw)-sem(graphs_all.faces_raw),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.fruit_raw)-sem(graphs_all.fruit_raw),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.places_raw)-sem(graphs_all.places_raw),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.bodyparts_raw)-sem(graphs_all.bodyparts_raw),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.objects_raw)-sem(graphs_all.objects_raw),'g-','LineWidth',0.5)
    xlabel('Time from stimulus onset (ms)','FontSize',8); xlim([xrange(1) xrange(end)]);
    ylabel('Firing rate (sp/s)','FontSize',8); set(gca,'FontSize',7);
    title([monkeyname,' Patch#',num2str(nl),' Firing Rate (n=',num2str(size(graphs_all.faces_raw,1)),')'],'FontSize',9,'FontWeight','Bold')

    subplot(4,3,11) % norm
    hold on
    plot(xrange,mean(graphs_all.faces_norm),'r-','LineWidth',1.5)
    plot(xrange,mean(graphs_all.fruit_norm),'m-','LineWidth',1.5)
    plot(xrange,mean(graphs_all.places_norm),'b-','LineWidth',1.5)
    plot(xrange,mean(graphs_all.bodyparts_norm),'y-','LineWidth',1.5)
    plot(xrange,mean(graphs_all.objects_norm),'g-','LineWidth',1.5)
    % error bars
    plot(xrange,mean(graphs_all.faces_norm)+sem(graphs_all.faces_norm),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.fruit_norm)+sem(graphs_all.fruit_norm),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.places_norm)+sem(graphs_all.places_norm),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.bodyparts_norm)+sem(graphs_all.bodyparts_norm),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.objects_norm)+sem(graphs_all.objects_norm),'g-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.faces_norm)-sem(graphs_all.faces_norm),'r-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.fruit_norm)-sem(graphs_all.fruit_norm),'m-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.places_norm)-sem(graphs_all.places_norm),'b-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.bodyparts_norm)-sem(graphs_all.bodyparts_norm),'y-','LineWidth',0.5)
    plot(xrange,mean(graphs_all.objects_norm)-sem(graphs_all.objects_norm),'g-','LineWidth',0.5)
    xlabel('Time from stimulus onset (ms)','FontSize',8); xlim([xrange(1) xrange(end)]);
    ylabel('Normalized firing rate (sp/s)','FontSize',8); set(gca,'FontSize',7);
    title('Normalized Firing Rate','FontSize',9,'FontWeight','Bold')

    subplot(4,3,12) % responses
    hold on
    pcolor(1:100,1:size(graphs_all.responses,1),graphs_all.responses)
    shading flat; % colorbar('SouthOutside');
    plot([20 20],[0 size(graphs_all.responses,1)],'k-','LineWidth',1)
    plot([40 40],[0 size(graphs_all.responses,1)],'k-','LineWidth',1)
    plot([60 60],[0 size(graphs_all.responses,1)],'k-','LineWidth',1)
    plot([80 80],[0 size(graphs_all.responses,1)],'k-','LineWidth',1)
    set(gca,'FontSize',7); box off; axis ij; ylim([0 size(graphs_all.responses,1)]);
    xlim([1 100]); set(gca,'XTick',[10 30 50 70 90],'XTickLabel',{'Faces','Fruit','Places','Bodyparts','Objects'})
    set(gca,'YTick',[1 size(graphs_all.responses,1)])
    ylabel('Neuron Number','FontSize',7)
    xlabel('Stimulus Identity','FontSize',1)

    jpgfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_AvgSpden_Patch',num2str(nl),'_',monkeyname,'.jpg']; print(gcf,jpgfigname,'-djpeg') % generates an JPEG file of the figure
    illfigname=[hmiconfig.rootdir,'rsvp500_postNeuron',filesep,'RSVP_AvgSpden_Patch',num2str(nl),'_',monkeyname,'.ai']; print(gcf,illfigname,'-dill') % generates an Adobe Illustrator file of the figure
    if hmiconfig.printer==1, print; end % prints the figure to the default printer (if printer==1)
end
