function td500_latency
global lsnconfig NeuronData qualityPointer
disp('Figure 1 - Latency')

% average latency, histogram, etc.
figure; clf; cla; %
set(gcf,'Units','Normalized','Position',[0.05 0.15 0.5 0.5]); set(gca,'FontName','Arial','FontSize',10);
hold on
pointer=intersect(qualityPointer,find(strcmp(NeuronData.pref_excite','Faces')==1));
[avgfunc, semfunc]=td500_avgSpden(pointer,3);
plx_shade_spden(-100:500,avgfunc(1,:),semfunc(1,:),0.1,'r') %
epsp_latencies=td500_EPSP_latency(pointer);
plot([nanmean(epsp_latencies(:,1),1) nanmean(epsp_latencies(:,1),1)],[5 40],'r-')

pointer=intersect(qualityPointer,find(strcmp(NeuronData.pref_excite','Fruit')==1));
[avgfunc, semfunc]=td500_avgSpden(pointer,3);
plx_shade_spden(-100:500,avgfunc(2,:),semfunc(2,:),0.1,'m') %
epsp_latencies=td500_EPSP_latency(pointer);
plot([nanmean(epsp_latencies(:,2),1) nanmean(epsp_latencies(:,2),1)],[5 40],'m-')

pointer=intersect(qualityPointer,find(strcmp(NeuronData.pref_excite','BodyParts')==1));
[avgfunc, semfunc]=td500_avgSpden(pointer,3);
plx_shade_spden(-100:500,avgfunc(3,:),semfunc(3,:),0.1,'c') %
epsp_latencies=td500_EPSP_latency(pointer);
plot([nanmean(epsp_latencies(:,3),1) nanmean(epsp_latencies(:,3),1)],[5 40],'c-')

pointer=intersect(qualityPointer,find(strcmp(NeuronData.pref_excite','Places')==1));
[avgfunc, semfunc]=td500_avgSpden(pointer,3);
plx_shade_spden(-100:500,avgfunc(4,:),semfunc(4,:),0.1,'b') %
epsp_latencies=td500_EPSP_latency(pointer);
plot([nanmean(epsp_latencies(:,4),1) nanmean(epsp_latencies(:,4),1)],[5 40],'b-')

pointer=intersect(qualityPointer,find(strcmp(NeuronData.pref_excite','Objects')==1));
[avgfunc, semfunc]=td500_avgSpden(pointer,3);
plx_shade_spden(-100:500,avgfunc(5,:),semfunc(5,:),0.1,'g') %
epsp_latencies=td500_EPSP_latency(pointer);
plot([nanmean(epsp_latencies(:,5),1) nanmean(epsp_latencies(:,5),1)],[5 40],'g-')


xdata=NeuronData.APindexNumber(qualityPointer);
ap=5;
pointer=intersect(qualityPointer,find(strcmp(NeuronData.preferred_category','Faces')==1 & ismember(NeuronData.APindexNumber',[ap ap+1 ap+2])==1));



subplot(2,2,1); hold on % mean latency in anterior vs. posterior locations
for ap = 5:3:19,
    pointer=intersect(qualityPointer,find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Faces')==1 & ismember(NeuronData.APindexNumber',[ap ap+1 ap+2])==1));
    [avgfunc, semfunc]=td500_avgSpden(pointer,3);
    plx_shade_spden(-100:500,avgfunc(1,:),semfunc(1,:),0.1,'k') % 
    xlim([-100 500]);
    mean_latency=mean(NeuronData.cat_latency(pointer,1)); 
    plot([mean_latency mean_latency],[0 40],'k:')

end
    
subplot(2,2,2); hold on % mean latency in anterior vs. posterior locations
for ap = 5:3:19,
    pointer=intersect(qualityPointer,find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Faces')==1 & ismember(NeuronData.APindexNumber',[ap ap+1 ap+2])==1));
    [avgfunc, semfunc]=td500_avgSpden(pointer,3);
    plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(end),avgfunc(5,:),semfunc(5,:),0.1,'g') % 
    xlim([-100 500]);
    mean_latency=mean(NeuronData.cat_latency(pointer,1));    
    plot([mean_latency mean_latency],[0 40],'k:')
end
    
subplot(2,2,3); hold on % mean latency in anterior vs. posterior locations
for ap = 5:3:19,
    pointer=intersect(qualityPointer,find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Places')==1 & ismember(NeuronData.APindexNumber',[ap ap+1 ap+2])==1));
    [avgfunc, semfunc]=td500_avgSpden(pointer,3);
    plx_shade_spden(-100:500,avgfunc(1,:),semfunc(1,:),0.1,'k') % 
    xlim([-100 500]);
    mean_latency=mean(NeuronData.cat_latency(pointer,1)); 
    plot([mean_latency mean_latency],[0 40],'k:')

end
    
subplot(2,2,4); hold on % mean latency in anterior vs. posterior locations
for ap = 5:3:19,
    pointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Objects')==1 & ismember(NeuronData.APindexNumber',[ap ap+1 ap+2])==1);
    [avgfunc, semfunc]=td500_avgSpden(pointer,3);
    plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(end),avgfunc(5,:),semfunc(5,:),0.1,'g') % 
    xlim([-100 500]);
    mean_latency=mean(NeuronData.cat_latency(pointer,1));    
    plot([mean_latency mean_latency],[0 40],'k:')
end






disp('Figure 2A,B - Effect of Facial Expression on Face Responses')

figure; clf; cla; %
set(gcf,'Units','Normalized','Position',[0.05 0.15 0.5 0.5]); set(gca,'FontName','Arial','FontSize',10);

subplot(2,2,1); hold on % faces
TEpointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Faces')==1 & NeuronData.APindexNumber' < 10);
[avgfunc, semfunc]=td500_avgSpden(NeuronData,TEpointer,lsnconfig);
plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(2),avgfunc(1,:),semfunc(1,:),0.25,'r') % 
TEOpointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Faces')==1 & NeuronData.APindexNumber' > 16);
[avgfunc, semfunc]=td500_avgSpden(NeuronData,TEOpointer,lsnconfig);
plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(2),avgfunc(1,:),semfunc(1,:),0.25,'k') % 

subplot(2,2,2); hold on % faces
TEpointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'BodyParts')==1 & NeuronData.APindexNumber' < 10);
[avgfunc, semfunc]=td500_avgSpden(NeuronData,TEpointer,lsnconfig);
plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(2),avgfunc(3,:),semfunc(1,:),0.25,'r') % 
TEOpointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'BodyParts')==1 & NeuronData.APindexNumber' > 16);
[avgfunc, semfunc]=td500_avgSpden(NeuronData,TEOpointer,lsnconfig);
plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(2),avgfunc(3,:),semfunc(1,:),0.25,'k') % 

subplot(2,2,3); hold on % faces
TEpointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Objects')==1 & NeuronData.APindexNumber' < 10);
[avgfunc, semfunc]=td500_avgSpden(NeuronData,TEpointer,lsnconfig);
plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(2),avgfunc(5,:),semfunc(1,:),0.25,'r') % 
TEOpointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Objects')==1 & NeuronData.APindexNumber' > 16);
[avgfunc, semfunc]=td500_avgSpden(NeuronData,TEOpointer,lsnconfig);
plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(2),avgfunc(5,:),semfunc(1,:),0.25,'k') % 

subplot(2,2,4); hold on % places
TEpointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Places')==1 & NeuronData.APindexNumber' < 10);
[avgfunc, semfunc]=td500_avgSpden(NeuronData,TEpointer,lsnconfig);
plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(2),avgfunc(4,:),semfunc(1,:),0.25,'r') % 
TEOpointer=find(strcmp(NeuronData.confSensory,'Sensory')==1 & strcmp(NeuronData.confPrefCat,'Places')==1 & NeuronData.APindexNumber' > 16);
[avgfunc, semfunc]=td500_avgSpden(NeuronData,TEOpointer,lsnconfig);
plx_shade_spden(lsnconfig.xScaleRange(1):lsnconfig.xScaleRange(2),avgfunc(4,:),semfunc(1,:),0.25,'k') % 





% latency vs. AP
xdata=NeuronData.APindexNumber(qualityPointer);
ydata=NeuronData.cat_latency(qualityPointer,1:5);

xdataticks=unique(xdata);
graphdata=[];
for ap=min(xdataticks):max(xdataticks),
    graphMean(ap-4,:)=mean(ydata(xdata==ap,:),1);
    graphSTD(ap-4,:)=std(ydata(xdata==ap,:),1);
end
hold on
errorbar(xdataticks,graphMean(:,1),graphSTD(:,1)/sqrt(length(graphSTD)))
errorbar(xdataticks,graphMean(:,2),graphSTD(:,2)/sqrt(length(graphSTD)))
errorbar(xdataticks,graphMean(:,3),graphSTD(:,3)/sqrt(length(graphSTD)))
errorbar(xdataticks,graphMean(:,4),graphSTD(:,4)/sqrt(length(graphSTD)))
errorbar(xdataticks,graphMean(:,5),graphSTD(:,5)/sqrt(length(graphSTD)))

catpointer=[ones(length(xdata),1)*1 ;ones(length(xdata),1)*2 ;ones(length(xdata),1)*3 ...
    ;ones(length(xdata),1)*4 ;ones(length(xdata),1)*5 ];

anovan(reshape(ydata,size(ydata,1)*size(ydata,2),1),{catpointer repmat(xdata,1,5)'},'model',2,'varnames',{'category' 'location'})


export_fig td500_Fig1_Latency.eps -eps -transparent -rgb



return



savefig('')