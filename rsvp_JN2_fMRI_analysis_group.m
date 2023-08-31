function rsvp_JN2_fMRI_analysis_group;

hmiconfig=generate_hmi_configplex;
figpath = ['\\.psf\Home\Documents\Manuscripts\RSVP_JN_RevisionNo2\Figure_Source_Images\'];

corrdataS=load([hmiconfig.rsvpanal,filesep,'Stewie_CorrelationData.mat']);
corrdataW=load([hmiconfig.rsvpanal,filesep,'Wiggum_CorrelationData.mat']);
cdS=corrdataS.corrdata; clear corrdataS
cdW=corrdataW.corrdata; clear corrdataW
cdG=[cdS cdW];

figure; clf; cla; set(gcf,'Units','Normalized'); set(gcf,'Position',[0.1 0.1 0.5 0.8]); set(gca,'FontName','Helvetica')
for x=1:3,
    if x==1,corrdata=cdS;
    elseif x==2, corrdata=cdW;
    else corrdata=cdG;
    end
    % Generate Correlation Plots
    % Compress Data
    total_patches=length(corrdata)
    corrdata_mini=struct('panelA_x',[],'panelB_x',[],'panelC_x',[],'panelD_x',[],'panelA_y',[],'panelB_y',[],'panelC_y',[],'panelD_y',[]);
    for tp=1:total_patches,
        corrdata_mini.panelA_x(tp,1:4)=corrdata(tp).fmri_CatSI;
        corrdata_mini.panelA_y(tp,1:4)=corrdata(tp).excite_mean_CatSI;
        corrdata_mini.panelB_x(tp,1:4)=corrdata(tp).fmri_Response;
        corrdata_mini.panelB_y(tp,1:4)=corrdata(tp).excite_mean_NormActivity;
        corrdata_mini.panelC_x(tp,1:4)=corrdata(tp).fmri_Response;
        corrdata_mini.panelC_y(tp,1:4)=corrdata(tp).excite_mean_Activity;
        corrdata_mini.panelD_x(tp,1:4)=corrdata(tp).fmri_Response;
        corrdata_mini.panelD_y(tp,1:4)=corrdata(tp).distribution;
    end

    % reshape
    panelA_x=reshape(corrdata_mini.panelA_x,4*size(corrdata_mini.panelA_x,1),1);
    panelA_y=reshape(corrdata_mini.panelA_y,4*size(corrdata_mini.panelA_y,1),1);
    panelB_x=reshape(corrdata_mini.panelB_x,4*size(corrdata_mini.panelB_x,1),1);
    panelB_y=reshape(corrdata_mini.panelB_y,4*size(corrdata_mini.panelB_y,1),1);
    panelC_x=reshape(corrdata_mini.panelC_x,4*size(corrdata_mini.panelC_x,1),1);
    panelC_y=reshape(corrdata_mini.panelC_y,4*size(corrdata_mini.panelC_y,1),1);
    panelD_x=reshape(corrdata_mini.panelD_x,4*size(corrdata_mini.panelD_x,1),1);
    panelD_y=reshape(corrdata_mini.panelD_y,4*size(corrdata_mini.panelD_y,1),1);

    subplot(3,4,(x-1)*4+1); hold on % compare FMRI SI with Cat SI
    plot(panelA_x,panelA_y,'ks','MarkerFaceColor',[0 0 0])
    polypar=polyfit(panelA_x,panelA_y,1);
    linefit=polyval(polypar,[-0.3 0.3 -0.2 0.2]);
    plot([-0.3 0.3 -0.2 0.2],linefit,'k-','LineWidth',1.5);
    [rho,pval]=corr(panelA_x,panelA_y); % Pearson (linear)
    [rhoK,pvalK]=corr(panelA_x,panelA_y,'type','Kendall'); % Kendall 
    [rhoS,pvalS]=corr(panelA_x,panelA_y,'type','Spearman'); % Spearman (rank)
    text(0.1,0.1,['K=',num2str(rhoK,'%0.2g'),'- (',num2str(pvalK,'%0.2g')])
    text(0.15,0.15,['S=',num2str(rhoS,'%0.2g'),'- (',num2str(pvalS,'%0.2g')])
    xlabel('fMRI Selectivity','FontSize',8); ylabel('Neuronal Selectivity','FontSize',8)
    xlim([-0.3 0.3]); ylim([-0.2 0.2]); axis square
    title({'fMRI vs. Neuronal Selectivity',['rho=',num2str(rho,'%0.2g'),', p=',num2str(pval,'%0.2g')]})

    subplot(3,4,(x-1)*4+2); hold on % compare FMRI activation with Normalized Firing Rate
    plot(panelB_x,panelB_y,'ks','MarkerFaceColor',[0 0 0])
    polypar=polyfit(panelB_x,panelB_y,1);
    linefit=polyval(polypar,[0 3 0.5 1.0]);
    plot([0 3 0.5 1.0],linefit,'k-','LineWidth',1.5);
    [rho,pval]=corr(panelB_x,panelB_y);
    [rhoK,pvalK]=corr(panelB_x,panelB_y,'type','Kendall'); % Kendall 
    [rhoS,pvalS]=corr(panelB_x,panelB_y,'type','Spearman'); % Spearman (rank)
    text(1,0.6,['K=',num2str(rhoK,'%0.2g'),'- (',num2str(pvalK,'%0.2g')])
    text(1,0.65,['S=',num2str(rhoS,'%0.2g'),'- (',num2str(pvalS,'%0.2g')])
    xlabel('fMRI Selectivity','FontSize',8); ylabel('Neuronal Selectivity','FontSize',8)
    xlabel('fMRI Activation (% Signal Change)','FontSize',8); ylabel('Normalized Neuronal Response (sp/s)','FontSize',8)
    xlim([0 3]); ylim([0.5 1.0]); axis square
    title({'fMRI vs. Neuronal Response',['rho=',num2str(rho,'%0.2g'),', p=',num2str(pval,'%0.2g')]})

    subplot(3,4,(x-1)*4+3); hold on % compare FMRI activation with Firing Rate
    plot(panelC_x,panelC_y,'ks','MarkerFaceColor',[0 0 0])
    polypar=polyfit(panelC_x,panelC_y,1);
    linefit=polyval(polypar,[0 3 0 50]);
    plot([0 3 0 50],linefit,'k-','LineWidth',1.5);
    [rho,pval]=corr(panelC_x,panelC_y);
    [rhoK,pvalK]=corr(panelC_x,panelC_y,'type','Kendall'); % Kendall 
    [rhoS,pvalS]=corr(panelC_x,panelC_y,'type','Spearman'); % Spearman (rank)
    text(0.1,0.1,['K=',num2str(rhoK,'%0.2g'),'- (',num2str(pvalK,'%0.2g')])
    text(0.15,0.15,['S=',num2str(rhoS,'%0.2g'),'- (',num2str(pvalS,'%0.2g')])
    xlabel('fMRI Selectivity','FontSize',8); ylabel('Neuronal Selectivity','FontSize',8)
    xlabel('fMRI Activation (% Signal Change)','FontSize',8); ylabel('Normalized Neuronal Response (sp/s)','FontSize',8)
    xlim([0 3]); ylim([0 25]); axis square
    title({'fMRI vs. Neuronal Response',['rho=',num2str(rho,'%0.2g'),', p=',num2str(pval,'%0.2g')]})
    
    subplot(3,4,(x-1)*4+4); hold on % compare FMRI activation with Proportion
    plot(panelD_x,panelD_y,'ks','MarkerFaceColor',[0 0 0])
    polypar=polyfit(panelD_x,panelD_y,1);
    linefit=polyval(polypar,[0 3 0 65]);
    plot([0 3 0 65],linefit,'k-','LineWidth',1.5);
    [rho,pval]=corr(panelD_x,panelD_y);
    [rhoK,pvalK]=corr(panelD_x,panelD_y,'type','Kendall'); % Kendall 
    [rhoS,pvalS]=corr(panelD_x,panelD_y,'type','Spearman'); % Spearman (rank)
    text(0.1,0.1,['K=',num2str(rhoK,'%0.2g'),'- (',num2str(pvalK,'%0.2g')])
    text(0.15,0.15,['S=',num2str(rhoS,'%0.2g'),'- (',num2str(pvalS,'%0.2g')])
    xlabel('fMRI Selectivity','FontSize',8); ylabel('Neuronal Selectivity','FontSize',8)
    xlabel('fMRI Activation (% Signal Change)','FontSize',8); ylabel('Neuronal Distribution','FontSize',8)
    xlim([0 3]); ylim([0 65]); axis square
    title({'fMRI vs. Neuronal Distribution',['rho=',num2str(rho,'%0.2g'),', p=',num2str(pval,'%0.2g')]})


end


print(gcf,[figpath,'GroupCorrelation.ai'],'-dill') % generates an AI file of the figure
print(gcf,[figpath,'GroupCorrelation.jpeg'],'-djpeg') % generates an JPEG file of the figure

return
