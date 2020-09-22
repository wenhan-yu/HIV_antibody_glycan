function plsda_plot(feanames,pls_model,pls_cv,D1,D2,plsda_group,plsda_group_name,xloadlim_ratio,yloadlim_ratio,fileName,FontSize,extraLabel,varargin)
    %parsing the variables
    p = inputParser;
    addRequired(p,'feanames',@iscell);
    addRequired(p,'pls_model',@isstruct);
    addRequired(p,'pls_cv',@isstruct);
    addRequired(p,'D1',@isnumeric);
    addRequired(p,'D2',@isnumeric);
    addRequired(p,'plsda_group',@ismatrix);
    addRequired(p,'plsda_group_name',@iscell);
    addRequired(p,'xloadlim_ratio',@isnumeric);
    addRequired(p,'yloadlim_ratio',@isnumeric);
    addRequired(p,'fileName',@ischar);
    addRequired(p,'FontSize',@isnumeric);
    addRequired(p,'extraLabel',@iscell);
    addParameter(p,'loadingColor',[],@ismatrix);
    addParameter(p,'errorElli_conf',0.95,@isnumeric);
    addParameter(p,'extraText','',@ischar);
    addParameter(p,'legendScale',{'1','2','3','4','5','6','7'},@iscell);
    addParameter(p,'legendText','# challenges',@ischar);
    addParameter(p,'pointsize',200,@isnumeric);
    p.KeepUnmatched = true;
    parse(p,feanames,pls_model,pls_cv,D1,D2,plsda_group,plsda_group_name,xloadlim_ratio,yloadlim_ratio,fileName,FontSize,extraLabel,varargin{:});
    
    figure;
    subplot(1,2,1);
    [unig,unig_ind]=unique(plsda_group,'legacy');
    cmap = colormap(parula(length(unig)*6));
    cmap = cmap((1:length(unig))*6-3,:);
    c=1;
    for e = unig'
        hold on;scatter(pls_model.T(plsda_group==e,D1),pls_model.T(plsda_group==e,D2),p.Results.pointsize,'Marker','o','Markerfacecolor',cmap(c,:),'MarkerEdgeColor',[0.8,0.8,0.8]);
        c=c+1;
    end
    if ~isempty(extraLabel)
        for i=1:size(pls_model.T,1)
            if isnumeric(extraLabel(i))
                hold on;text(pls_model.T(i,D1)+range(pls_model.T(:,D1))/50,pls_model.T(i,D2)+range(pls_model.T(:,D2))/50,num2str(extraLabel(i)),'FontSize',FontSize-2);
            else
                hold on;text(pls_model.T(i,D1)+range(pls_model.T(:,D1))/50,pls_model.T(i,D2)+range(pls_model.T(:,D2))/50,extraLabel(i),'FontSize',FontSize-2);
            end
        end
    end
    amin=min(pls_model.T(:,D1));
    amax=max(pls_model.T(:,D1));
    xrange=[amin-(amax-amin)*0.2,amax+(amax-amin)*0.2];
    amin=min(pls_model.T(:,D2));
    amax=max(pls_model.T(:,D2));
    yrange=[amin-(amax-amin)*0.2,amax+(amax-amin)*0.2];
    xlim(xrange);
    ylim(yrange);
    %title(fileName);
    hline(0,'--black');
    vline(0,'--black');
    
    if D1==1
        XD1Var=pls_model.VLvar(D1,1);
        YD1Var=pls_model.VLvar(D1,2);
    else
        XD1Var=pls_model.VLvar(D1,1)-pls_model.VLvar(D1-1,1);
        YD1Var=pls_model.VLvar(D1,2)-pls_model.VLvar(D1-1,2);
    end
    if D2==1
        XD2Var=pls_model.VLvar(D2,1);
        YD2Var=pls_model.VLvar(D2,2);
    else
        XD2Var=pls_model.VLvar(D2,1)-pls_model.VLvar(D2-1,1);
        YD2Var=pls_model.VLvar(D2,2)-pls_model.VLvar(D2-1,2);
    end
    xlabel({strcat(['Scores on LV',num2str(D1)]),strcat(['X-Variance: ',num2str(sprintf('%0.1f',XD1Var)),'% ',' Y-Variance: ',num2str(sprintf('%0.1f',YD1Var)),'%'])},...
        'FontName','Arial','FontSize',12,'FontWeight','Bold');
    ylabel({strcat(['Scores on LV',num2str(D2)]),strcat(['X-Variance: ',num2str(sprintf('%0.1f',XD2Var)),'% ',' Y-Variance: ',num2str(sprintf('%0.1f',YD2Var)),'%'])},...
        'FontName','Arial','FontSize',12,'FontWeight','Bold');
    if isfield(pls_cv,'RMSEcv')
        hcb=colorbar('Ticks',0:1/(length(p.Results.legendScale)-1):1,'TickLabels',p.Results.legendScale,'FontSize',12);
        hcb.Label.String = p.Results.legendText;
        text(xrange(1)+diff(xrange)*0.02,yrange(1)+diff(yrange)*1.05,strcat(['Calibration RMSE: ',num2str(pls_model.RMSEc)]),'FontWeight','bold');
        text(xrange(1)+diff(xrange)*0.02,yrange(1)+diff(yrange)*1.03,strcat(['Calibration R ^2: ',num2str(pls_model.R2c)]),'FontWeight','bold');
        text(xrange(1)+diff(xrange)*0.02,yrange(1)+diff(yrange)*1.00,strcat(['Leave-1-out CV RMSE: ',num2str(min(pls_cv.RMSEcv))]),'FontWeight','bold');
        text(xrange(1)+diff(xrange)*0.02,yrange(1)+diff(yrange)*0.98,strcat(['Leave-1-out CV R ^2: ',num2str(max(pls_cv.R2cv))]),'FontWeight','bold');
        %text(xrange(1)+diff(xrange)*0.1,yrange(1)+diff(yrange)*0.78,p.Results.extraText,'FontWeight','bold');
    else
        lh=legend(plsda_group_name(unig_ind),'Location','southeast');
        text(xrange(1)+diff(xrange)*0.1,yrange(1)+diff(yrange)*1.05,strcat(['Calibration success: ',num2str(pls_model.Succ),'%']),'FontWeight','bold');
        text(xrange(1)+diff(xrange)*0.1,yrange(1)+diff(yrange)*1.03,strcat(['Leave-1-out CV success: ',num2str(max(pls_cv.Succv(1:size(pls_model.VLvar,1)))),'%']),'FontWeight','bold');
        %Show only in PLSDA model
        %Plot error ellipse
        c=1;
        for e = unig'
            covariance=cov(pls_model.T(plsda_group==e,[D1,D2]));
            [elliX,elliY]=error_ellipse(covariance,'conf',p.Results.errorElli_conf,'mu',mean(pls_model.T(plsda_group==e,D1:D2)));
            hold on;plot(elliX,elliY,'--','Color',cmap(c,:));
            c=c+1;
        end
    end
    
    %%%%%%%%%
    subplot(1,2,2);
    if isempty(p.Results.loadingColor)
        scatter(pls_model.P(:,D1),pls_model.P(:,D2),p.Results.pointsize,'Marker','o','Markerfacecolor','black','MarkerEdgeColor',[0.8,0.8,0.8]);
    else
        for i=1:length(pls_model.P(:,1))
            hold on;scatter(pls_model.P(i,D1),pls_model.P(i,D2),p.Results.pointsize,'Marker','o','Markerfacecolor',p.Results.loadingColor(i,:),'MarkerEdgeColor',[0.8,0.8,0.8]);
        end
    end
    for i=1:length(pls_model.P(:,1))
        hold on;text(pls_model.P(i,D1)+range(pls_model.P(:,D1))/50,pls_model.P(i,D2)+range(pls_model.P(:,D2))/50,strrep(feanames(i),'_','\_'),'FontSize',FontSize);
    end
    xlim(xrange*xloadlim_ratio);
    
    ylim(yrange*yloadlim_ratio);
    hline(0,'--black');
    vline(0,'--black');
    xlabel(strcat('Loadings on LV',num2str(D1)),'FontName','Arial','FontSize',12,'FontWeight','Bold');
    ylabel(strcat('Loadings on LV',num2str(D2)),'FontName','Arial','FontSize',12,'FontWeight','Bold');
    %set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    set(gcf,'Units', 'Inches','PaperUnits','inches','PaperPosition', [0 0 12 8],'PaperSize', [12,8]);
    print(gcf, '-dpdf', strcat(fileName,'.pdf'));
end