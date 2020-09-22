function pca_plot_df2(feanames,fscore,fcoefforth,fexplained,D1,D2,pca_group,pca_group_name,pca_group2,pca_group_name2,xloadlim_ratio,yloadlim_ratio,fileName,FontSize,extraLabel,legendList,varargin)
    %parsing the variables
    p = inputParser;
    addRequired(p,'feanames',@iscell);
    addRequired(p,'fscore',@ismatrix);
    addRequired(p,'fcoefforth',@ismatrix);
    addRequired(p,'fexplained',@ismatrix);
    addRequired(p,'D1',@isnumeric);
    addRequired(p,'D2',@isnumeric);
    addRequired(p,'pca_group',@ismatrix);
    addRequired(p,'pca_group_name',@iscell);
    addRequired(p,'pca_group2',@ismatrix);
    addRequired(p,'pca_group_name2',@iscell);
    addRequired(p,'xloadlim_ratio',@isnumeric);
    addRequired(p,'yloadlim_ratio',@isnumeric);
    addRequired(p,'fileName',@ischar);
    addRequired(p,'FontSize',@isnumeric);
    addRequired(p,'extraLabel',@iscell);
    addRequired(p,'legendList',@iscell);
    addParameter(p,'plot_errorElli',1,@isnumeric);
    addParameter(p,'errorElli_conf',0.95,@isnumeric);
    addParameter(p,'LoadingWeight',[],@ismatrix);  
    addParameter(p,'threeD',0,@isnumeric);
    addParameter(p,'legendScale',{'1','2','3','4','5','6','7'},@iscell);
    addParameter(p,'legendText','# challenges',@ischar);
    addParameter(p,'rescaleX',0,@isnumeric);
    addParameter(p,'rescaleY',0,@isnumeric);
    p.KeepUnmatched = true;
    parse(p,feanames,fscore,fcoefforth,fexplained,D1,D2,pca_group,pca_group_name,pca_group2,pca_group_name2,xloadlim_ratio,yloadlim_ratio,fileName,FontSize,extraLabel,legendList,varargin{:});
	
    if p.Results.rescaleX==1
        fscore_ori=fscore;
        fscore(:,1)=log2(abs(fscore(:,1))+1);%if outliert exist in pca plot
        tmp=ones(size(fscore_ori,1),1);
        tmp(fscore_ori(:,1)<0)=-1;
        fscore(:,1)=fscore(:,1).*tmp;
    end
    if p.Results.rescaleY==1
        fscore_ori=fscore;
        fscore(:,2)=log2(abs(fscore(:,2))+1);%if outliert exist in pca plot
        tmp=ones(size(fscore_ori,1),1);
        tmp(fscore_ori(:,2)<0)=-1;
        fscore(:,2)=fscore(:,2).*tmp;
    end
    
    figure();
    subplot(1,2,1);
    [unig,unig_ind]=unique(pca_group,'legacy');
    cmap = parula(length(unig));
    [unig2,unig_ind2]=unique(pca_group2,'legacy');
    markers={'o','^','s','p','h'};
    le={};
    c=1;
    %exL_sort=legendList;
    for e = unig'
        m=1;
        for s=unig2'
            show=intersect(find(pca_group==e),find(ismember(pca_group2,s)));
            if ~isempty(show)>0
                if p.Results.threeD==1
                    hold on;scatter3(fscore(show,1),fscore(show,2),fscore(show,3),200,'Marker',markers{m},'Markerfacecolor',cmap(c,:),'MarkerEdgeColor',[0.8,0.8,0.8]);
                else
                    hold on;scatter(fscore(show,D1),fscore(show,D2),200,'Marker',markers{m},'Markerfacecolor',cmap(c,:),'MarkerEdgeColor',[0.8,0.8,0.8]);
                end
                %le=[le,strcat([pca_group_name{unig_ind(c)},' ',pca_group_name2{unig_ind2(m)}])];
                le=[le,pca_group_name2{unig_ind2(m)}];
                m=m+1;
            end
        end
        c=c+1;
    end
    if ~isempty(extraLabel)
        for i=1:length(fscore(:,1))
            hold on;text(fscore(i,D1)+range(fscore(:,D1))/50,fscore(i,D2)+range(fscore(:,D2))/50,extraLabel{i},'FontSize',FontSize);
        end
    end
    if p.Results.threeD==1
        amin=min(fscore(:,1:3));
        amax=max(fscore(:,1:3));
        xrange=[amin(1)-(amax(1)-amin(1))*0.25,amax(1)+(amax(1)-amin(1))*0.25];
        yrange=[amin(2)-(amax(2)-amin(2))*0.25,amax(2)+(amax(2)-amin(2))*0.25];
        zrange=[amin(3)-(amax(3)-amin(3))*0.25,amax(3)+(amax(3)-amin(3))*0.25];
        xlim(xrange);
        ylim(yrange);
        zlim(zrange);
        xlabel(strcat('Scores on LV',num2str(1),' Variance:',num2str(sprintf('%0.1f',fexplained(1))),'%'));
        ylabel(strcat('Scores on LV',num2str(2),' Variance:',num2str(sprintf('%0.1f',fexplained(2))),'%'));
        zlabel(strcat('Scores on LV',num2str(3),' Variance:',num2str(sprintf('%0.1f',fexplained(3))),'%'));
        hold on;plot3(xrange,[0,0],[0,0],':','color',[0.8,0.8,0.8]);
        hold on;plot3([0,0],yrange,[0,0],':','color',[0.8,0.8,0.8]);
        hold on;plot3([0,0],[0,0],zrange,':','color',[0.8,0.8,0.8]);
        axis tight
    else
        amin=min(fscore(:,[D1,D2]));
        amax=max(fscore(:,[D1,D2]));
        xrange=[amin(1)-(amax(1)-amin(1))*0.1,amax(1)+(amax(1)-amin(1))*0.1];
        yrange=[amin(2)-(amax(2)-amin(2))*0.1,amax(2)+(amax(2)-amin(2))*0.1];
        xlim(xrange);
        ylim(yrange);
        xlabel(strcat('Scores on LV',num2str(D1),' Variance:',num2str(sprintf('%0.1f',fexplained(D1))),'%'));
        ylabel(strcat('Scores on LV',num2str(D2),' Variance:',num2str(sprintf('%0.1f',fexplained(D2))),'%'));
        hline(0,'--black');
        vline(0,'--black');
    end
    %title('PCA');
    if ~isempty(pca_group_name)
        le=unique(le,'stable');
        lh=legend(le,'Location','northeast');
    else
        hcb=colorbar('Ticks',0:1/(length(p.Results.legendScale)-1):1,'TickLabels',p.Results.legendScale);
        hcb.Label.String = p.Results.legendText;
    end
    %Plot error ellipse
    if p.Results.plot_errorElli
        c=1;
        for e = unique(pca_group,'stable')'
            if p.Results.threeD==1
                covariance=cov(fscore(pca_group==e,1:3));
                [elliX,elliY]=error_ellipse(covariance,'conf',p.Results.errorElli_conf,'mu',mean(fscore(pca_group==e,1:3)));
                hold on;plot(elliX,elliY,'-','Color',cmap(c,:));
            else
                covariance=cov(fscore(pca_group==e,[D1,D2]));
                [elliX,elliY]=error_ellipse(covariance,'conf',p.Results.errorElli_conf,'mu',mean(fscore(pca_group==e,[D1,D2])));
                hold on;plot(elliX,elliY,'--','Color',cmap(c,:));
            end
            c=c+1;
        end
    end
    %%%%%%%%%
    subplot(1,2,2);
    cmap = hsv(8);
    clist=ones(length(fcoefforth(:,1)),1);
    for i=1:length(fcoefforth(:,1))
        clist(i)=grouping(feanames{i});
    end
    for i=unique(clist,'legacy')'
        hold on;scatter(fcoefforth(clist==i,D1),fcoefforth(clist==i,D2),150,'Marker','o','Markerfacecolor',cmap(i,:),'MarkerEdgeColor',[0.8,0.8,0.8]);
    end
    for i=1:length(fcoefforth(:,1))
        if clist(i)<3 || clist(i)>5
            hold on;text(fcoefforth(i,D1)+range(fcoefforth(:,D1))/50,fcoefforth(i,D2)+range(fcoefforth(:,D2))/50,strrep(feanames(i),'_','\_'),'FontSize',FontSize);
        end
    end
    if ~isempty(p.Results.LoadingWeight)
        FeaWeight=p.Results.LoadingWeight;
        cmap = colormap(hot(length(FeaWeight.Properties.VariableNames)+5));
        for i=1:length(FeaWeight.Properties.VariableNames)
            ind=find(ismember(feanames,FeaWeight.Properties.VariableNames(i)));
            hold on;scatter(fcoefforth(ind,D1),fcoefforth(ind,D2),'Markerfacecolor',cmap(i,:));
        end 
    end
    xlim(xrange*xloadlim_ratio);
    ylim(yrange*yloadlim_ratio);
    hline(0,'--black');
    vline(0,'--black');
    xlabel(strcat('Loadings on LV',num2str(D1)));
    ylabel(strcat('Loadings on LV',num2str(D2)));
    l={'Ab functions','Fc glycans','Human Fc gamma R','Rhesus Fc gamma R','C1q','ELISA titer','ELISPOT','Isotypes'};
    lh=legend(l(unique(clist,'legacy')),'Location','northeast');
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1],'PaperType','uslegal');
    print(gcf, '-dpdf', strcat(fileName,'.pdf'));
end

function g=grouping(SS)
    if ~isempty(regexp(SS,'(CD107|CD107a|IFN_g|IFNy|MIP1b|MIP1_b|ADCD|ADCC|ADCP|ADNP|Polyfunctionality)','ONCE'))
        g=1;
    elseif ~isempty(regexp(SS,...
    '(^G0|^G1|^G2|Fucose|Bisecting|Di_Sialic_Acid|Mono_Sialic_Acid|Total_Sialic_Acid|Sialic|G0F|G0FB|G1_1|G1_2|G1F_1|G1F_2|G1S1F|G1FB|G2F|G2FB|G2S1|G2S1F|G2S1B|G2S1FB|G2S2|G2S2F|G2S2B|G2S2FB)','ONCE'))
        g=2;
    elseif ~isempty(regexp(SS,...
    '(G0\(FB\)|G1\(-\)|G1\(S\)|G1\(F\)|G1\(B\)|G2\(S\)|G2\(S1\)|G2\(S2\)|G2S\(-\)|G2\(Fx\)|G2\(F-S\)|G2\(F-B\)|G2\(Bx\)|G2\(-B\))','ONCE'))
        g=2;
    elseif ~isempty(regexp(SS,'(^FcgR2AH_|^FcgR2AR_|^FcgR2B_|^FcgR3AV_|^FcgR3B_NA1_|^FcgR3AF_)','ONCE'))
        g=3;
    elseif ~isempty(regexp(SS,'(^R2A_2_|^R2A_3_|^R2A_4_|^R2B_1_|^R3A_1_|^R3A_3_)','ONCE'))
        g=4;
    elseif ~isempty(regexp(SS,'^C1q_','ONCE'))
        g=5;
    elseif ~isempty(regexp(SS,'(^aRhIgG_|^aRhIgA|^ELISA_|^IgG_breadth|^IgA_breadth)','ONCE'))
        g=6;
    elseif ~isempty(regexp(SS,'(^ELISPOT_)','ONCE'))
        g=7;
    elseif ~isempty(regexp(SS,'(IgGtot|IgG1|IgG2|IgG3|IgG4)','ONCE'))
        g=8;
    end
end
