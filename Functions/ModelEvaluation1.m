%Generate PLSDA/R plot, VIP score, prediction plot and correlation plot

function [VIPs,colorMap,pls_model]=ModelEvaluation1(ModelType, bestfeaName, FeaturesT, FeaturesZT, class, className, ExtraLabel,FileName,varargin)
    %parsing the variables
    p = inputParser;
    addRequired(p,'ModelType',@ischar);
    addRequired(p,'bestfeaName',@iscell);
    addRequired(p,'FeaturesT',@ismatrix);
    addRequired(p,'FeaturesZT',@ismatrix);
    addRequired(p,'class',@ismatrix);
    addRequired(p,'className',@iscell);
    addRequired(p,'ExtraLabel',@iscell);
    addRequired(p,'FileName',@ischar);
    addParameter(p,'outcomeAdj',1,@isnumeric);%Sometime outcome display n the plot may be to adjust
    addParameter(p,'group',[],@ismatrix);
    addParameter(p,'groupName',{},@iscell);
    addParameter(p,'PLSDAD1',1,@isnumeric);
    addParameter(p,'PLSDAD2',2,@isnumeric);
    addParameter(p,'PLSDA_xr',15,@isnumeric);
    addParameter(p,'PLSDA_yr',10,@isnumeric);
    addParameter(p,'errorElli_conf',0.95,@isnumeric);
    addParameter(p,'legendScale',{},@iscell);
    addParameter(p,'legendText','',@ischar);
    addParameter(p,'cvfold',10,@isnumeric);
    addParameter(p,'cviteration',100,@isnumeric);
    addParameter(p,'TargetScale',1,@isnumeric);%dependent variable rescale
    addParameter(p,'feaName_adjust',{},@iscell);%Specify match and replace of the string ex.{{'-','_'},{'\-','\_'}}
    addParameter(p,'plsVL',0,@isnumeric);%the default pls latent variable number
     
    %Plot parameter
    addParameter(p,'Cxlab','',@ischar);
    addParameter(p,'Cylab','',@ischar);
    addParameter(p,'km',1,@isnumeric);%Kmplot
    addParameter(p,'KMxlab','',@ischar);
    addParameter(p,'KMylab','',@ischar);
    addParameter(p,'kmAdj',1,@isnumeric);
    addParameter(p,'KMxtick',1,@isnumeric);
    
    addParameter(p,'altVarName',{},@iscell);
    p.KeepUnmatched = true;
    parse(p,ModelType,bestfeaName, FeaturesT, FeaturesZT, class, className, ExtraLabel,FileName,varargin{:});
    rng(1); %Avoid to repeat a resultsfrom previous matlab session
    
    Oric=class;
    FeaturesZ=table2array(FeaturesZT(:,bestfeaName));
    %Convert co-correlate labels if necessary
    bestfeaNameM=regexprep(bestfeaName,'_D_',' \/ ');
    bestfeaNameM=regexprep(bestfeaNameM,'_S_',' \* ');
    
    if ~isempty(p.Results.feaName_adjust)
        for i=1:length(p.Results.feaName_adjust{1})
            bestfeaNameM=strrep(bestfeaNameM,p.Results.feaName_adjust{1}{i},p.Results.feaName_adjust{2}{i});
        end
    end
    if ~isempty(p.Results.altVarName)
        bestfeaNameM=p.Results.altVarName;
    end
    %% Generate PLSDA/R plot
    if strcmp(ModelType,'da')
        pls_cv=PLSCV(FeaturesZ,class,size(FeaturesZ,2),'da');
        if p.Results.plsVL>0
            num=p.Results.plsVL;
        else
            [~,n]=max(pls_cv.Succv);
            if n<3; num=3; else num=n; end%consider the LV # from 1 to at least 3
            if size(FeaturesZ,2)<3;num=size(FeaturesZ,2);end %if fea # is less than LV
        end
        pls_model = PLS(FeaturesZ,class,num,'da');
        %VIP scores
        [VIPs,~,~,bfInd,colorMap]=VIPscore(pls_model,bestfeaNameM,FeaturesZ,class,'correlation',strcat(FileName,'_VIP'));
        plsda_plot(bestfeaNameM,pls_model,pls_cv,p.Results.PLSDAD1,p.Results.PLSDAD2,class,...
        className,p.Results.PLSDA_xr,p.Results.PLSDA_yr,strcat(FileName,'_plsda'),10,ExtraLabel,'errorElli_conf',p.Results.errorElli_conf,'loadingColor',colorMap);
    else
        if p.Results.TargetScale==1
            adjclass=zscore(class);
        else
            adjclass=class;
        end
        
        pls_cv=PLSCV(FeaturesZ,adjclass,size(FeaturesZ,2));
        if p.Results.plsVL>0
            num=p.Results.plsVL;
        else
            [~,n]=min(pls_cv.RMSEcv);
            if n<3; num=3; else num=n; end%consider the LV # from 1 to at least 3
            if size(FeaturesZ,2)<3;num=size(FeaturesZ,2);end %if fea # is less than LV
        end
        pls_cv=PLSCV(FeaturesZ,adjclass,num);
        pls_model = PLS(FeaturesZ,adjclass,num);
        %VIP scores
        [VIPs,~,~,bfInd,colorMap]=VIPscore(pls_model,bestfeaNameM,FeaturesZ,class,'correlation',strcat(FileName,'_VIP'));
        disp(strcat(FileName,'_plsr'))
        plsda_plot(bestfeaNameM,pls_model,pls_cv,p.Results.PLSDAD1,p.Results.PLSDAD2,class,...
        className,p.Results.PLSDA_xr,p.Results.PLSDA_yr,strcat(FileName,'_plsr'),10,ExtraLabel,'extraText','',...
        'legendScale',p.Results.legendScale,'legendText',p.Results.legendText,'loadingColor',colorMap);
    end
    bestfeaName=bestfeaName(bfInd);
    %% ROC or correlation/survive curve plot
    %Detemine # of variable
    if strcmp(ModelType,'da')
        m=PLSCV(FeaturesZ,class,size(FeaturesZ,2),'da');
        [~,num]=max(m.Succv);
    else
        if p.Results.TargetScale==1%Rescale outcome 
            class=zscore(class);
        end
        pls_leave1out=PLSCV(FeaturesZ,class,size(FeaturesZ,2));
        [~,num]=min(pls_leave1out.RMSEcv);
    end
    cvfold=p.Results.cvfold;
    if length(class)<cvfold %Special case with small amount of samples
        cvfold=ceil(length(class)/2);%to ensure at least two samples in testing => (#samples/fold) = test size
    end  
    %%CV for plotting confidence
    Ypred=[];ROCx=[];ROCy=[];AUC=[];accuracy=[];
    for i=1:p.Results.cviteration
        ind=crossvalind('Kfold', length(class), cvfold);
        tmp=class;
        for eind=unique(ind)'
            Ttind=find(ind==eind);
            Trind=setdiff(1:length(class),Ttind)';
            predR = plspred(FeaturesZ(Ttind,:),PLS(FeaturesZ(Trind,:),class(Trind),num),class(Ttind));
            tmp(Ttind)=predR.Yp;
        end
        if strcmp(ModelType,'da')
            roc = roccurve(class,tmp,length(class),0);
            ROCx=[ROCx,1-roc.value(:,1)];ROCy=[ROCy,roc.value(:,2)];AUC=[AUC,roc.AUC];accuracy=[accuracy,roc.accuracy];
        end
        if ~strcmp(ModelType,'da') && p.Results.TargetScale==1
            tmp=tmp.*std(Oric)+mean(Oric);
        end
        Ypred=[Ypred,tmp];
    end
    class=Oric;
    figure;
    if strcmp(ModelType,'da')
        subplot(1,2,1);%ROC
        ROCcurve(ROCx,ROCy,AUC,accuracy,cvfold);
        subplot(1,2,2);%boxplot
        class_boxplot(class,className,Ypred,1,'Classification prediction',1);
    else
        subplot(1,2,1);%Scatter plot
        %plotCorrelation('Correlation',p.Results.Cxlab,p.Results.Cylab,Ypred,class*p.Results.outcomeAdj,0,p.Results.group,p.Results.groupName,cvfold);
        plotCorrelation('Correlation',p.Results.Cxlab,p.Results.Cylab,pls_leave1out.Ycv(:,num).*std(Oric)+mean(Oric),...
            class*p.Results.outcomeAdj,1,p.Results.group,p.Results.groupName,cvfold);
        
        if p.Results.km
            subplot(1,2,2);%Kaplan-Meier survival plot
            KM_plot(class,Ypred,p.Results.KMxlab,p.Results.KMylab,p.Results.kmAdj,p.Results.KMxtick,p.Results.group,p.Results.groupName);
        end
    end
    %set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1]);
    set(gcf,'Units', 'Inches','PaperUnits','inches','PaperPosition', [0 0 12 8],'PaperSize', [12,8]);
    print(gcf, '-dpdf', strcat(FileName,'_PredRobust.pdf'));
    %% Correlation plot against outcome
    %LV1 and top three features against outcome
    if strcmp(ModelType,'da')
        figure;
        subplot(2,2,1);
        class_boxplot(class,className,pls_model.T(:,1),0,'LV1',0);
        for i=1:3
            if i<=length(bestfeaName)
                subplot(2,2,1+i);
                bfN=strrep(regexprep(regexprep(bestfeaName{i},'_D_',' \/ '),'_S_',' \* '),'_','\_');
                class_boxplot(class,className,table2array(FeaturesT(:,bestfeaName{i})),0,bfN,0);
            end
        end
        set(gcf,'Units', 'Inches','PaperUnits','inches','PaperPosition', [0 0 12 10],'PaperSize', [12,10]);
        print(gcf, '-dpdf', strcat(FileName,'_PredCorrelation.pdf'));
    else
        figure;
        subplot(2,2,1);
        plotCorrelation(strcat([p.Results.Cxlab,' vs LV1']),p.Results.Cxlab,'LV1',pls_model.T(:,1),class*p.Results.outcomeAdj,1,p.Results.group,p.Results.groupName,'');
        for i=1:3
            if i<=length(bestfeaName)
                subplot(2,2,1+i);
                bfN=strrep(regexprep(regexprep(bestfeaName{i},'_D_',' \/ '),'_S_',' \* '),'_','\_');
                plotCorrelation(strcat([p.Results.Cxlab,' vs ',bfN]),p.Results.Cxlab,bfN,...
                    table2array(FeaturesT(:,bestfeaName{i})),class*p.Results.outcomeAdj,1,p.Results.group,p.Results.groupName,'');
            end
        end
        set(gcf,'Units', 'Inches','PaperUnits','inches','PaperPosition', [0 0 12 10],'PaperSize', [12,10]);
        print(gcf, '-dpdf', strcat(FileName,'_PredCorrelation.pdf'));
        if length(bestfeaName)>3
            cc=1;ext=1;
            for i=4:length(bestfeaName)
                if cc==1
                    figure;
                end
                subplot(2,2,cc);
                bfN=strrep(regexprep(regexprep(bestfeaName{i},'_D_',' \/ '),'_S_',' \* '),'_','\_');
                plotCorrelation(strcat([p.Results.Cxlab,' vs ',bfN]),p.Results.Cxlab,bfN,...
                    table2array(FeaturesT(:,bestfeaName{i})),class*p.Results.outcomeAdj,1,p.Results.group,p.Results.groupName,'');
                if cc==4 || i==length(bestfeaName)
                    set(gcf,'Units', 'Inches','PaperUnits','inches','PaperPosition', [0 0 12 10],'PaperSize', [12,10]);
                    print(gcf, '-dpdf', strcat(FileName,'_PredCorrelation_',num2str(ext),'.pdf'));
                    cc=1;ext=ext+1;
                else
                    cc=cc+1;
                end
            end
        end
    end
    %% Show bf heatmap corresponding to the clinical outcome
    cmap = colormap(parula(256));
    [~,b]=sort(class,'descend');
    SortFeaZ=FeaturesZT{b,bestfeaName};
    SortFeaZ=(SortFeaZ-repmat(min(SortFeaZ),size(SortFeaZ,1),1))./repmat(max(SortFeaZ)-min(SortFeaZ),size(SortFeaZ,1),1);
    SortFeaZ=(SortFeaZ-0.5)/0.5;

    figure;heatmap(SortFeaZ,bestfeaName,strrep(strcat(FeaturesZT.Properties.RowNames(b),':',num2str(class(b))),'_','\_'), {}, 'TickAngle', 90,...
            'ShowAllTicks',true,'TickFontSize',10,'Colormap', cmap, 'GridLines','none','Colorbar',true);
    set(gca, 'Ticklength', [0 0]);
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1],'PaperType','uslegal','PaperPositionMode', 'auto');
    print(gcf, '-dpdf', strcat(FileName,'_bf_Heatmap.pdf'));

end
function ROCcurve(ROCx,ROCy,AUC,accuracy,cvfold)
    %Calculate 95% confidence interval
    SEM = std(ROCy')/sqrt(size(ROCy,2));% Standard Error
    ts = tinv([0.025  0.975],size(ROCy,2)-1);% T-Score for 95% interval
    y_b = repmat(ts,size(ROCy,1),1).*repmat(SEM',1,2);       
    x_m=mean(ROCx,2);y_m=mean(ROCy,2);
    
    plot(x_m, y_m,'linestyle','-','color','blue','LineWidth',2)
    hold on;boundedline(x_m, y_m, y_b(:,2), 'alpha','transparency', 0.3);
    hold on;line([0,1],[0,1],'linestyle','--','color','black','LineWidth',1);
    ylim([0,1]);
    title(strcat(['ROC ',num2str(cvfold),'-fold CV']),'FontName','Arial','FontSize',14,'FontWeight','Bold');
    xlabel('1-specificity','FontName','Arial','FontSize',10,'FontWeight','Bold');
    ylabel('sensitivity','FontName','Arial','FontSize',10,'FontWeight','Bold');
    legend({strcat('Accuracy:',num2str(sprintf('%0.1f',mean(accuracy)*100)),'% AUC:',num2str(sprintf('%0.2f',mean(AUC))))},...
        'Location','SouthEast');
    axis square
end
function class_boxplot(class,className,Ypred,conf,ylab,vl)
    if size(Ypred,2)>1%iteration results 
        SEM = std(Ypred')/sqrt(size(Ypred,2));% Standard Error
        ts = tinv([0.025  0.975],size(Ypred,2)-1);% T-Score for 95% interval
        y_b = repmat(ts,size(Ypred,1),1).*repmat(SEM',1,2);       
        Y=mean(Ypred,2);
    else
        Y=Ypred;
    end
    
    [unig,unig_ind]=unique(class,'legacy');
    cmap = colormap(parula(length(unig)*6));
    cmap = cmap((1:length(unig))*6-3,:);
    cN=className(unig_ind);
	
    grp1=[];grp2=[];grp2_r=[];
    for i=1:length(unig)
        cc=unig(i);
        grp1=[grp1,Y(find(class==cc))'];
        grp2=[grp2,ones(1,length(find(class==cc)))+i];
        r=ones(1,length(find(class==cc)))+i-1+(rand(1,length(find(class==cc)))-0.5)*0.2;
        grp2_r=[grp2_r,r];
        hold on;scatter(r,Y(find(class==cc))',40,'Marker','o','Markerfacecolor',cmap(i,:),'MarkerEdgeColor',[0.8,0.8,0.8]);
    end
    boxplot(grp1,grp2,'labels',cN');
    if vl==1%Draw vline 
        hold on;plot(get(gca,'xlim'),[mean(unig) mean(unig)],'linestyle','--','color','black','LineWidth',1);
    end
    ylabel(ylab);
    
    %ylim([-0.05,1.05]);
    if length(unig)==2 && conf==1
        lable=[];%[True, false]
        ycoord=[];
        for i=1:length(unig)
            cc=unig(i);
            v1=length(find(Y(find(class==cc))< mean(unig)));
            v2=length(find(Y(find(class==cc))>= mean(unig)));
            lable=[lable;[v1,v2]];
            ycoord=[ycoord;mean(Y(find(class==cc)))];
        end
        text(1.15,mean(unig)-0.1,{strcat(['True labels = ',num2str(lable(1,1))])},'FontWeight','bold');
        text(1.15,mean(unig)+0.1,{strcat(['False labels = ',num2str(lable(1,2))])},'FontWeight','bold');
        text(2.15,mean(unig)-0.1,{strcat(['False labels = ',num2str(lable(2,1))])},'FontWeight','bold');
        text(2.15,mean(unig)+0.1,{strcat(['True labels = ',num2str(lable(2,2))])},'FontWeight','bold');
    end
    axis square
end
function plotCorrelation(tl,xn,yn,Y,X,fitline,group,groupName,cvfold)
    sg={};
    if ~isempty(group)
        [unig,unig_ind]=unique(group,'legacy');
        gN=groupName(unig_ind);
        for i=unig'
            sg=[sg,find(group==i)];
        end
        cmap = colormap(parula(length(unig)*6));
        cmap = cmap((1:length(unig))*6-3,:);
    else
        sg=[sg,1:length(X)];
        cmap=[0 0 1];
    end
    shape={'o','s','^','h'};
    cmap = colormap(colorcube(6));
    if size(Y,2)>1%iteration results
        %Calculate 95% confidence interval
        SEM = std(Y')/sqrt(size(Y,2));% Standard Error
        ts = tinv([0.025  0.975],size(Y,2)-1);% T-Score for 95% interval
        y_b = repmat(ts,size(Y,1),1).*repmat(SEM',1,2);       
        Y=mean(Y,2);
        for i=1:length(sg)
            X_i=X(sg{i});
            Y_i=Y(sg{i},:);
            y_b_i=y_b(sg{i},2);
            hold on;errorbar(X_i,Y_i,y_b_i,'Marker','o','Markerfacecolor',cmap(i,:),'MarkerEdgeColor',[0.5,0.5,0.5],'LineStyle','none');
        end
    else
        for i=1:length(sg)
            X_i=X(sg{i});
            Y_i=Y(sg{i});
            %hold on;scatter(X_i,Y_i,40,'Marker','o','Markerfacecolor',cmap(i,:),'MarkerEdgeColor',[0.5,0.5,0.5]);
            hold on;scatter(X_i,Y_i,100,cmap(i,:),'filled',shape{i});
            %hold on;scatter(X_i,Y_i,100,'Marker',shape{i},'Markerfacecolor',[0 0 0],'MarkerEdgeColor',[0.5,0.5,0.5]);
        end
    end
    [pr,pp]=corr(X,Y,'type','Spearman');
    xrange=[min(X)-(max(X)-min(X))*0.1,max(X)+(max(X)-min(X))*0.1];
    yrange=[min(Y)-(max(Y)-min(Y))*0.1,max(Y)+(max(Y)-min(Y))*0.1];
    xlim(xrange);
    ylim(yrange);
    axis square
    hold on;text(xrange(1)+diff(xrange)*0.2,yrange(1)+diff(yrange)*0.8,...
        {strcat(['\it r = ',num2str(sprintf('%0.2f',pr))]);strcat(['\it p = ',num2str(sprintf('%0.2e',pp))])},'FontSize',10);
    if isempty(cvfold)
        title(tl,'FontName','Arial','FontSize',14,'FontWeight','Bold');
    else
        title(strcat([tl,' ',num2str(cvfold),'-fold cv']),'FontName','Arial','FontSize',14,'FontWeight','Bold');
    end
    xlabel(xn,'FontName','Arial','FontSize',12,'FontWeight','Bold');
    ylabel(yn,'FontName','Arial','FontSize',12,'FontWeight','Bold');
    
    if fitline==1
        %Fitting
        f=fit(X,Y,'poly1');
        nx=linspace(min(X),max(X),100);
        bound= predint(f,nx,0.95,'functional','off');
        hold on;line(nx,f(nx),'LineWidth',1,'Color',[0.5,0.5,0.5]);
        hold on, line(nx,bound(:,1),'LineWidth',1,'LineStyle','--','Color','blue');
        hold on, line(nx,bound(:,2),'LineWidth',1,'LineStyle','--','Color','blue');
    end
    if ~isempty(group)
        legend(gN,'Location','southeast');
    end
end
function KM_plot(class,Ypred,xlab,ylab,kmAdj,sAdj,group,groupName)
    %kmAdj:outcome adjustment for KM
    %sAdj:Xlabel tick adjutment  
    sg={};
    if ~isempty(group)
        [unig,unig_ind]=unique(group,'legacy');
        gN=groupName(unig_ind);
        for i=unig'
            sg=[sg,find(group==i)];
        end
        cmap = colormap(parula(length(unig)*6));
        cmap = cmap((1:length(unig))*6-3,:);
    else
        sg=[sg,1:length(class)];
        cmap=[0 0 1];
    end
    t=[];
    for i=1:length(sg)
        class_i=class(sg{i});
        Ypred_i=Ypred(sg{i},:);
        [t1,T1_act]=KMpro(class_i,kmAdj,class_i);
        T1_preds=[];
        for j=1:size(Ypred_i,2)
            [~,T1_i]=KMpro(Ypred_i(:,j),kmAdj,class_i);
            T1_preds=[T1_preds;T1_i];
        end
        %Calculate 95% confidence interval
        SEM = std(T1_preds)/sqrt(size(T1_preds,1));% Standard Error
        ts = tinv([0.025  0.975],size(T1_preds,1)-1);% T-Score for 95% interval
        T1_preds_b = repmat(ts,size(T1_preds,2),1).*repmat(SEM',1,2);
        T1_preds_m=mean(T1_preds);
        %t1 adjust for x ticks
        t1=t1*(sAdj/kmAdj);
        hold on;
        h((i-1)*4+1)=stairs(t1,T1_act,'LineWidth',2,'LineStyle','-','Color',cmap(i,:)); %Kaplan-Meier survival function actual
        hold on;
        h((i-1)*4+2)=stairs(t1,T1_preds_m,'LineWidth',2,'LineStyle','--','Color',cmap(i,:)); %mean prediction
        hold on;
        h((i-1)*4+3)=stairs(t1,T1_preds_m+T1_preds_b(:,1)','LineWidth',1,'LineStyle',':','Color',[0.5 0.5 0.5]); %95% confidence
        hold on;
        h((i-1)*4+4)=stairs(t1,T1_preds_m+T1_preds_b(:,2)','LineWidth',1,'LineStyle',':','Color',[0.5 0.5 0.5]); %95% confidence
        t=[t,t1];
    end
    %set the axis properly
    xmax=max(t)+1;
    xlim([1,7]);
    axis([0 xmax 0 1.2]);
    axis square
    %add labels and legend
    title('Kaplan-Meier function','FontName','Arial','FontSize',14,'FontWeight','Bold'); 
    ylabel(ylab,'FontName','Arial','FontSize',12,'FontWeight','Bold'); 
    xlabel(xlab,'FontName','Arial','FontSize',12,'FontWeight','Bold'); 
    if ~isempty(group)
        ind=[];le={};c=1;
        for e=gN'
            e=e{:};
            ind=[ind,[c,c+1]];
            c=c+4;
            le=horzcat(le,{strcat([e,' observed']),strcat([e,' predicted'])});
        end
        ind=[ind,c-1];le=horzcat(le,'95% confidence');
        legend(h(ind),le);
    end
end
function [t1,T1]=KMpro(X,kmAdj,ref)
    %Ref: keep X interval consistent
    ref=round(ref*kmAdj);
    X=round(X*kmAdj);
    if min(ref)>0
        t1=[0,min(ref):max(ref)];%this is the x variable (time);
    else
        t1=[min(ref)-1,min(ref):max(ref)]; 
    end
    T1=ones(1,length(t1));
    for i=2:length(t1)-1
       T1(i)=sum(X>t1(i))/length(ref);
    end
    T1(length(t1))=T1(length(t1)-1);%The last point keep same number as n-1 point
end
