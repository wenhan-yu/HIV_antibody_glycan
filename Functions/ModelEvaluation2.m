%Null testing
function [TruePred_dis,PermutPred1,PermutPred2]=ModelEvaluation2(modelType,cycle,class,feaTable,SelfeaName,SaveName,varargin)
    %parsing the variables
    p = inputParser;
    addRequired(p,'modelType',@ischar);%modelType='da'(classification) or 'r'(regression)
    addRequired(p,'cycle',@isnumeric);%number of random model
    addRequired(p,'class',@ismatrix);
    addRequired(p,'feaTable',@istable);
    addRequired(p,'SelfeaName',@iscell);
    addRequired(p,'SaveName',@ischar);
    addParameter(p,'cvIteration',100,@isnumeric);%numer of cv cycles
    addParameter(p,'plot',1,@isnumeric);
    addParameter(p,'use_r',0,@isnumeric);
    addParameter(p,'cvfold',10,@isnumeric);
    addParameter(p,'plsVL',0,@isnumeric);%the default pls latent variable number
    addParameter(p,'TargetScale',1,@isnumeric);%dependent variable rescale
    p.KeepUnmatched = true;
    parse(p,modelType,cycle,class,feaTable,SelfeaName,SaveName,varargin{:});
    rng(1); %Avoid to repeat a resultsfrom previous matlab session
    
    if p.Results.TargetScale==1
        class=zscore(class);
    else
        class=class;
    end
    features=table2array(feaTable(:,SelfeaName));
    %Detemine # of variable
    if p.Results.plsVL>0
        num=p.Results.plsVL;
    else
        if strcmp(modelType,'da')
            m=PLSCV(features,class,size(features,2),'da');
            [~,num]=max(m.Succv);
        else
            m=PLSCV(features,class,size(features,2));
            [~,num]=min(m.RMSEcv);
        end
    end
    disp(strcat('Best # of latent variables:',num2str(num)));
    pls_cv = PLSCV(features,class,num);%leave-one-out cv
    if strcmp(modelType,'da')
        [~,b]=max(pls_cv.R2cv);
        TrueROC = roccurve(pls_cv.Ycv(:,b),class,length(class),0);
        TruePred=TrueROC.AUC;
        xlab='cross-validation AUC';
    else
        if p.Results.use_r==1
            TruePred=sqrt(max(pls_cv.R2cv));
            xlab='cross-validation correlation coefficient r';
        else
            TruePred=max(pls_cv.R2cv);
            xlab='cross-validation correlation coefficient R ^2';
        end
    end
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    %Comparison between the true/permutated model distribution
    cvfold=p.Results.cvfold;
    if length(class)<cvfold %Special case with small amount of samples
        cvfold=ceil(length(class)/2);%to ensure at least two samples in testing => (#samples/fold) = test size
    end
    TruePred_dis=[];PermutPred1=[];PermutPred2=[];
    %%True model
    for i=1:p.Results.cvIteration
        ind=crossvalind('Kfold', length(class), cvfold);
        TruePred_dis=[TruePred_dis,Kfold_plscv(modelType,features,class,num,ind,p.Results.use_r)];
    end
    
    %Two Null model::
    Exd=find(sum(feaTable{:,:})==0);
    Ind=1:size(feaTable,2);Ind(Exd)=[];
    feaTable=feaTable(:,Ind);%Remove no change features
    for j=1:cycle
        if~rem(j,50);disp(j);end;
        %1>Random mode 1::Randomly selecting same number of features as correlate for modeling
        sel=randperm(size(table2array(feaTable),2),length(SelfeaName));
        feaRond=table2array(feaTable(:,sel));
        %2>Randome model 2: Using the defined correlates but shuffling the index of the outcome
        classRond=class(randperm(length(class)));
        
        for i=1:p.Results.cvIteration/10%Reduce claculation time
            ind=crossvalind('Kfold', length(class), cvfold);
            %tmp=Kfold_plscv(modelType,feaRond,class,num,ind,p.Results.use_r);
            %if tmp>0.9
            %    disp(feaTable(:,sel).Properties.VariableNames);
            %end
            PermutPred1=[PermutPred1,Kfold_plscv(modelType,feaRond,class,num,ind,p.Results.use_r)];
            PermutPred2=[PermutPred2,Kfold_plscv(modelType,features,classRond,num,ind,p.Results.use_r)];
        end
    end
    [~,pvalue1]=ttest2(TruePred_dis,PermutPred1,'Tail','right');
    [~,pvalue2]=ttest2(TruePred_dis,PermutPred2,'Tail','right');
    if p.Results.plot==1
        %The purpose of plotting violin plot
        TruePred_dis=[TruePred_dis,repmat(mean(TruePred_dis)-0.01,1,2)];
        figure;
        grpNew=[ones(1,length(TruePred_dis))+(rand(1,length(TruePred_dis))-0.5)*0.2,...
            ones(1,length(PermutPred1))+1+(rand(1,length(PermutPred1))-0.5)*0.2,...
            ones(1,length(PermutPred2))+2+(rand(1,length(PermutPred2))-0.5)*0.2];
        Y{:,1}=TruePred_dis';Y{:,2}=PermutPred1';Y{:,3}=PermutPred2';
        violin(Y,'xlabel',{'True model','Features random selection','Outcome shuffling'},...
            'facecolor',[0.6,0.6,0.6],'edgecolor','none','mc',[]);
        ylabel(xlab);
        l=max([TruePred_dis,PermutPred1,PermutPred2]);
        hold on;scatter(grpNew,[TruePred_dis,PermutPred1,PermutPred2],7,'blue','filled');
        hold on;line([1,2],[l+0.05,l+0.05],'Color',[0.6 0.6 0.6],'LineStyle','--');
        hold on;line([1,3],[l+0.15,l+0.15],'Color',[0.6 0.6 0.6],'LineStyle','--');
        text(1.4,l+0.07, strcat(['\it p = ',num2str(sprintf('%0.2e',pvalue1))]));
        text(1.9,l+0.17, strcat(['\it p = ',num2str(sprintf('%0.2e',pvalue2))]));
        ylim([min([TruePred_dis,PermutPred1,PermutPred2])-0.1,l+0.2]);
        
        set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1],'PaperType','uslegal');
        print(gcf, '-dpdf', strcat(SaveName,'.pdf'));
    end
end

function outcome=Kfold_plscv(modelType,features,class,Vnum,Kindex,use_r)
    Yobe=class;Ypred=class;
    for eind=unique(Kindex)'
        Ttind=find(Kindex==eind);
        Trind=setdiff(1:length(class),Ttind)';
        predR = plspred(features(Ttind,:),PLS(features(Trind,:),class(Trind),Vnum),class(Ttind));
        Ypred(Ttind)=predR.Yp;
    end
    if strcmp(modelType,'da')
        roc = roccurve(Ypred,Yobe,length(class),0);
        outcome=roc.AUC;
    else
        if use_r==1
            outcome=corr(Ypred,Yobe);
        else
            outcome=corr(Ypred,Yobe)^2;
        end
    end
end