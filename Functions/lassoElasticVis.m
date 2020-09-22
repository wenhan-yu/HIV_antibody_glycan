function [bestnames,pls_model_s,feaWeight,num]=lassoElasticVis(modelType,class,class_name,alt_class,alt_class_name,Rc_features_cz,feanames_c,SaveName,varargin)
    %parsing the variables
    p = inputParser;
    addRequired(p,'modelType',@ischar);%modelType='da'(classification) or 'r'(regression)
    addRequired(p,'class',@ismatrix);
    addRequired(p,'class_name',@iscell);
    addRequired(p,'alt_class',@ismatrix);
    addRequired(p,'alt_class_name',@iscell);
    addRequired(p,'Rc_features_cz',@ismatrix);
    addRequired(p,'feanames_c',@iscell);
    addRequired(p,'SaveName',@ischar);
    addParamValue(p,'plsVL',0,@isnumeric);%the default pls latent variable number
    addParamValue(p,'compute_cycle',0,@isnumeric);%0 means that # cycles to cover all combinations will be used
    addParameter(p,'Kfold',10,@isnumeric);
    addParamValue(p,'rObj','MSE',@ischar);%MSE or R2 for plsr model
    addParamValue(p,'PCAD1',1,@isnumeric);
    addParamValue(p,'PCAD2',2,@isnumeric);
    addParamValue(p,'PLSDAD1',1,@isnumeric);
    addParamValue(p,'PLSDAD2',2,@isnumeric);
    addParamValue(p,'PCA_xr',1,@isnumeric);
    addParamValue(p,'PCA_yr',1,@isnumeric);
    addParamValue(p,'PLSDA_xr',10,@isnumeric);
    addParamValue(p,'PLSDA_yr',10,@isnumeric);
    addParamValue(p,'ElasticLASSO','ElasticNet',@ischar);
    addParamValue(p,'LASSORelTol',1e-4,@isnumeric);%Convergence threshold for LASSO/ElasticNet
    addParamValue(p,'plotall',0,@isnumeric);
    addParamValue(p,'optimalModel',1,@isnumeric);
    addParamValue(p,'SetFea','min',@ischar);
    addParamValue(p,'legendScale',{},@iscell);
    addParamValue(p,'legendText','',@ischar);
    addParamValue(p,'feaName_adjust',{},@iscell);%Specify match and replace of the string ex.{{'-','_'},{'\-','\_'}}
    p.KeepUnmatched = true;
    parse(p,modelType,class,class_name,alt_class,alt_class_name,Rc_features_cz,feanames_c,SaveName,varargin{:});
    rng(1);
    modelSel=upper(p.Results.ElasticLASSO);
    %% %%%%%%%%%%%%%%%%%%%%%%%%%
    %Determine number of cycles to cover all possible combinations for
    %randomly splitting the data into testing and training
    TestSize=ceil(length(class)/p.Results.Kfold);
    m=length(class);
    n=TestSize;
    g=length(unique(class));
    if strcmp(modelType,'r')%contineous number
        %C(m,n)= m!/(n!*(m-n)!), m=total sample size, n=testing size
        ReqSize=factorial(m)/(factorial(n)*factorial(m-n));
    else%category
        %m=total sample size, n=testing size, Gm=the sample size in group m
        %C(m,n)-(C(Gm,n), m=1...i, if Gm>n)
        %All possible combinations minus the combinations from only one group
        ReqSize=factorial(m)/(factorial(n)*factorial(m-n));
        for e=unique(class)'
            if length(find(class==e))>= n && n>1
                ReqSize=ReqSize-factorial(length(find(class==e)))/(factorial(n)*factorial(length(find(class==e))-n));
            end
        end
    end
    disp(strcat(['To cover all combinations needs ',num2str(ReqSize),' iterations.'])); 
	%% %%%%%%%%%%%%%%%%%%%%%%%%%
    %Determine the combination of training and testing
    TrainList=[];TestList=[];%%Split training testing data
    if p.Results.compute_cycle==0 || ReqSize<p.Results.compute_cycle %Determine # iterations if cycle are not given or ReqSize small than the given cycle
        cycle=ReqSize;
    else
        cycle=p.Results.compute_cycle;
    end
    if strcmp(modelType,'da')%
        while size(TestList,2)<cycle
            ind=sort(randperm(m,n))';
            if length(unique(class(ind)))==g || n==1 %Need to cover all classes
                TestList=[TestList,ind];
            end
            TestList=unique(TestList','rows')';
        end
    else
        while size(TestList,2)<cycle
            TestList=[TestList,sort(randperm(m,n))'];
            TestList=unique(TestList','rows')';
        end
    end
    for i=1:cycle %Capture Training list
        TrainList=[TrainList,setdiff(1:m,TestList(:,i))'];
    end
    disp(strcat(['The model will take ',num2str(cycle),' iterations.']));   
    %% %%%%%%%%%%%%%%%%%%%%%%%%%
    %Determine alpha
    if strcmp(modelSel,'ELASTICNET')
        alpha=EvaAlpha(Rc_features_cz,class,p.Results.Kfold);%Determine alpha
    elseif strcmp(modelSel,'LASSO')
        alpha=1;
    else
        alpha=0.01;
    end
    %% %%%%%%%%%%%%%%%%%%%%%%%%%
    %Determine latent vairable #
    if size(Rc_features_cz,2)>30; tt=30; else tt=size(Rc_features_cz,2); end;%For the purpose of reducing calculation time
    if p.Results.plsVL==0
        if strcmp(modelType,'da')
            m=PLSCV(Rc_features_cz,class,tt,'da');
            [~,num]=max(m.Succv);
        else
            m=PLSCV(Rc_features_cz,class,tt);
            [~,num]=min(m.RMSEcv);
        end
        if num<3;num=3;end %To have at least three dimensions
    else
        num=p.Results.plsVL;
    end
    disp(strcat(['PLS LVs in the model=',num2str(num),' Elistic-Net alpha:',num2str(alpha)]));
    %% %%%%%%%%%%%%%%%%%%%%%%%%%
    %Start the iterations
    featureWeight=zeros(cycle,length(feanames_c));
    freq=zeros(cycle,length(feanames_c));
    for c=1:cycle
        if~rem(c,1);disp(strcat('cycle:',num2str(c)));end;
        Train=TrainList(:,c);Test=TestList(:,c);
        %Inner CV
        %Elastic-net/LASSO Feature reduction
        if length(Train)<p.Results.Kfold %For some small sample size
            lasscv_fold=length(Train);
        else
            lasscv_fold=p.Results.Kfold;
        end
        [b,fitinfo]=lasso(Rc_features_cz(Train,:),class(Train),'CV',lasscv_fold,'Alpha',alpha,'RelTol',p.Results.LASSORelTol);
        %if p.Results.plotall
        %    lassoPlot(b,fitinfo,'PlotType','Lambda','XScale','log');
        %    lassoPlot(b,fitinfo,'PlotType','CV');
        %end
        if p.Results.optimalModel==0%Optimal model based on LASSO-MSE
            bestmodel=fitinfo.IndexMinMSE;
        else
            set='max';%if minimal cross-validation error has more than one models, choose the model has the max/min features
            bestmodel=estimateBestmodel(p.Results.rObj,b,fitinfo,Rc_features_cz(Train,:),class(Train),set,modelType,num);
        end
        %identify columns containing the best features
        col_bestmodel=b(:,bestmodel);
        r = find(abs(col_bestmodel)>0);
        % pull out names of best features 
        In_bestnames=feanames_c(r);
        In_bestfeatures=Rc_features_cz(:,r');
        
        %Outer CV
        if num>length(r)%Just in case for the best of features # is small than required LVs #
            numM=length(r);
        else
            numM=num;
        end
        if strcmp(modelType,'da')
            pls_model = PLS(In_bestfeatures(Train,:),class(Train),numM,'da');
            predict = plspred(In_bestfeatures(Test,:),pls_model,class(Test));
            %Variable Importance in Projection (VIP) score for each features 
            %Yvar=[pls_model.VLvar(1,2);diff(pls_model.VLvar(:,2))];
            %VIP=pls_model.P.^2*Yvar;
            %VIP=VIP./max(VIP);
            %In_error=100-predict.Sucp;
            %Calculate the feature weight from each cycle
            %featureWeight(c,r')=VIP'*(predict.Sucp/100);
            featureWeight(c,r')=predict.Sucp/100;
            freq(c,r')=1;
            %disp(strcat('TestNum:',num2str(length(Test)),' Error:',num2str(In_error),'%'));
        else
            pls_model = PLS(In_bestfeatures(Train,:),class(Train),numM);
            predict = plspred(In_bestfeatures(Test,:),pls_model,class(Test));
            %Variable Importance in Projection (VIP) score for each features 
            %Yvar=[pls_model.VLvar(1,2);diff(pls_model.VLvar(:,2))];
            %VIP=pls_model.P.^2*Yvar;
            %VIP=VIP./max(VIP);
            if strcmp(p.Results.rObj,'MSE')
                %featureWeight(c,r')=VIP'*(1/min(predict.RMSEp));
                featureWeight(c,r')=1/min(predict.RMSEp);
            else
                if ~isnan(max(predict.R2p))%If class(Test) are the same number, predict.R2p would be NaN 
                    %featureWeight(c,r')=VIP'*max(predict.R2p);
                    featureWeight(c,r')=max(predict.R2p);
                end
            end
            freq(c,r')=1;
            %disp(strcat('TestNum:',num2str(length(Test)),' RMSE:',num2str(min(predict.RMSEp)),' R2:',num2str(max(predict.R2p))));
        end
    end
    totalfw=sum(featureWeight);
    totalfq=sum(freq);
    [FeaSort,Rank]=sort(totalfw,'descend');
    [bestnames,bestfeatures,error]=DetermineModel(p.Results.rObj,modelType,class,Rc_features_cz,feanames_c,Rank,num,p.Results.SetFea);
    feaWeight=array2table([FeaSort;error;totalfq(Rank)],'VariableNames',feanames_c(Rank),'RowNames',{'FeaWeight','CV_errors','Frequency'});
    bestnames_M=bestnames;
    if ~isempty(p.Results.feaName_adjust)
        for i=1:length(p.Results.feaName_adjust{1})
            bestnames_M=strrep(bestnames_M,p.Results.feaName_adjust{1}{i},p.Results.feaName_adjust{2}{i});
        end
    end
    %PCA 
    %if p.Results.RefSeqTable
    %    [fwcoeff_s,fscore_s,flatent_s,ftsquared_s,fexplained_s,fmu_s] = pca(bestfeatures);
    %    fc3_s = fwcoeff_s(:,1:3);
    %    fcoefforth_s = inv(diag(std(bestfeatures)))*fwcoeff_s;
    %    fI_s = fc3_s'*fc3_s;
    %    D1=p.Results.PCAD1;D2=p.Results.PCAD2;
    %    pca_plot(bestnames_M,fscore_s,fcoefforth_s,fexplained_s,D1,D2,class,class_name,p.Results.PCA_xr,p.Results.PCA_yr,strcat(SaveName,'_pca'),10,alt_class_name);
    %    gridPlot(class,class_name,fscore_s,strcat(SaveName,'_pca_grid'),0,'PCA',alt_class_name);
    %end
    %PLSDA
    if num<3; num=3; end%for drawing purpose (3D)
    if length(bestnames_M)<num; num=length(bestnames_M); end%for drawing purpose (3D)
    if strcmp(modelType,'da') 
        pls_model_s = PLS(bestfeatures,class,num,'da');
        pls_cv_s = PLSCV(bestfeatures,class,num,'da');%Leave-one-out cross-validation for PLS regression or discriminant analysis
    else
        pls_model_s = PLS(bestfeatures,class,num);
        pls_cv_s = PLSCV(bestfeatures,class,num);%Leave-one-out cross-validation for PLS regression or discriminant analysis
    end
    %disp('Feature weights at each PC:');
    %disp(array2table(pls_model_s.W,'RowNames',bestnames_M));
    %disp('Feature loading at each PC:');
    %disp(array2table(pls_model_s.P,'RowNames',bestnames_M));
    
    if p.Results.plotall
        D1=p.Results.PLSDAD1;D2=p.Results.PLSDAD2;
        if p.Results.legendText
            plsda_plot(bestnames_M,pls_model_s,pls_cv_s,D1,D2,class,class_name,p.Results.PLSDA_xr,p.Results.PLSDA_yr,strcat(SaveName,'_pls',modelType),12,alt_class_name,...
                'legendScale',p.Results.legendScale,'legendText',p.Results.legendText);
        else
            plsda_plot(bestnames_M,pls_model_s,pls_cv_s,D1,D2,class,class_name,p.Results.PLSDA_xr,p.Results.PLSDA_yr,strcat(SaveName,'_pls',modelType),12,alt_class_name);
        end
        gridPlot(class,class_name,pls_model_s.T,strcat(SaveName,'_plsda_grid'),1,'PLSDA',alt_class_name);
    end
end
function alpha=EvaAlpha(Rc_features_cz,class,lasscv_fold)
    alphaTest=0.1:0.1:1;
    MinMSE=zeros(1,length(alphaTest));
    c=1;
    for i=alphaTest
        [b fitinfo]=lasso(Rc_features_cz,class,'CV',lasscv_fold,'Alpha',i);
        MinMSE(c)=fitinfo.MSE(fitinfo.IndexMinMSE);
        c=c+1;
    end
    [a,b]=min(MinMSE);
    alpha=alphaTest(b);
end
function bestmodel=estimateBestmodel(rObj,b,fitinfo,Rc_features_cz,class,set,modelType,num)
    if max(fitinfo.DF)<=num
        %No other models qualified, because not enough features to biold
        %model with required LV numbers
        bestmodel=fitinfo.IndexMinMSE;
    else
        LASSObestmodel=fitinfo.IndexMinMSE;
        bmSE=fitinfo.SE(LASSObestmodel);
        bmMSE=fitinfo.MSE(LASSObestmodel);
        testModel=[];tmp=1;
        while isempty(testModel)%In centain case, it could have no features selected and need to increase bmSE
            SDrange=find(fitinfo.MSE <= bmMSE+(bmSE*tmp/2) & fitinfo.MSE >= bmMSE);%only consider the models whose MSE in the minimal MSE+(SE/2)
            testModel=intersect(find(sum(abs(b)>0)>=num),SDrange);%And its # of features more than PLS dimensions
            tmp=tmp+1;
        end
        error=zeros(1,length(testModel));
        c=1;
        for e=testModel
            % identify columns containing the best features
            col_bestmodel=b(:,e);
            r = find(abs(col_bestmodel)>0);
            bestfeatures=Rc_features_cz(:,r');
            if strcmp(modelType,'da')
                pls_cv = PLSCV(bestfeatures,class,num,'da');
                error(c)=100-max(pls_cv.Succv);
            else
                pls_cv = PLSCV(bestfeatures,class,num);
                if strcmp(rObj,'MSE')
                    error(c)=min(pls_cv.RMSEcv);
                else
                    error(c)=1-max(pls_cv.R2cv);
                end
            end
            c=c+1;
        end
        if strcmp(set,'min')
            bestmodel=testModel(max(find(ismember(error,min(error)))));
        elseif strcmp(set,'max')
            bestmodel=testModel(min(find(ismember(error,min(error)))));
        end
        %disp(strcat('Original BestModel:',num2str(LASSObestmodel)));
        %disp([testModel;error]);
        %disp(strcat('Selected BestModel:',num2str(bestmodel)));
    end
end
function [bestnames,bestfeatures,error]=DetermineModel(rObj,modelType,class,Rc_features_cz,feanames_c,Rank,num,SetFea)
    %Given a rank of the features, define the model by adding the
    %features and find the minimal errors
    Ind=Rank(1:num-1);
    error=zeros(1,length(Rank));error(1:num-1)=nan;
    for c=num:length(Rank)
        Ind=[Ind,Rank(c)];
        if strcmp(modelType,'da') 
            pls_cv = PLSCV(Rc_features_cz(:,Ind),class,num,'da');
            error(c)=100-max(pls_cv.Succv);
        else
            pls_cv = PLSCV(Rc_features_cz(:,Ind),class,num);
            if strcmp(rObj,'MSE')
                error(c)=min(pls_cv.RMSEcv);
            else
                error(c)=1-max(pls_cv.R2cv);
            end
        end
    end
    figure;
    plot(num:length(Rank),error(num:length(Rank)),'-o','LineWidth',2,'MarkerSize',3,'MarkerEdgeColor','b','MarkerFaceColor','red');
    
    [minValue,minInd1]=min(error(num:length(Rank)));
    if strcmp(SetFea,'max')
        minInd=max(find(ismember(error,minValue)));
    else
        minInd=minInd1+num-1;
    end
    bestnames=feanames_c(sort(Rank(1:minInd)));
    bestfeatures=Rc_features_cz(:,sort(Rank(1:minInd)));
end