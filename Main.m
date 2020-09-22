%Preacquisition : MATLAB Bioinformatics Toolbox.

%The variables saved in Data.mat
%Raz_OutTable: The outcomes
%Raz_FeaTable: All the assay measurements used to build the multivariate model
%Raz_ZFeaTable: Zscore of all the measurements

clear;
load('Data.mat');
SaveFolder='Results';
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%PCA: considering all time points
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rT_index=logical(~isnan(table2array((Raz_OutTable(:,'reboundtime')))));%rT_index(3)=0;
OutTable=[];
InputTable=[];
List={'_v1$','_v6$','_v9$','_v12$'};
for i=1:length(List)
    each=List{i};
    feaIndex=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames, each)));
    
    tmp=Raz_OutTable(rT_index,{'reboundtime_binary','reboundtime_binary_name','reboundtime','reboundtime_name'});
    tmp.tp=repmat(i,size(tmp,1),1);
    tmp.tpName=repmat({regexprep(each,'_|\$','')},size(tmp,1),1);
    tmp.Pat=tmp.Properties.RowNames;
    tmp.Properties.RowNames=strcat(tmp.Properties.RowNames,'_',num2str(i));
    OutTable=vertcat(OutTable,tmp);
    
    tmp=Raz_ZFeaTable(rT_index,feaIndex);
    tmp.Properties.RowNames=strcat(tmp.Properties.RowNames,'_',num2str(i));
    tmp.Properties.VariableNames=strrep(tmp.Properties.VariableNames,regexprep(each,'\$',''),'');
    InputTable=vertcat(InputTable,tmp);
end
Input=table2array(InputTable);%Input=zscore(table2array(InputTable));

[fwcoeff,fscore,flatent,ftsquared,fexplained] = pca(Input);
fcoefforth = inv(diag(std(Input)))*fwcoeff;

Name=strrep(strrep(strrep(InputTable.Properties.VariableNames,'pp','('),'qq',')'),'ww','-');
pca_df(Name,fscore,fcoefforth,fexplained,1,3,table2array(OutTable(:,'reboundtime')),table2array(OutTable(:,'reboundtime_name')),...
    table2array(OutTable(:,'tp')),table2array(OutTable(:,'tpName')),0.1,0.1,strcat(SaveFolder,'/ProfileDynamics_normFirst_Ras'),10,...
    strcat(table2array(OutTable(:,'Pat')),':',table2array(OutTable(:,'reboundtime_name'))),{},'threeD',0,'plot_errorElli',0,'rescaleX',1);  

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build the prediction models based on all time points 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Considering all measurements
rT_index=find(~isnan(table2array((Raz_OutTable(:,'reboundtime')))));
feaIndex=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames,'_v1$|_v6$|_v9$|_v12$|_d_v6v1$|_d_v9v1$|_d_v12v1$')));

Cxlab='Rebound time (days)';Cylab='Prediction';legendScale={'10','20','30','40','50'};legendText='Rebound time (days)';
%PLSR modeling 
[Rbf_r,Rplsr,RfW_r]=lassoElasticVis('r',zscore(table2array(Raz_OutTable(rT_index,'reboundtime'))),{},[],{},...
table2array(Raz_ZFeaTable(rT_index,feaIndex)),Raz_ZFeaTable(rT_index,feaIndex).Properties.VariableNames,strcat(SaveFolder,'/Raz_allTime'),...
'ElasticLASSO','LASSO','Kfold',5,'plsVL',3);
%Visualization of the model
ModelEvaluation1('r',Rbf_r,Raz_FeaTable(rT_index,:),Raz_ZFeaTable(rT_index,:),table2array(Raz_OutTable(rT_index,'reboundtime')),...
    table2array(Raz_OutTable(rT_index,'reboundtime_name')),{},strcat(SaveFolder,'/Raz_allTime_test'),...
    'PLSDAD1',1,'PLSDAD2',2,'PLSDA_xr',10,'PLSDA_yr',2,'Cxlab',Cxlab,'Cylab',Cylab,'KMxlab','Days after ATI','KMylab','Unrebound percentage',...
    'feaName_adjust',{{'pp','qq','ww'},{'(',')','-'}},'legendScale',legendScale,'legendText',legendText,'plsVL',3);   
%Permutation test
ModelEvaluation2('r',200,zscore(table2array(Raz_OutTable(rT_index,'reboundtime'))),Raz_ZFeaTable(rT_index,feaIndex),...
        Rbf_r,strcat(SaveFolder,'/Raz_allTime_test_Null'),'cvIteration',200,'plsVL',3);   

%Only considering antibody function measurements
feaIndex_func=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames, 'ADCD|ADCP|CD107|MIP1b|IFNy|ADCC|ADNP')));
feaIndexFunc=intersect(feaIndex,feaIndex_func);
Cxlab='Rebound time (days)';Cylab='Prediction';legendScale={'10','20','30','40','50'};legendText='Rebound time (days)';
[Rbf_func_r,Rplsr_func,RfW_func_r]=lassoElasticVis('r',zscore(table2array(Raz_OutTable(rT_index,'reboundtime'))),{},[],{},...
table2array(Raz_ZFeaTable(rT_index,feaIndexFunc)),Raz_ZFeaTable(rT_index,feaIndexFunc).Properties.VariableNames,strcat(SaveFolder,'/Raz_allTime_func'),...
'ElasticLASSO','LASSO','Kfold',5,'plsVL',3);
ModelEvaluation1('r',Rbf_func_r,Raz_FeaTable(rT_index,:),Raz_ZFeaTable(rT_index,:),table2array(Raz_OutTable(rT_index,'reboundtime')),...
    table2array(Raz_OutTable(rT_index,'reboundtime_name')),{},strcat(SaveFolder,'/Raz_allTime_func'),...
    'PLSDAD1',1,'PLSDAD2',2,'PLSDA_xr',5,'PLSDA_yr',2,'Cxlab',Cxlab,'Cylab',Cylab,'KMxlab','Days after ATI','KMylab','Unrebound percentage',...
    'feaName_adjust',{{'pp','qq','ww'},{'(',')','-'}},'legendScale',legendScale,'legendText',legendText,'plsVL',3);   
ModelEvaluation2('r',200,zscore(table2array(Raz_OutTable(rT_index,'reboundtime'))),Raz_ZFeaTable(rT_index,feaIndexFunc),...
        Rbf_func_r,strcat(SaveFolder,'/Raz_allTime_func_Null'),'cvIteration',200,'plsVL',3);   
    

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Build the prediction models based on rebount time on either
%functions or glycans or all
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rT_index=find(~isnan(table2array((Raz_OutTable(:,'reboundtime')))));
Cxlab='Rebound time group';Cylab='Prediction';legendScale={'10','20','30','40','50'};legendText='Rebound time (days)';
rebT=log2(table2array(Raz_OutTable(rT_index,'reboundtime')));

for type={'function','glycan','All'}
    type=type{:};
    for each={'_b$','_auc$','all','_v1$','_v6$','_v9$','_v12$','_d_v6v1$','_d_v9v1$','_d_v12v1$'}
        each=each{:};
        model=regexprep(each,'\_d\_|\_|\$','');disp(strcat('Build the model:',model));
        if strcmp(each,'all')
            feaIndex1=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames, '_v1$|_v6$|_v9$|_v12$|_b$|_auc$|d_v6v1$|d_v9v1$|d_v12v1$')));
        else
            feaIndex1=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames, each)));
        end
        feaIndex_func=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames, 'ADCD|ADCP|CD107|MIP1b|IFNy|ADCC|ADNP')));
        feaIndex_glycan=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames,...
        '^G0|^G1|^G2|Sialic|Fucose|Bisecting|G0F|G0ppFBqq|G1_1|G1F_1|G1ppwwqq|G1ppSqq|G1ppFqq|G1ppBqq|G2ppSqq|G2ppS1qq|G2ppS2qq|G2S1|G2Sppwwqq|G2ppFxqq|G2ppFwwSqq|G2ppFwwBqq|G2F|G2ppBxqq|G2ppwwBqq|G2S1F|G2FB')));

        if strcmp(type,'glycan')
            feaIndex=intersect(feaIndex1,feaIndex_glycan);
        elseif strcmp(type,'function')
            feaIndex=intersect(feaIndex1,feaIndex_func);
        else
            feaIndex=feaIndex1;
        end

        [Rbf_r,Rplsr,RfW_r]=lassoElasticVis('r',zscore(rebT),{},[],{},...
        table2array(Raz_ZFeaTable(rT_index,feaIndex)),Raz_ZFeaTable(rT_index,feaIndex).Properties.VariableNames,strcat(SaveFolder,'/Raz_',type,'_',model),...
        'ElasticLASSO','LASSO','Kfold',5,'plsVL',3);
        ModelEvaluation1('r',Rbf_r,Raz_FeaTable(rT_index,:),Raz_ZFeaTable(rT_index,:),rebT,...
        table2array(Raz_OutTable(rT_index,'reboundtime_name')),{},strcat(SaveFolder,'/Raz_',type,'_',model),...
        'PLSDAD1',1,'PLSDAD2',2,'PLSDA_xr',5,'PLSDA_yr',2,'Cxlab',Cxlab,'Cylab',Cylab,'KMxlab','Days after ATI','KMylab','Unrebound percentage',...
        'feaName_adjust',{{'pp','qq','ww'},{'(',')','-'}},'legendScale',legendScale,'legendText',legendText,'plsVL',3);
        ModelEvaluation2('r',200,zscore(rebT),Raz_ZFeaTable(rT_index,feaIndex),Rbf_r,strcat(SaveFolder,'/Raz_',type,'_',model,'_Null'),'cvIteration',200,'plsVL',3);
        save(strcat(SaveFolder,'/Model_',type,'_',strrep(each,'$','')),'Rbf_r','Rplsr','RfW_r');
        close all;
    end
end


%%
%Build the prediction models based on DNA
Cxlab='delta\_DNA';Cylab='Prediction';legendScale={'-1','0','1'};legendText='delta\_DNA';
DNA=zscore(table2array(Raz_OutTable(:,'totDNA_b')));

for type={'glycan','function','all'}
    type=type{:};
    for each={'_v1$','_v6$','_v9$'}
        each=each{:};
        model=regexprep(each,'\_d\_|\_|\$','');disp(strcat('Build the model:',model));
        feaIndex1=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames, each)));
        feaIndex_func=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames, 'ADCD|ADCP|CD107|MIP1b|IFNy|ADCC|ADNP')));
        feaIndex_glycan=find(~cellfun(@isempty,regexp(Raz_ZFeaTable.Properties.VariableNames,...
        '^G0|^G1|^G2|Sialic|Fucose|Bisecting|G0F|G0ppFBqq|G1_1|G1F_1|G1ppwwqq|G1ppSqq|G1ppFqq|G1ppBqq|G2ppSqq|G2ppS1qq|G2ppS2qq|G2S1|G2Sppwwqq|G2ppFxqq|G2ppFwwSqq|G2ppFwwBqq|G2F|G2ppBxqq|G2ppwwBqq|G2S1F|G2FB')));

        if strcmp(type,'glycan')
            feaIndex=intersect(feaIndex1,feaIndex_glycan);
        elseif strcmp(type,'function')
            feaIndex=intersect(feaIndex1,feaIndex_func);
        else
            feaIndex=feaIndex1;
        end

        [Rbf_r,Rplsr,RfW_r]=lassoElasticVis('r',DNA,{},[],{},...
        table2array(Raz_ZFeaTable(:,feaIndex)),Raz_ZFeaTable(:,feaIndex).Properties.VariableNames,strcat(SaveFolder,'/Raz_',type,'_',model),...
        'ElasticLASSO','LASSO','Kfold',5,'plsVL',3);
        ModelEvaluation1('r',Rbf_r,Raz_FeaTable,Raz_ZFeaTable,DNA,...
        table2array(Raz_OutTable(:,'totDNA_b_binary_name')),{},strcat(SaveFolder,'/Raz_',type,'_',model),...
        'PLSDAD1',1,'PLSDAD2',2,'PLSDA_xr',10,'PLSDA_yr',10,'Cxlab',Cxlab,'Cylab',Cylab,'KMxlab','DNA','KMylab','Unrebound percentage',...
        'feaName_adjust',{{'pp','qq','ww'},{'(',')','-'}},'legendScale',legendScale,'legendText',legendText);
        ModelEvaluation2('r',200,DNA,Raz_ZFeaTable(:,feaIndex),Rbf_r,strcat(SaveFolder,'/Raz_',type,'_',model,'_Null'),'cvIteration',200);
        save(strcat(SaveFolder,'/Model_DNA_',type,'_',strrep(each,'$','')),'Rbf_r','Rplsr','RfW_r');
        close all;
    end
end

