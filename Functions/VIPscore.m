%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%Variable Importance in Projection (VIP) score
%
%Interpreting a PLS model with many components and a multitude of responses can be a
%complex task. A parameter which summarizes the importance of the X-variables, both
%for the X- and Y-models, is called the variable influence on projection, VIP. VIP is a
%weighted sum of squares of the PLS weights, w*, taking into account the amount of
%explained Y-variance in each dimension. Its attraction lies in its intrinsic parsimony;
%for a given model and problem there will always be only one VIP-vector, summarizing
%all components and Y-variables. One can compare the VIP of one term to the others.
%Terms with large VIP, larger than 1, are the most relevant for explaining Y. 
%The VIP values reflect the importance of terms in the model both with respect to Y, i.e.
%its correlation to all the responses, and with respect to X (the projection). With
%designed data, i.e. close to orthogonal X, the VIP values mainly reflect the correlation
%of the terms to all the responses. VIP values are computed, by default, from all
%extracted components. 

%dcont:the direction of VIP determine by 1]correlation 2]loading score

function [VIP,YXcorr,FeaName_sort,Fea_ind,colorMap]=VIPscore(Model,FeaName,Feaprofile,outcome,dcont,fileName)
% o Model: The model returned from PLS 
% o FeaName: The feature names used in PLS model (1 X vl)
% o Feaprofile: The corresponding feature profiles (samples X vl)
% o outcome: the outcome(samples X 1)
    
    %based on all components 
    %Yvar=[Model.VLvar(1,2);diff(Model.VLvar(:,2))];
    %VIP=Model.W.^2*Yvar;
    %VIP=Model.P.^2*Yvar;
    
    %based on the frist two components 
    Yvar=[Model.VLvar(1,2);diff(Model.VLvar(1:2,2))];
    VIP=Model.W(:,1:2).^2*Yvar;
    %VIP=Model.P(:,1:2).^2*Yvar;
    
    %YXcorr=corr(Feaprofile,outcome,'type','Spearman');
    YXcorr=corr(Feaprofile,outcome);
    %sort
    [tmp,ind]=sort(VIP,'ascend');
    VIPorder=VIP(ind);
    YXcorr_order=YXcorr(ind);
    loadingOrder=Model.P(ind,1);%First dim
    
    cmap = redblue(1024);
    colorMap=zeros(length(FeaName),3);
    figure;
    if strcmp(dcont,'none')
        barh(1:length(VIPorder),VIPorder);
    elseif strcmp(dcont,'correlation')
        for i=1:length(VIPorder)
            if YXcorr_order(i)> 0%Positive correlate to outcome
                col=cmap(round(((VIPorder(i)+max(VIPorder))/(max(VIPorder)*2))*1024),:);
                h = barh(i,VIPorder(i));
            else%Negative correlate to outcome
                col=cmap(round(((max(VIPorder)-VIPorder(i))/(max(VIPorder)*2))*1024)+1,:);
                h = barh(i,-VIPorder(i));
            end
            colorMap(ind(i),:)=col;
            set(h, 'FaceColor', col);
            hold on;
        end
    elseif strcmp(dcont,'loading')
        for i=1:length(VIPorder)
            if loadingOrder(i)> 0%Positive Loading score
                col=cmap(round(((VIPorder(i)+max(VIPorder))/(max(VIPorder)*2))*1024),:);
                h = barh(i,VIPorder(i));
            else%Negative Loading score
                col=cmap(round(((max(VIPorder)-VIPorder(i))/(max(VIPorder)*2))*1024)+1,:);
                h = barh(i,-VIPorder(i));
            end
            colorMap(ind(i),:)=col;
            set(h, 'FaceColor', col);
            hold on;
        end
        
    end
    Labels = strrep(FeaName(ind),'_','\_');
    set(gca, 'YTick', 1:length(FeaName), 'YTickLabel', Labels);
    title('VIP scores');
    set(gcf,'PaperOrientation','landscape','PaperUnits','normalized','PaperPosition', [0 0 1 1],'PaperType','uslegal');
    print(gcf, '-dpdf', strcat(fileName,'.pdf'));
    FeaName_sort=fliplr(FeaName(ind));
    Fea_ind=flip(ind);%The index of fea sort by VIP
end
