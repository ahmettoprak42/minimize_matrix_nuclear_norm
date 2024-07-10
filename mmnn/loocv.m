clc, clear,

%%
% A = importdata('miRNA_disease.txt');
% rna = importdata('miRNASimilarity.txt');
% dss = importdata('diseaseSemantic1.txt');
% 
% interaction = A;
% disSim = rna;
% microSim = dss;

%%
interaction_ori = interaction; 
n = 1;
m = 1;

%%
[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number
[dn,dr] = size(interaction);

maxiter = 300;
alpha = 1;
beta = 1;
tol1 = 2*1e-3;
tol2 = 1*1e-5;  

[GD,GM] = GIPSim(interaction,1,1);
HD=zeros(nd,nd);
HM=zeros(nc,nc);

%%
for i=1:nd
    for j=1:nd
        if disSim(i,j)>0
            HD(i,j)=(disSim(i,j)+GD(i,j))/2;
        else
            HD(i,j)=GD(i,j);
        end
    end
end

%%
for i=1:nc
    for j=1:nc
        if microSim(i,j)>0
            HM(i,j)=(microSim(i,j)+GM(i,j))/2;
        else
            HM(i,j)=GM(i,j);
        end
    end
end

%%
T = [HM, interaction'; interaction, HD];
[t1, t2] = size(T);
trIndex = double(T ~= 0);
[WW,iter] = DRMNN(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
F_ori = WW((t1-dn+1) : t1, 1 : dr);

%% LOOCV
index = find(interaction_ori==1);

for u=1:length(index)
    u
    interaction(index(u))=0;
    [GD,GM] = GIPSim(interaction,1,1);
    HD=zeros(nd,nd);
    HM=zeros(nc,nc);
    for i=1:nd
        for j=1:nd
            if disSim(i,j)>0
                HD(i,j)=(disSim(i,j)+GD(i,j))/2;
            else
                HD(i,j)=GD(i,j);
            end
        end
    end
    for i=1:nc
        for j=1:nc
            if microSim(i,j)>0
                HM(i,j)=(microSim(i,j)+GM(i,j))/2;
            else
                HM(i,j)=GM(i,j);
            end
        end
    end
    T = [HM, interaction'; interaction, HD];
    [t1, t2] = size(T);
    trIndex = double(T ~= 0);
    [WW,iter] = DRMNN(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);
    F = WW((t1-dn+1) : t1, 1 : dr);
    F_ori(index(u))=F(index(u));
    interaction = interaction_ori;   
end

%%
pre_label_score = F_ori(:);

% kaydetmek için
score_loocv = pre_label_score;
save score_loocv

%% ROC eğrisi
label_y = interaction_ori(:);
auc = roc_1(pre_label_score,label_y,'red');

%% ACC,PRE,SEN,F1_score,MCC değerleri
% [ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.99' );

%%
% [AUC,AUPR,Acc,Sen,Spe,Pre] = ROCcompute(F_ori,interaction,1); 
% % aupr aynısı
% aupr = pr_cure(pre_label_score,label_y,'red');
