clc, clear,

%%
A = importdata('miRNA_disease.txt');
rna = importdata('miRNASimilarity.txt');
dss = importdata('diseaseSemantic1.txt');

interaction = A;
disSim = rna;
microSim = dss;

%%
interaction_ori = interaction;

[nd,nc] = size(interaction); % nd:diseases number, nc:microbes number
[dn,dr] = size(interaction);

maxiter = 300;
%%
% for alpha = [0.1,1,10,100]
% for beta = [0.1,1,10,100,1000]

alpha = 1;
beta = 1;

%%
tol1 = 2*1e-3;
tol2 = 1*1e-5;  
Kd = 32;
Km = 5;
[GD,GM] = GIPSim(interaction,1,1);
HD = zeros(nd,nd);
HM = zeros(nc,nc);
    
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

%%
[WW,iter] = DRMNN(alpha, beta, T, trIndex, tol1, tol2, maxiter, 0, 1);

%%
F_ori = WW((t1-dn+1) : t1, 1 : dr);
index = find(interaction_ori==1);

%%  5-fold
auc = zeros();
aupr = zeros();

rng(0);
for k = 1:100 %iterasyon sayýsý
    indices = crossvalind('Kfold', length(index), 5);
    for cv = 1:5
        % k
        cv
        index2 = find(cv == indices);
        interaction(index(index2)) = 0;

        %% 
        [GD,GM] = GIPSim(interaction,1,1);
        HD = zeros(nd,nd);
        HM = zeros(nc,nc);
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
        F_ori(index(index2)) = F(index(index2));
        interaction = interaction_ori;
    end
    
    pre_label_score = F_ori(:);
    % resim 1 için
    score_5fold = pre_label_score;
    save score_5fold
    
    label_y = interaction_ori(:);
    auc(k) = roc_1(pre_label_score,label_y,'red');

    % [ACC,PRE,SEN,F1_score,MCC] = myACC_1( F_ori(:),interaction_ori(:),'sp0.95' );
    % figure;
    % aupr(k) = pr_cure(pre_label_score,label_y,'red');
end

%%
% [AUC1,AUPR1,Acc,Sen,Spe,Pre] = ROCcompute(F_ori,interaction,1);
% auc_ave = mean(auc);
% auc_std = std(auc);