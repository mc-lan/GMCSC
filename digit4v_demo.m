clear;
addpath(genpath('measure/'));
addpath(genpath('code/'));
addpath(genpath('misc/'));
datasetdir='data/';
resultdir = 'Results/';

dataname='digit4v';    %alpha=0.3 beta=1e2 lambda=1e1 0.2

% Parameter to set partial/incomplete example ratio
PERi = [0.3];

% Parameters for the model
alpha = [0.3];
option.beta = 1e2;
option.lambda = 1e1;
option.latentdim = 100;
for per = 1 : length(PERi)
    PER = PERi(per);
    option.alpha = alpha*(1+PER);
    nmi_All = [];
    acc_All = [];
    Purity_All = [];
    
    folds = 1:1;
    for i = 1:length(folds)
        fold = folds(i);
        [X,Xs,Xc,M,truth] = DataCreate(dataname,PER,fold);
        num_views = length(X);
        numClust = length(unique(truth));
        [Z,P,acc,nmi,Purity] = MyPVCclust(X,Xc,Xs,M,numClust,truth,option);
        
        nmi_All = [nmi_All nmi];
        acc_All = [acc_All acc];
        Purity_All = [Purity_All Purity];
    end
    Result(per,1) = mean(nmi_All);
    Result(per,2) = std(nmi_All);
    Result(per,3) = mean(acc_All);
    Result(per,4) = std(acc_All);
    Result(per,5) = mean(Purity_All);
    Result(per,6) = std(Purity_All);
end
% save([resultdir,dataname,'_result.mat'],'Result');