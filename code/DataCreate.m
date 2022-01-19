function [X,Xs,Xc,M,truth] = DataCreate(dataname,PER,fold)
datasetdir='D:\Program Files\MATLAB\R2020b\bin\lan\PartialMV\Experiment\data\';
pairPortion = 1 - PER;
load(strcat(datasetdir,dataname,'RnSp.mat'));
load(strcat(datasetdir,dataname,'Folds.mat'));
[~,numInst] = size(folds);
num_views = length(X);
for v =1 :num_views
    X{v} = normc(X{v});
end
instanceIdx = folds(fold,:);
commonInst = floor(numInst*pairPortion);  % number of paired instances  have complete views
common_list = sort(instanceIdx(1:commonInst));
singleInst_missing = floor((numInst-commonInst)/num_views);
for v = 1 :num_views
    X_copy{v} = X{v};
    Xc{v} = X{v}(:,common_list);  %common instances
    paritial_list = setdiff(instanceIdx,instanceIdx(1:commonInst));
    X{v}(:,paritial_list) = 0; %single instances with mising 0
    if v<num_views
        single_list{v} = instanceIdx(commonInst+singleInst_missing*(v-1)+1:commonInst+singleInst_missing*v);
        X{v}(:,single_list{v}) = X_copy{v}(:,single_list{v});  %recover single data in each views
        temp = X{v};
        temp(:,find(sum(abs(temp),1)==0))=[];
        Xs{v} = temp;   %single instances only data
        MM{v} = zeros(numInst,1);
        MM{v}([common_list single_list{v}]) = 1;
        M{v} = diag(MM{v});
    else
        single_list{v} = instanceIdx(commonInst+singleInst_missing*(v-1)+1:end);
        X{v}(:,single_list{v}) = X_copy{v}(:,single_list{v});  %recover instances with mising 0
        temp = X{v};
        temp(:,find(sum(abs(temp),1)==0))=[];
        Xs{v} = temp;   %single instances only data
        MM{v} = zeros(numInst,1);
        MM{v}([common_list single_list{v}]) = 1;
        M{v} = diag(MM{v});
    end
end
end