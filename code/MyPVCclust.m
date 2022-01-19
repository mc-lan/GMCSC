function [Z,P,ACC,NMI,Purity,H] = MyPVCclust(X,Xc,Xs,M,numClust,truth,option)
%%
alpha = option.alpha;
beta = option.beta;
lambda = option.lambda;
k = option.latentdim;
if (min(truth)==0)
    truth = truth + 1;
end
num_views = length(X);
n = size(X{1},2);
nc = size(Xc{1},2);
Hc = ones(nc,nc)*(1/nc)*(-1) + eye(nc);
for v = 1:num_views
    d(v) = size(X{v},1);
    ns(v) = size(Xs{v},2);
    Hs{v} = ones(ns(v),ns(v))*(1/ns(v))*(-1) + eye(ns(v));
    if size(X{v},1)<100
        k = size(X{v},1);
    else
        k = option.latentdim;
    end
    [P{v},~] = eigs(Xs{v}*Hs{v}*Xs{v}',k,'LR');
end

T = min(ceil(n/numClust),30);
iter_out_max = 4;
H = [];
for iter_out = 1:iter_out_max
    for v = 1:num_views
        Y{v} = P{v}'*X{v};
        H = [H;Y{v}];
    end
    [Z,~] = LRS(X,M,alpha,T,0);
    if iter_out<iter_out_max
        [P] = solveP(Z,P,X,Xc,Xs,M,option);
    end
end

% fprintf('running spectral clustering...\n');
kmeans_avg_iter = 5;
for i=1: kmeans_avg_iter
    W = 1/2*(abs(Z)+abs(Z)');
    C = SpectralClustering(W,numClust);
    [~, nmii(i),~] = compute_nmi(truth,C);
    acci(i) = Accuracy(C,truth);
    %[Fi(i),Pi(i),Ri(i)] = compute_f(truth,C);
    %[ARi(i),RIi(i),MIi(i),HIi(i)]=RandIndex(truth,C);
    purityi(i) = compute_Purity(truth, C);
    
end
NMI = mean(nmii);
ACC = mean(acci);
Purity = mean(purityi);
fprintf('nmi: %f(%f)\n', NMI, std(nmii));
end