function P = solveP(Z,P,X,Xc,Xs,M,option)
alpha = option.alpha;
beta = option.beta;
lambda = option.lambda;
k = option.latentdim;

num_views = length(X);
n = size(X{1},2);
nc = size(Xc{1},2);
Hc = ones(nc,nc)*(1/nc)*(-1) + eye(nc);
for v = 1:num_views
    [d{v},ns(v)] = size(Xs{v});
    Hs{v} = ones(ns(v),ns(v))*(1/ns(v))*(-1) + eye(ns(v));
end
for v = 1:num_views
    Z_new = Z*M{v};
    Zc{v} = Z(find(sum(Z_new,1)~=0),find(sum(Z_new,1)~=0));
    Cs{v} = Xs{v}*Hs{v}*Xs{v}';
end

for i = 1 : 2
    for v = 1:num_views
        temp = P{v}'*(Xs{v}-Xs{v}*Zc{v});
        D = diag(0.5./sqrt(sum(temp.*temp))+eps);
        D(isnan(D))=0;
        D(isinf(D))=1e5;
        S{v} = (Xs{v}-Xs{v}*Zc{v})*D*(Xs{v}-Xs{v}*Zc{v})';
        
        K_complement = zeros(nc,nc);
        for v1 = 1:num_views
            K_complement =  K_complement + Hc*Xc{v1}'*P{v1}*P{v1}'*Xc{v1}*Hc;
        end
        T = K_complement - Hc*Xc{v}'*P{v}*P{v}'*Xc{v}*Hc;
        A = - alpha*S{v} + beta*Cs{v} + lambda*Xc{v}*T*Xc{v}';
        A = (A+A')/2;
        [V,D] = eig(A);
        [~, ind] = sort(diag(D),'descend');
        if size(X{v},1)<100
            k = size(X{v},1);
        else
            k = option.latentdim;
        end
        P{v} = V(:, ind(1:k));
        P{v} = real(P{v});
    end
end
end