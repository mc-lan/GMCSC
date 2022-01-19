function [ Z,E] = LRS(Y,M,alpha,T,DEBUG)
maxIter = 80;
rho = 1.2;
mu = 1e-2;
max_mu = 1e6;
tol = 1e-6;
numOfView = length(Y);
numOfSamples = size(Y{1},2);
%% Initializing optimization variables
% intializing
Z = zeros(numOfSamples, numOfSamples);
W = zeros(numOfSamples, numOfSamples);
J = zeros(numOfSamples, numOfSamples);
Y2 = zeros(numOfSamples, numOfSamples);
Y3 = zeros(numOfSamples, numOfSamples);
for i=1:numOfView    
    E{i} = zeros(size(Y{i},1), numOfSamples);
    Y1{i} = zeros(size(Y{i},1), numOfSamples);
end

%% Start main loop
iter = 0;
if DEBUG
    disp(['initial,rank(Z)=' num2str(rank(Z))]);
end
while iter<maxIter
    iter = iter + 1;
    
    %solving J
    M1 = Z+Y2/mu;
    [U1, S1, V1] = svd((M1+eps),'econ');
    %     try
    %         [V,D] = eig(M);
    %         [U1, S1, V1] = svd((M1+eps),'econ');
    %     catch ME
    %         if (strcmpi(ME.identifier,'MATLAB:eig:NoConvergence'))
    %             [V,D] = eig(M, eye(size(M)));
    %         else
    %             rethrow(ME);
    %         end
    %     end
    
    S1 = diag(S1);
    svp_J = length(find(S1>1/mu));
    if svp_J>=1
        S1 = S1(1:svp_J)-1/mu;
    else
        svp_J = 1;
        S1 = 0;
    end
    J = U1(:, 1:svp_J)*diag(S1)*V1(:, 1:svp_J)';
    
    %solving Z
    Z = solveZ_block(Y,M,J,W,E,Y1,Y2,Y3,mu);
    %Z = solveZ_single(Y,M,J,W,E,Y1,Y2,Y3,mu);
 
    %solving W
    temp_W = Z+Y3/mu;
    for i=1:numOfSamples
        [~,index] = sort(abs(temp_W(:,i)));
        temp_W(index(1:(numOfSamples-T)),i) = 0;
    end
    W = temp_W;
    
    %solving E
    for i=1:numOfView
        temp_E = Y{i}-Y{i}*Z*M{i}+Y1{i}/mu;
        for j = 1 : size(temp_E, 2)
            if norm(temp_E(:,j), 2) > alpha/mu
                E{i}(:,j) =( norm(temp_E(:,j), 2) - alpha/mu)/norm(temp_E(:,j),2)*temp_E(:,j);
            else
                E{i}(:,j) = 0;
            end
        end
    end
    
    %computing function value
    %     a1 = sum(svd(Z));
    %     for i=1:numOfView
    %         for j=1:size(E{i},2)
    %             E_v(i,j) = norm(E{i}(:,j),2);
    %         end
    %         E_v_j(i) = sum(E_v(i,:));
    %     end
    %     a3 = alpha*sum(E_v_j);
    %     VALUE(iter) = a1+a3;
    
    
    recErr1 = zeros(1,numOfView);
    for i=1:numOfView
        dY1{i} = Y{i}-Y{i}*Z*M{i}-E{i};
        recErr1(i) = norm(dY1{i},'inf');
    end
    recerr1 = max(recErr1);
    dY2 =  Z - J;
    recErr2 = norm(dY2,'inf');
    dY3 =  Z - W;
    recErr3 = norm(dY3,'inf');
    recErr = max(recErr3,max(recerr1, recErr2));
    convergenced = recErr <tol||iter>maxIter;
    
    if DEBUG
        if iter==1 || mod(iter,5)==0 || convergenced
            disp(['iter ' num2str(iter) ',mu=' num2str(mu) ...
                ',rank(Z)=' num2str(rank(Z)) ...
                ',recErr=' num2str(recErr)]);
        end
    end
    if convergenced
        break;
    else
        for i=1:numOfView
            Y1{i} = Y1{i}+mu*dY1{i};
        end
        Y2 = Y2 + mu*dY2;
        Y3 = Y3 + mu*dY3;
        mu = min(max_mu, mu*rho);
    end
end
end