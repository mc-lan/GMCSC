function Z = solveZ_block(Y,M,J,W,E,Y1,Y2,Y3,mu)
numOfView = length(Y);
numOfSamples = size(Y{1},2);

for v = 1 : numOfView
    A{v} = Y{v}-E{v}+Y1{v}/mu;
end
B = J-Y2/mu;
C = W-Y3/mu;

%common part
Z = zeros(numOfSamples,numOfSamples);
for v = 1:numOfView
    YtY{v} = Y{v}'*Y{v};
end
index = ones(numOfSamples,1);
for v = 1 :numOfView
    index = index.*diag(M{v});
end
index_common = find(index==1);
for v = 1 :numOfView
    Ac{v} = A{v}(:,index_common);
end
Bc = B(:,index_common);
Cc = C(:,index_common);
temp1 = 0;
temp2 = 0;
for v = 1 :numOfView
    temp1 = temp1 + YtY{v};
    temp2 = temp2 + Y{v}'*Ac{v};
end
Zc = (temp1+2*eye(numOfSamples,numOfSamples))\(temp2+Bc+Cc);
Z(:,index_common) = Zc;
%single part
for v = 1 : numOfView
    index1 = find(diag(M{v})==1);
    index_single{v} = setdiff(index1,index_common);
    As{v} = A{v}(:,index_single{v});
    Bs = B(:,index_single{v});
    Cs = C(:,index_single{v});
    Zs = (YtY{v}+2*eye(numOfSamples,numOfSamples))\(Y{v}'*As{v}+Bs+Cs);
    Z(:,index_single{v}) = Zs;
end
end