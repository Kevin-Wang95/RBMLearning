function Z = AIS(WB, aA, aB, bA, bB, v0)
%AIS �˴���ʾ�йش˺�����ժҪ
%   v is a col vector
%   a* and b* are row vectors
%   WA is zeros matrix

    WA = zeros(size(WB));
    beta = sort(rand(1,30000));
    beta(1) = 0; beta(end) = 1;
    v = v0;
    Z = 1;
    for i = 1:length(aA)
        Z = Z*(1+exp(aA(i)));
    end
    for i = 1:length(bA)
        Z = Z*(1+exp(bA(i)));
    end
    
    Z = log(Z);
    
    for k = 1:length(beta)-1
        Pk = exp((1-beta(k))*(bA*v));
        for j = 1:length(aA)
            Pk = Pk * (1+exp((1-beta(k))*(WA(:,j)'*v+aA(j))));
        end
        Pk = Pk * exp(beta(k)*(bB*v));
        for j = 1:length(aB)
            Pk = Pk * (1+exp(beta(k)*(WB(:,j)'*v+aB(j))));
        end
        Pk_1 = exp((1-beta(k+1))*(bA*v));
        for j = 1:length(aA)
            Pk_1 = Pk_1 * (1+exp((1-beta(k+1))*(WA(:,j)'*v+aA(j))));
        end
        Pk_1 = Pk_1 * exp(beta(k+1)*(bB*v));
        for j = 1:length(aB)
            Pk_1 = Pk_1 * (1+exp(beta(k+1)*(WB(:,j)'*v+aB(j))));
        end
        Z = Z + log(Pk_1/Pk);
        pA = zeros(1,size(WA,2));
        hA = zeros(1,size(WA,2));
        for j = 1:size(WA,2)
            pA(j) = sigm((1-beta(k))*(WA(:,j)'*v+aA(j)));            
            u = rand(1,1);
            if pA(j) > u
                hA(j) = 1;
            end
        end
        pB = zeros(1,size(WB,2));
        hB = zeros(1,size(WB,2));
        for j = 1:size(WB,2)
            pB(j) = sigm(beta(k)*(WB(:,j)'*v+aB(j)));            
            u = rand(1,1);
            if pB(j) > u
                hB(j) = 1;
            end
        end
        pV = zeros(1,size(WB,1));
        for i = 1:length(v)
            pV(i) = sigm((1-beta(k))*(WA(i,:)*hA'+bA(i))+beta(k)*(WB(i,:)*hB'+bB(i)));
            u = rand(1,1);
            if pV(i) > u
                v(i) = 1;
            else
                v(i) = 0;
            end
        end        
    end
end
