function Z = AIS(WB, aA, aB, bA, bB, v0)
%AIS 此处显示有关此函数的摘要
%   v is a col vector
%   a* and b* are row vectors
%   WA is zeros matrix

    WA = zeros(size(WB));
    beta = sort(rand(1,20000));
    beta(1) = 0; beta(end) = 1;
    v = v0;
    Z = 1;
    for i = 1:length(aA)
        Z = Z*(1+exp(aA(i)));
    end
    PA = 1;
    for i = 1:length(bA)
        Z = Z*(1+exp(bA(i)));
        PA = PA/(1+exp(-bA(i)));
    end
    
    Z = log10(Z);
    
    for k = 1:length(beta)
        pA = zeros(1,size(WA,2));
        hA = zeros(1,size(WA,2));
        for j = 1:size(WA,2)
            pA(j) = g((1-beta(k))*(WA(:,j)'*v+aA(j)));            
            u = rand(1,1);
            if pA(j) > u
                hA(j) = 1;
            end
        end
        pB = zeros(1,size(WB,2));
        hB = zeros(1,size(WB,2));
        for j = 1:size(WB,2)
            pB(j) = g(beta(k)*(WB(:,j)'*v+aB(j)));            
            u = rand(1,1);
            if pB(j) > u
                hB(j) = 1;
            end
        end
        pV = zeros(1,size(WB,1));
        for i = 1:length(v)
            pV(i) = g((1-beta(k))*(WA(i,:)*hA'+bA(i))+beta(k)*(WB(i,:)*hB'+bB(i)));
            u = rand(1,1);
            if pV(i) > u
                v(i) = 1;
            else
                v(i) = 0;
            end
        end
        Pk = exp((1-beta(k))*(bA*v));
        for j = 1:length(aA)
            Pk = Pk * (1+exp((1-beta(k))*(WA(:,j)'*v+aA(j))));
        end
        Pk = Pk * exp(beta(k)*(bB*v));
        for j = 1:length(aB)
            Pk = Pk * (1+exp(beta(k)*(WB(:,j)'*v+aB(j))));
        end
        if(k == 1)
            Z = Z + log10(Pk/PA);
        else
            Pk_1 = exp((1-beta(k-1))*(bA*v));
            for j = 1:length(aA)
                Pk_1 = Pk_1 * (1+exp((1-beta(k-1))*(WA(:,j)'*v+aA(j))));
            end
            Pk_1 = Pk_1 * exp(beta(k-1)*(bB*v));
            for j = 1:length(aB)
                Pk_1 = Pk_1 * (1+exp(beta(k-1)*(WB(:,j)'*v+aB(j))));
            end
            Z = Z + log10(Pk/Pk_1);
        end
    end
end

function y = g(x)
    y = 1/(1+exp(-x));
end