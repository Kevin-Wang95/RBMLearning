function Z = AIS(WB, aA, aB, bA, bB, v0, K)
%AIS 此处显示有关此函数的摘要
%   v is a col vector
%   a* and b* are row vectors
%   WA is zeros matrix

    WA = zeros(size(WB));
    beta = sort(rand(1,K));
    beta(1) = 0; beta(end) = 1;
    v = v0;
    Z = 0;
    Z = Z+sum(log(1+exp(aA)));
    Z = Z+sum(log(1+exp(bA)));
   
    for k = 1:length(beta)-1
        Pk = (1-beta(k))*(bA*v);
        Pk = Pk + sum(log(1+exp((1-beta(k))*(WA'*v+aA'))));
        Pk = Pk + (beta(k)*(bB*v));
        Pk = Pk + sum(log(1+exp(beta(k)*(WB'*v+aB'))));
        Pk = exp(Pk);
        Pk_1 = (1-beta(k+1))*(bA*v);
        Pk_1 = Pk_1 + sum(log(1+exp((1-beta(k+1))*(WA'*v+aA'))));
        Pk_1 = Pk_1 + (beta(k+1)*(bB*v));
        Pk_1 = Pk_1 + sum(log(1+exp(beta(k+1)*(WB'*v+aB'))));
        Pk_1 = exp(Pk_1);
        Z = Z + log(Pk_1/Pk);
        
        hA = zeros(1,size(WA,2));
        pA = sigm((1-beta(k))*(WA'*v+aA'));            
        u = rand(size(pA));
        hA(pA>u) = 1;

        hB = zeros(1,size(WB,2));
        pB = sigm(beta(k)*(WB'*v+aB'));            
        u = rand(size(pB));
        hB(pB>u) = 1;
        
        pV = sigm((1-beta(k))*(WA*hA'+bA')+beta(k)*(WB*hB'+bB'));
        u = rand(size(pV));
        v(pV<=u) = 0; v(pV>u) = 1;
    end
end
