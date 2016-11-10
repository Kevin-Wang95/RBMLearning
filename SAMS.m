function Zout = SAMS(W, a, bA, bB)
%SAMS 此处显示有关此函数的摘要
%   此处显示详细说明
%   a and b are row vectors
    N = 50000;
    alpha = 0.85;
    beta_certain = 0.65;
    K = 100;
    Z = zeros(1,K);
    WA = zeros(size(W));
    WB = W;
    aA = a; aB = a;
%     bA = b; bB = b;
    ZA = 1;
    for i = 1:length(aA)
        ZA = ZA*(1+exp(aA(i)));
    end
    for i = 1:length(bA)
        ZA = ZA*(1+exp(bA(i)));
    end
    beta = sort(rand(1,K));
    beta(1) = 0;beta(end) = 1;
    beta_s = beta(randi([1,length(beta)]));

    for n = 1:N
%         beta = sort(rand(1,K));
%         beta(1) = 0;beta(end) = 1;
%         beta = linspace(0,1,K);
        v = zeros(size(W,1),1); %   v is a col vector
        pA = zeros(1,size(WA,2));
        hA = zeros(1,size(WA,2));
        for j = 1:size(WA,2)
            pA(j) = sigm((1-beta_s)*(WA(:,j)'*v+aA(j)));            
            u = rand(1,1);
            if pA(j) > u
                hA(j) = 1;
            end
        end
        pB = zeros(1,size(WB,2));
        hB = zeros(1,size(WB,2));
        for j = 1:size(WB,2)
            pB(j) = sigm(beta_s*(WB(:,j)'*v+aB(j)));            
            u = rand(1,1);
            if pB(j) > u
                hB(j) = 1;
            end
        end
        pV = zeros(1,size(WB,1));
        for i = 1:length(v)
            pV(i) = sigm((1-beta_s)*(WA(i,:)*hA'+bA(i))+beta_s*(WB(i,:)*hB'+bB(i)));
            u = rand(1,1);
            if pV(i) > u
                v(i) = 1;
            else
                v(i) = 0;
            end
        end                    
        
        Q = zeros(1,K);
        for k = 1:K
            Q(k) = exp((1-beta(k))*(bA*v));
            for j = 1:length(aA)
                Q(k) = Q(k) * (1+exp((1-beta(k))*(WA(:,j)'*v+aA(j))));
            end
            Q(k) = Q(k) * exp(beta(k)*(bB*v));
            for j = 1:length(aB)
                Q(k) = Q(k) * (1+exp(beta(k)*(WB(:,j)'*v+aB(j))));
            end
            Q(k) = Q(k)/exp(Z(k))/ZA^(1-beta(k));
        end
        Q = Q./sum(Q);
        s = find(beta==beta_s);
        if s == 1
            beta_s = beta(2);
        else
            if s == length(beta)
                beta_s = beta(end-1);
            else
                if Q(s-1)< Q(s+1)
                    beta_s = beta(s-1);
                else
                    beta_s = beta(s+1);
                end
            end
        end
        if n < floor(alpha*N)
            t = min(1/K, n^(-beta_certain));
        else
            t = min(1/K, (n - floor(alpha*N) + floor(alpha*N)^beta_certain)^-1);
        end
        
        for k = 1:K
            Z(k) = Z(k) + t*Q(k)/(1/K);
        end
        
        Z = Z - Z(1);
    end
    Zout = Z(end);

end

