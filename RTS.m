function Zout = RTS(W, a, bA, bB)
%RTS 此处显示有关此函数的摘要
%   a and b are row vectors
    aA = a;aB = a;
%     bA = b;bB = b;
    WA = zeros(size(W));
    WB = W;
    K = 100;
    N = 500;
    Z = ones(1,K);
    c = 0;
    ZA = 1;
    for i = 1:length(aA)
        ZA = ZA*(1+exp(aA(i)));
    end
    for i = 1:length(bA)
        ZA = ZA*(1+exp(bA(i)));
    end
    while max(abs(1/K-c)) > 0.2/K
%         beta = sort(rand(1,K));
%         beta(1) = 0;beta(end) = 1;
        beta = linspace(0,1,K);
        beta_s = beta(randi([1,length(beta)]));
        v = zeros(size(W,1),1); %   v is a col vector
        c = zeros(1,K);
        for n = 1:N
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
                Q(k) = Q(k)/Z(k)/ZA^(1-beta(k));
            end
            Q = Q./sum(Q);
            u = rand(1,1);
            s = 1;
            while u > 0
                u = u - Q(s);
                s = s + 1; 
            end
            beta_s = beta(s-1);
            
            for k = 1:K
                c(k) = c(k) + Q(k)/N;
            end
        end
%         max(abs(1/K-c))
        for k = 2:K
            Z(k) = Z(k)*c(k)/c(1);
        end
    end
    Zout = log(Z(end));
end

