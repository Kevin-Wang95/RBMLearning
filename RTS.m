function Zout = RTS(W, a, bA, bB, T)
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
    ZA = sum(log(1+exp(aA)));
    ZA = ZA + sum(log(1+exp(bA)));
    beta = linspace(0,1,K);
%     figure;
    i = 1;
    while max(abs(1/K-c)) >= T/K | i==1
%         beta = sort(rand(1,K));
%         beta(1) = 0;beta(end) = 1;
        beta_s = beta(randi([1,length(beta)]));
        v = zeros(size(W,1),1); %   v is a col vector
        c = zeros(1,K);
        for n = 1:N
            hA = zeros(1,size(WA,2));
            pA = logsig((1-beta_s)*(WA'*v+aA'));            
            u = rand(size(pA));
            hA(pA>u) = 1;

            hB = zeros(1,size(WB,2));
            pB = logsig(beta_s*(WB'*v+aB'));            
            u = rand(size(pB));
            hB(pB>u) = 1;

            pV = logsig((1-beta_s)*(WA*hA'+bA')+beta_s*(WB*hB'+bB'));
            u = rand(size(pV));
            v(pV<=u) = 0; v(pV>u) = 1;
            
            Q = zeros(1,K);
            for k = 1:K
                Q(k) = (1-beta(k))*(bA*v);
                Q(k) = Q(k) + sum(log(1+exp((1-beta(k))*(WA'*v+aA'))));
                Q(k) = Q(k) + beta(k)*(bB*v);
                Q(k) = Q(k) + sum(log(1+exp(beta(k)*(WB'*v+aB'))));
                Q(k) = exp(Q(k))/Z(k)/exp(ZA)^(1-beta(k));
            end
            Q = Q./sum(Q);
            u = rand(1,1);
            s = 1;
            while u > 0
                u = u - Q(s);
                s = s + 1; 
            end
            beta_s = beta(s-1);
            
            c = c + Q/N;
        end
        max(abs(1/K-c))
%         scatter(i,max(abs(1/K-c))); hold on;
        i = i+1;
        Z = Z.*c/c(1);
    end
    Zout = log(Z(end));
end

