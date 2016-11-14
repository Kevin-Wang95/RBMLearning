function Zout = SAMS(W, a, b)
%SAMS 此处显示有关此函数的摘要
%   此处显示详细说明
%   a and b are row vectors
    N = 60000;
    alpha = 0.25;
    beta_certain = 0.8;
    K = 100;
    Z = zeros(1,K);
    WA = zeros(size(W));
    WB = W;
    aA = a; aB = a;
    bA = b; bB = b;
    ZA = sum(log(1+exp(aA)));
    ZA = ZA + sum(log(1+exp(bA)));
    beta = sort(rand(1,K));
    beta(1) = 0;beta(end) = 1;
    beta_s = beta(randi([1,length(beta)]));
    v = zeros(size(W,1),1); %   v is a col vector

    for n = 1:N
%         beta = sort(rand(1,K));
%         beta(1) = 0;beta(end) = 1;
%         beta = linspace(0,1,K);
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

        Q = sum(log(1+exp((1-beta')*(WA'*v+aA')')),2);
        Q = Q + bB*v;
        Q = Q + sum(log(1+exp(beta'*(WB'*v+aB')')),2);
        Q = exp(Q'-mean(Q'))./exp(Z);
        Q = Q./sum(Q);
        u = rand(1,1);
%         s = 1;
%         while u > 0
%             u = u - Q(s);
%             s = s + 1; 
%         end
%         beta_s = beta(s-1);
                        s = find(beta==beta_s);
            if s == 1
                beta_s = beta(2);
            else
                if s == length(beta)
                    beta_s = beta(end-1);
                else
                    u = rand();
                    if u < Q(s-1)/(Q(s-1)+Q(s+1))
                        beta_s = beta(s-1);
                    else
                        beta_s = beta(s+1);
                    end
                end
            end
        if n < floor(alpha*N)
            t = min(1/K,n^(-beta_certain));
        else
            t = min(1/K, (n - floor(alpha*N) + floor(alpha*N)^beta_certain)^(-1));
        end
%         t = n^-1;
        Z = Z + t*Q/(1/K);
        Z = Z - Z(1);
    end
    Zout = Z(end) + ZA;

end

