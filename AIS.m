function Z = AIS( WA, WB, aA, aB, bA, bB)
%AIS 此处显示有关此函数的摘要
%   此处显示详细说明
    beta = sort(rand(1,20000));
    beta(1) = 0; beta(end) = 1;
    
    

end

function y = g(x)
    y = 1/(1+exp(-x));
end