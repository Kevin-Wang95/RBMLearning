function Z = AIS( WA, WB, aA, aB, bA, bB)
%AIS �˴���ʾ�йش˺�����ժҪ
%   �˴���ʾ��ϸ˵��
    beta = sort(rand(1,20000));
    beta(1) = 0; beta(end) = 1;
    
    

end

function y = g(x)
    y = 1/(1+exp(-x));
end