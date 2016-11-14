function Z = TAP3(W, a, b, Tlimit)
%TAP 此处显示有关此函数的摘要
%   a and b are row vectors
%     figure;
    mh = rand(1,size(W,2));
    mv = rand(1,size(W,1));
    flag = true;
    
    T = 1000;
    T_pre = 0;
    i = 1;
    while flag
        d = abs(T-T_pre);
        if d >Tlimit
            mh = logsig(b + mv*W - (mv-mv.^2)*W.^2.*(mh-0.5) + ...
                (2*(mh.^2-mh)+1/3).*(((0.5-mv).*(mv-mv.^2))*(W.^3)));
            mv = logsig(a + (W*mh')' - (W.^2*(mh-mh.^2)'.*(mv-0.5)')' + ...
                (2*(mv.^2-mv)+1/3).*(((0.5-mh).*(mh-mh.^2))*(W.^3)'));
                T_pre = T;
                S = - (mh*log(mh') + (1-mh)*log((1-mh)') + mv*log(mv') + (1-mv)*log((1-mv)')); 
                T = -S - a*mv' - b*mh' - mv*W*mh' - 0.5*(mv-mv.^2)*W.^2*(mh-mh.^2)' - ...
                    2/3*(((mv-mv.^2).*(0.5-mv))*W.^3*((mh-mh.^2).*(0.5-mh))');
                if isnan(T)
                    flag = false;
                    T = T_pre;
                end
        else
            flag = false;
        end
%         scatter(i,d);hold on;
        i = i + 1;
    end
   
    Z = -T;
end
